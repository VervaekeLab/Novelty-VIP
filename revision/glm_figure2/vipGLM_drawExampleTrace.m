function ok = vipGLM_drawExampleTrace(axTop, axBot, sData, cellID, tWin, lpfSec, FSlab)
% VIPGLM_DRAWEXAMPLETRACE  Draw observed vs held-out GLM dF/F (top axes) and
% velocity with reward markers (bottom axes) for one VIP cell, into the supplied
% axes. Non-interactive version of vipGLM_exampleCellTrace, for embedding in the
% Fig.1 strip (panel-k reserved space).
%
%   axTop, axBot : target axes (top = traces, bottom = velocity)
%   sData        : loaded session struct
%   cellID       : row index into roiSignals(2).dff
%   tWin         : [t0 t1] seconds (within the masked familiar stream) to show
%   lpfSec       : low-pass smoothing window (s) applied to the OBSERVED trace
%   FSlab        : label font size (optional)
%
% Returns ok=false if the cell could not be fit. Fit settings match the cohort
% (familiar context, binary position, 100 bins, lambda=0.5, lags -3..3 s, 7 pts).

    if nargin<7||isempty(FSlab), FSlab=10; end
    ok=false;
    lambdaVal=0.5; lagMin=-3; lagMax=3; nLags=7; nPosRBF=100;

    try
        dff_all=sData.imdata.roiSignals(2).dff;
        dff_len=size(dff_all,2); time_len=length(sData.behavior.signals.timeDs);
        pos_len=length(sData.behavior.signals.corridorPositionDs);
        vel_len=length(sData.behavior.signals.velocityDs);
        N=min([dff_len,time_len,pos_len,vel_len]);
        time=sData.behavior.signals.timeDs(1:N); time=time(:);
        Fs=1/mean(diff(time),'omitnan');
        if isnan(Fs)||Fs<=0, return; end

        dff =double(dff_all(cellID,1:N)'); dff=dff(:);
        pos =double(sData.behavior.signals.corridorPositionDs(1:N)); pos=pos(:);
        vel =double(sData.behavior.signals.velocityDs(1:N)); vel=vel(:);
        lick=double(sData.behavior.signals.lickIfreqDs(1:N)); lick=lick(:);
        if length(vel)>2, vel(1:2)=vel(3); end
        acc=gradient(vel)*Fs; acc=acc(:);

        valve=zeros(N,1);
        if isfield(sData.behavior,'signals')&&isfield(sData.behavior.signals,'rewardDs')
            vt=double(sData.behavior.signals.rewardDs);
            valve(1:min(N,length(vt)))=vt(1:min(N,length(vt)))>0;
        elseif isfield(sData,'daqdata')&&isfield(sData.daqdata,'frameIndex')&&isfield(sData.daqdata,'waterValve')
            vr=double(sData.daqdata.waterValve>0.5);
            ft=find(diff(sData.daqdata.frameIndex)==1)+1; lim=min(N,length(ft)-1);
            for i=1:lim, if any(vr(ft(i):ft(i+1)-1)), valve(i)=1; end, end
            if N<=length(ft)&&any(vr(ft(N):end)), valve(N)=1; end
        end
        valve=valve(:);

        pupil=nan(N,1);
        try
            if isfield(sData.behavior,'pupildata')&&isfield(sData.behavior.pupildata,'diameter')
                pd=double(sData.behavior.pupildata.diameter); pd(isnan(pd))=inf;
                pfi=find(diff(sData.behavior.signals.pupilCamFrameIndexCorr)==1)+1;
                ifi=find(diff(sData.daqdata.frameIndex)==1);
                pUs=nan(size(sData.behavior.signals.pupilCamFrameIndexCorr));
                vl=min(length(pfi),length(pd)); pUs(pfi(1:vl))=pd(1:vl);
                pUs=fillmissing(pUs,'nearest'); pUs(isinf(pUs))=nan;
                pr=double(pUs(ifi));
                if length(pr)>=N, pupil=pr(1:N); else, pupil=[pr(:);nan(N-length(pr),1)]; end
                pupil=pupil(:);
            end
        catch, pupil=nan(N,1); end
        if all(isnan(pupil)), pupil_glm=zeros(N,1);
        else, pupil_glm=fillmissing(pupil,'nearest'); pupil_glm(isnan(pupil_glm))=mean(pupil_glm,'omitnan'); end

        ctx_raw=mode(sData.behavior.trialMatrices.contextMatrix,2);
        ctx=round(fillmissing(ctx_raw,'nearest')); ctx=ctx(:);
        trialIdx=sData.behavior.signals.trialIndexesDs(1:N); trialIdx=trialIdx(:);
        vtr=trialIdx>0&trialIdx<=length(ctx);
        cctx=nan(N,1); cctx(vtr)=ctx(trialIdx(vtr)); cctx=fillmissing(cctx,'nearest');
        fam_id=ctx(1); mask=(cctx==fam_id); mask=mask(:);
        if sum(mask)<50, return; end

        N_mask=sum(mask); k_folds=10; chunk_size=round(Fs*10);
        num_chunks=ceil(N_mask/chunk_size);
        fa=repmat(1:k_folds,1,ceil(num_chunks/k_folds)); cv_idx=zeros(N_mask,1);
        for i=1:num_chunks
            cv_idx(((i-1)*chunk_size+1):min(i*chunk_size,N_mask))=fa(i);
        end

        zc=@(v)(v-mean(v,'omitnan'))./max(std(v,'omitnan'),1e-8);
        vel_m=vel(mask); acc_m=acc(mask); lick_m=lick(mask);
        valve_m=valve(mask); pos_m=pos(mask); pupil_m=pupil_glm(mask);
        dff_m=dff(mask); Nm=sum(mask);
        dff_z=zc(dff_m); dff_z(isnan(dff_z))=0;

        sigma_frames=round(0.125*Fs); kernel=gausswin(round(sigma_frames*3)); kernel=kernel/sum(kernel);
        shifts=round(linspace(lagMin,lagMax,nLags)*Fs);
        X_vel  =make_shifted_matrix(conv(zc(vel_m),  kernel,'same'),shifts);
        X_acc  =make_shifted_matrix(conv(zc(acc_m),  kernel,'same'),shifts);
        X_lick =make_shifted_matrix(conv(zc(lick_m), kernel,'same'),shifts);
        X_rew  =make_shifted_matrix(conv(zc(valve_m),kernel,'same'),shifts);
        X_pupil=make_shifted_matrix(conv(zc(pupil_m),kernel,'same'),shifts);
        pos_edges=linspace(min(pos),max(pos),nPosRBF+1);
        pb=discretize(pos_m,pos_edges); X_pos=zeros(Nm,nPosRBF);
        for b=1:nPosRBF, X_pos(:,b)=double(pb==b); end

        Xm=[X_vel,X_acc,X_lick,X_rew,X_pos,X_pupil]; Xm(isnan(Xm))=0; ym=dff_z;
        yp=nan(N_mask,1);
        for fo=1:k_folds
            te=(cv_idx==fo); tr=(cv_idx~=fo & cv_idx>0);
            if sum(tr)==0||sum(te)==0, continue; end
            Xtr=[ones(sum(tr),1),Xm(tr,:)]; Xte=[ones(sum(te),1),Xm(te,:)];
            pen=lambdaVal*eye(size(Xtr,2)); pen(1,1)=0;
            b=(Xtr'*Xtr+pen)\(Xtr'*ym(tr)); yp(te)=Xte*b;
        end
    catch
        return;
    end

    idx=find(mask); t=(0:N_mask-1)'/Fs;
    velP=vel(idx); rew=valve(idx); rewT=t(rew>0);

    % low-pass the observed trace only
    yObs=ym;
    if lpfSec>0.01, w=max(1,round(lpfSec*Fs)); yObs=movmean(yObs,w,'omitnan'); end

    if isempty(tWin), tWin=[t(1) min(t(1)+30,t(end))]; end
    sel=t>=tWin(1)&t<=tWin(2);
    LW=1.4;   % match the rest of the figure's axis line width

    % ---- top: observed vs predicted ----
    axes(axTop); cla(axTop); hold(axTop,'on');
    hObs =plot(axTop,t(sel),yObs(sel),'k-','LineWidth',0.9);
    hPred=plot(axTop,t(sel),yp(sel),'-','Color',[0.85 0.33 0.10],'LineWidth',1.2);
    xlim(axTop,tWin);
    set(axTop,'TickDir','out','LineWidth',LW,'FontSize',FSlab-2, ...
        'XTick',[],'XColor','none');     % no x-axis on top panel (velocity panel carries it)
    ylabel(axTop,'\DeltaF/F (z)','FontSize',FSlab-1); box(axTop,'off');
    % tight y-limits so the traces fill the panel (larger visual amplitude); only a
    % small pad. Legend is placed ABOVE the axes (northoutside) so it needs no
    % in-panel headroom and sits high, near the panel letter.
    dlo=min([yObs(sel);yp(sel)]); dhi=max([yObs(sel);yp(sel)]);
    pad=(dhi-dlo)*0.08; ylim(axTop,[dlo-pad dhi+pad]);
    legend(axTop,[hObs hPred],{'Observed','Predicted'},'Box','off', ...
        'FontSize',FSlab-3,'Orientation','horizontal','Location','northoutside');

    % ---- bottom: velocity + reward dashes ----
    axes(axBot); cla(axBot); hold(axBot,'on');
    plot(axBot,t(sel),velP(sel),'-','Color',[0 0.45 0.74],'LineWidth',0.9);
    yl=[min(velP(sel)) max(velP(sel))]; if ~all(isfinite(yl))||yl(1)==yl(2), yl=[0 1]; end
    rwin=rewT(rewT>=tWin(1)&rewT<=tWin(2));
    for rt=rwin', plot(axBot,[rt rt],yl,'--','Color',[0.30 0.30 0.30],'LineWidth',0.7); end
    ylim(axBot,yl); xlim(axBot,tWin);
    % coarse y-ticks (0 and 30 cm/s) so the labels don't overlap in the short panel
    vt=[0 30]; vt=vt(vt>=yl(1)&vt<=yl(2));
    if isempty(vt), vt=round([yl(1) yl(2)]); end
    set(axBot,'TickDir','out','LineWidth',LW,'FontSize',FSlab-2,'YTick',vt);
    box(axBot,'off');
    % single-line short label: the velocity panel is very short, so a 2-line label
    % overflows vertically into the dF/F label above. Keep it to one short token.
    ylabel(axBot,'cm/s','FontSize',FSlab-2);
    xlabel(axBot,'Time (s)','FontSize',FSlab-1);

    ok=true;
end

function X_shifted = make_shifted_matrix(vec, shifts)
    N=length(vec); X_shifted=zeros(N,length(shifts));
    for i=1:length(shifts)
        s=shifts(i);
        if s==0, X_shifted(:,i)=vec;
        elseif s>0, X_shifted(:,i)=[zeros(s,1);vec(1:end-s)];
        else, sa=abs(s); X_shifted(:,i)=[vec(sa+1:end);zeros(sa,1)]; end
    end
end
