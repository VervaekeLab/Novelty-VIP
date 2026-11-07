function vipGLM_exampleCellTrace()
% VIPGLM_EXAMPLECELLTRACE  Self-contained viewer for a single VIP cell:
% observed dF/F vs the GLM's HELD-OUT (cross-validated test) prediction, in a
% scrollable ~30 s window, with velocity + reward locations below.
%
% The window stays open: use "Load session/cell" to pick another session and
% VIP cell without re-running. A low-pass filter slider smooths the OBSERVED
% trace only (the GLM prediction is left untouched). Filtering is for display.
%
% Design matrix + 10-fold/10-s-chunk ridge CV are copied verbatim from
% vipGLM_fitCell (familiar context, binary position, 100 bins, lambda=0.5,
% lags -3..3 s, 7 pts) so the held-out prediction matches the published GLM.

    win=30;   % seconds visible

    f=figure('Color','w','Position',[80 80 1120 620], ...
        'Name','VIP example cell: observed vs GLM (held-out)','NumberTitle','off');

    % ---- persistent state ----
    st=struct('t',[],'ym',[],'yp',[],'vel',[],'rewT',[],'Fs',30,'win',win, ...
        'sfn','','cellID',NaN,'R2',NaN);
    guidata(f,st);

    % ---- axes ----
    axTr=axes('Parent',f,'Position',[0.08 0.46 0.86 0.46],'Tag','axTr'); hold(axTr,'on');
    axB =axes('Parent',f,'Position',[0.08 0.20 0.86 0.18],'Tag','axB'); hold(axB,'on');

    % ---- controls ----
    uicontrol('Parent',f,'Style','pushbutton','Units','normalized', ...
        'Position',[0.08 0.93 0.16 0.05],'String','Load session/cell', ...
        'FontWeight','bold','Callback',@(~,~) loadData(f));
    uicontrol('Parent',f,'Style','text','Units','normalized', ...
        'Position',[0.55 0.935 0.18 0.04],'String','Low-pass smoothing:', ...
        'BackgroundColor','w','HorizontalAlignment','right');
    % filter slider: 0 = none, up to ~2 s smoothing window
    uicontrol('Parent',f,'Style','slider','Units','normalized','Tag','lpf', ...
        'Position',[0.74 0.935 0.20 0.04],'Min',0,'Max',2,'Value',0, ...
        'Callback',@(~,~) redraw(f));
    addlistener(findobj(f,'Tag','lpf'),'ContinuousValueChange',@(~,~) redraw(f));

    % scroll slider (created/updated after data load)
    uicontrol('Parent',f,'Style','slider','Units','normalized','Tag','scroll', ...
        'Position',[0.08 0.05 0.86 0.04],'Min',0,'Max',1,'Value',0,'Visible','off', ...
        'Callback',@(~,~) doScroll(f));
    addlistener(findobj(f,'Tag','scroll'),'ContinuousValueChange',@(~,~) doScroll(f));

    loadData(f);   % prompt immediately on open
end

% ===================== load a session + cell, fit, store =====================
function loadData(f)
    lambdaVal=0.5; lagMin=-3; lagMax=3; nLags=7; nPosRBF=100; posBasis='binary';
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData');
    if ~exist(dfSess,'dir'), dfSess=userpath; end

        uiwait(msgbox({'Select ONE SESSION DATA file.','','A single m8xxx session .mat file.'},'Session data (one file)','help'));
        [sfn,sp]=uigetfile('*.mat','ONE SESSION DATA file (m8xxx .mat)',dfSess);
    if isequal(sfn,0), return; end
    S=load(fullfile(sp,sfn));
    if ~isfield(S,'sData'), errordlg('No sData in that file.'); return; end
    sData=S.sData;
    if ~isfield(sData,'imdata')||~isfield(sData.imdata,'roiSignals') ...
            || numel(sData.imdata.roiSignals)<2 || ~isfield(sData.imdata.roiSignals(2),'dff')
        errordlg('No roiSignals(2).dff in that session.'); return;
    end
    dff_all=sData.imdata.roiSignals(2).dff; nCells=size(dff_all,1);

    labels=arrayfun(@(c) sprintf('Cell %d',c),1:nCells,'UniformOutput',false);
    [sel,ok]=listdlg('PromptString','Select a VIP cell:','SelectionMode','single', ...
        'ListString',labels,'ListSize',[200 400]);
    if ~ok, return; end
    cellID=sel;

    hw=waitbar(0.3,sprintf('Fitting cell %d (cross-validated)...',cellID));
    try
        dff_len=size(dff_all,2); time_len=length(sData.behavior.signals.timeDs);
        pos_len=length(sData.behavior.signals.corridorPositionDs);
        vel_len=length(sData.behavior.signals.velocityDs);
        N=min([dff_len,time_len,pos_len,vel_len]);
        time=sData.behavior.signals.timeDs(1:N); time=time(:);
        Fs=1/mean(diff(time),'omitnan');
        if isnan(Fs)||Fs<=0, closew(hw); errordlg('bad Fs'); return; end

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
        else
            closew(hw); errordlg('no reward data'); return;
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
        fam_id=ctx(1);
        mask=(cctx==fam_id); mask=mask(:);
        if sum(mask)<50, closew(hw); errordlg('too few familiar frames'); return; end

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

        num_bases=nPosRBF;
        if num_bases>0 && strcmpi(posBasis,'binary')
            pos_edges=linspace(min(pos),max(pos),num_bases+1);
            pb=discretize(pos_m,pos_edges);
            X_pos=zeros(Nm,num_bases);
            for b=1:num_bases, X_pos(:,b)=double(pb==b); end
        else
            X_pos=zeros(Nm,0);
        end

        Xm=[X_vel,X_acc,X_lick,X_rew,X_pos,X_pupil]; Xm(isnan(Xm))=0;
        ym=dff_z;

        yp=nan(N_mask,1);
        for fo=1:k_folds
            te=(cv_idx==fo); tr=(cv_idx~=fo & cv_idx>0);
            if sum(tr)==0||sum(te)==0, continue; end
            Xtr=[ones(sum(tr),1),Xm(tr,:)]; Xte=[ones(sum(te),1),Xm(te,:)];
            pen=lambdaVal*eye(size(Xtr,2)); pen(1,1)=0;
            b=(Xtr'*Xtr+pen)\(Xtr'*ym(tr));
            yp(te)=Xte*b;
        end
        vi=~isnan(yp)&cv_idx>0;
        sst=sum((ym(vi)-mean(ym(vi))).^2); sse=sum((ym(vi)-yp(vi)).^2);
        R2=1-(sse/max(sst,eps));
    catch ME
        closew(hw); errordlg(sprintf('Fit failed: %s',ME.message)); return;
    end
    closew(hw);

    idx=find(mask);
    st=guidata(f);
    st.t=(0:N_mask-1)'/Fs; st.ym=ym; st.yp=yp;
    st.vel=vel(idx); rew=valve(idx); st.rewT=st.t(rew>0);
    st.Fs=Fs; st.sfn=sfn; st.cellID=cellID; st.R2=R2;
    guidata(f,st);

    set(f,'Name',sprintf('%s  |  Cell %d  |  held-out R^2=%.3f',sfn,cellID,R2));

    % configure scroll slider for this session length
    sc=findobj(f,'Tag','scroll'); tmax=st.t(end);
    if tmax>st.win
        set(sc,'Min',0,'Max',tmax-st.win,'Value',0,'Visible','on', ...
            'SliderStep',[st.win/4/(tmax-st.win) st.win/(tmax-st.win)]);
    else
        set(sc,'Visible','off');
    end
    redraw(f);
end

% ===================== redraw (filter + view) =====================
function redraw(f)
    st=guidata(f); if isempty(st.t), return; end
    axTr=findobj(f,'Tag','axTr'); axB=findobj(f,'Tag','axB');

    % low-pass: slider value = smoothing window in seconds (0 = none)
    lpfSec=get(findobj(f,'Tag','lpf'),'Value');
    yObs=st.ym;
    if lpfSec>0.01
        w=max(1,round(lpfSec*st.Fs));
        yObs=movmean(yObs,w,'omitnan');
    end

    cla(axTr); hold(axTr,'on');
    h1=plot(axTr,st.t,yObs,'k-','LineWidth',1.0);
    h2=plot(axTr,st.t,st.yp,'-','Color',[0.85 0.33 0.10],'LineWidth',1.4);
    ylabel(axTr,'\DeltaF/F (z)'); box(axTr,'off'); set(axTr,'TickDir','out');
    legend(axTr,[h1 h2],{'Observed','GLM (held-out)'},'Box','off','Location','northeast');
    if lpfSec>0.01, ftxt=sprintf(' (LP %.2f s)',lpfSec); else, ftxt=''; end
    title(axTr,sprintf('Cell %d: observed%s vs cross-validated GLM prediction',st.cellID,ftxt));

    cla(axB); hold(axB,'on');
    plot(axB,st.t,st.vel,'-','Color',[0 0.45 0.74],'LineWidth',1.0);
    yl=[min(st.vel) max(st.vel)]; if ~all(isfinite(yl))||yl(1)==yl(2), yl=[0 1]; end
    for rt=st.rewT', plot(axB,[rt rt],yl,'--','Color',[0.30 0.30 0.30],'LineWidth',0.8); end
    ylim(axB,yl);
    ylabel(axB,'Velocity (cm/s)'); xlabel(axB,'Time (s)'); box(axB,'off'); set(axB,'TickDir','out');
    text(axB,0.99,0.92,'dashed = reward','Units','normalized','HorizontalAlignment','right', ...
        'FontSize',8,'Color',[0.3 0.3 0.3]);

    linkaxes([axTr axB],'x');
    doScroll(f);
end

% ===================== scroll =====================
function doScroll(f)
    st=guidata(f); if isempty(st.t), return; end
    axTr=findobj(f,'Tag','axTr'); tmax=st.t(end);
    sc=findobj(f,'Tag','scroll');
    if strcmp(get(sc,'Visible'),'on'), x0=get(sc,'Value'); else, x0=0; end
    xlim(axTr,[x0 min(x0+st.win,tmax)]);
end

% ---- helpers ----
function X_shifted = make_shifted_matrix(vec, shifts)
    N=length(vec); X_shifted=zeros(N,length(shifts));
    for i=1:length(shifts)
        s=shifts(i);
        if s==0, X_shifted(:,i)=vec;
        elseif s>0, X_shifted(:,i)=[zeros(s,1);vec(1:end-s)];
        else, sa=abs(s); X_shifted(:,i)=[vec(sa+1:end);zeros(sa,1)]; end
    end
end

function closew(h), if ishandle(h), close(h); end, end
