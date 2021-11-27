%% compare arousal delay high and low sleep pressure conditions
clear all
%% load emg data and spectra etc 
paths.extractedSignals ='F:\Tomoko_OptoStim\EMGAnalysis_matfiles';
hsp=load(fullfile(paths.extractedSignals,'HSP'));
lsp=load(fullfile(paths.extractedSignals,'LSP'));


% colormap 
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 240 228 66;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue,yellow
npgCMap = [npgCMap;npgCMap];
HSPColor = [0 192 0]./255;
%% extract features: arousal delay, preceeding power in delta, theta, sigma, time to last waking
%% extract features: arousal delay, preceeding power in delta, theta, sigma, time to last waking
EMGThresholds= [1,2,3,4,6,8]'; % EMG threshold
nEpochsSWA= ([1])'; % EMG threshold
hsp.settings.EMGThresholds = EMGThresholds;
hsp.settings.nEpochsSWA = nEpochsSWA;
hsp.settings.minDur = 4;
lsp.settings.minDur = 4;
hsp.settings.tPre = 15;
lsp.settings.tPre = 15;

lsp.settings.EMGThresholds = EMGThresholds;
lsp.settings.nEpochsSWA = nEpochsSWA;

% detect arousals based on thresholds as before

try lsp.settings=rmfield(lsp.settings,'forceThreshold'); hsp.settings=rmfield(hsp.settings,'forceThreshold'); catch; end
opto_HSP=extract_features_from_opto_evoked_arousals(hsp.EMGs,hsp.stims,hsp.scorings,hsp.EMGVarRatios,hsp.EMGVarBinning,hsp.settings); % runs the script
opto_LSP=extract_features_from_opto_evoked_arousals(lsp.EMGs,lsp.stims,lsp.scorings,lsp.EMGVarRatios,lsp.EMGVarBinning,lsp.settings); % runs the script

% for each EMG threshold, get the higher value between hsp and lsp 
thresholdsHSP=opto_HSP.thresholdsUsed(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
thresholdsLSP=opto_LSP.thresholdsUsed(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
higherThresholds = cell(size(thresholdsHSP));
for i = 1:numel(thresholdsHSP)
    if thresholdsHSP{i}>thresholdsLSP{i}
        higherThresholds{i}=thresholdsHSP{i};
    else
        higherThresholds{i}=thresholdsLSP{i};
        
    end
end

%rerun threshold based detection using the same thresholds for both conditions
% detect arousals based on thresholds as before
hsp.settings.forceThreshold = higherThresholds;
lsp.settings.forceThreshold = higherThresholds;
opto_HSP=extract_features_from_opto_evoked_arousals(hsp.EMGs,hsp.stims,hsp.scorings,hsp.EMGVarRatios,hsp.EMGVarBinning,hsp.settings); % runs the script
opto_LSP=extract_features_from_opto_evoked_arousals(lsp.EMGs,lsp.stims,lsp.scorings,lsp.EMGVarRatios,lsp.EMGVarBinning,lsp.settings); % runs the script


%% comapre arousal delays hsp vs lsp 
f1=figure('name','ignoring mouse ID');
f2=figure('name','using mean of each mouse');
% f3=figure('name','distribution of arousal delay in each mouse');
tArousNR_lsp=opto_LSP.tArousal_NR(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
tArousNR_hsp=opto_HSP.tArousal_NR(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
pValsPairedT = nan(numel(EMGThresholds),1);
pValsLME = nan(numel(EMGThresholds),1);
nMice = size(tArousNR_lsp,2);

for EMGCnt = 1:numel(EMGThresholds) % cycle through various emg thresholds
    %extract the arousalT based on the current EMG threshold
    acrossMiceNR_LSP = cell2mat(tArousNR_lsp(EMGCnt,:));
    withinMiceNR_LSP = cellfun(@mean,tArousNR_lsp(EMGCnt,:));
    acrossMiceNR_HSP = cell2mat(tArousNR_hsp(EMGCnt,:));
    withinMiceNR_HSP = cellfun(@nanmean,tArousNR_hsp(EMGCnt,:));

    % plot difference HSP LSP across all samples, ignoring mouse ID
    figure(f1)
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    combined = [acrossMiceNR_LSP,acrossMiceNR_HSP]';
    groups = [ones(size(acrossMiceNR_LSP)),2.*ones(size(acrossMiceNR_HSP))]';
    boxplot(combined,groups,'notch','on');
    h = findobj(gca,'Tag','Box');
    h2= findobj(gca,'Tag','Outliers');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),npgCMap(j,:),'FaceAlpha',.4);
        h(j).Color = 'k'; % set box edge color
        h2(j).MarkerEdgeColor=[0.8 0.8 0.8];
    end
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]); 
    ylabel('arousal delay [s]'); ax = gca; ax.XTickLabel ={'LSP','HSP'};
    set(gca,'TickDir','out'); box off;set(gcf,'Position',[205 484 1264 415]);
    linkaxes

    %do stats
    nStimsLSP=arrayfun(@(x) x.*ones(numel(tArousNR_lsp{EMGCnt,x}),1),1:nMice,'un',0);
    nStimsHSP=arrayfun(@(x) x.*ones(numel(tArousNR_hsp{EMGCnt,x}),1),1:nMice,'un',0);
    mouseIDs = [cell2mat(nStimsLSP');cell2mat(nStimsHSP')];
    statsTable=table(combined,mouseIDs,groups);
    statsTable.mouseIDs=nominal(mouseIDs);
    statsTable.groups=nominal(groups);
    withGroupLME =fitlme(statsTable,'combined~1+groups+(1|mouseIDs)');
    LMEnull = fitlme(statsTable,'combined~1+(1|mouseIDs)');
    tmp=withGroupLME.anova;
    pValsLME(EMGCnt) = tmp.pValue(2); % get p value for group effect (same as llr test in this case)
    text(1,max(get(gca,'YLim'))-1,['p= ',num2str(pValsLME(EMGCnt))])

   % plot difference HSP LSP for each mouse sepparately
    figure(f2)
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    combined = [withinMiceNR_LSP;withinMiceNR_HSP];
    boxplot([withinMiceNR_LSP,withinMiceNR_HSP],[ones(size(withinMiceNR_LSP)),2.*ones(size(withinMiceNR_HSP))]);hold on;
    plot(1:2,[withinMiceNR_LSP;withinMiceNR_HSP],'-ok','MarkerFaceColor','k','MarkerSize',4); hold on;
%     plot(2,withinMiceNR_HSP,'-ok','MarkerFaceColor','k','MarkerSize',4);
    h = findobj(gca,'Tag','Box');
    h2= findobj(gca,'Tag','Outliers');
    h3= findobj(gca,'Tag','Median');

    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',.4,'EdgeColor','none');
        h(j).Color = 'none'; % set box edge color
        h2(j).MarkerEdgeColor='none'; h2(j).MarkerFaceColor='none';
        h3(j).LineWidth=1.2;
    end
    
    
%     
%     plot(1:2,combined,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3)
%     errorbar(1:2,nanmean(combined,2),std(combined,[],2)./sqrt(size(combined,2)),'ok','MarkerFaceColor','k')
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('arousal delay [s]'); ax = gca; ax.XTickLabel ={'LSP','HSP'};
    set(gca,'TickDir','out'); box off
    ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415])
    [h,pValsPairedT(EMGCnt)]=ttest(withinMiceNR_LSP,withinMiceNR_HSP,'tail','both');
    linkaxes
%     text(1,max(get(gca,'YLim'))-1,['p= ',num2str(pValsPairedT(EMGCnt))])

%        % plot distributions of arousal delay in HSP and LSP 
%     figure(f3)
%     subplot(1,numel(EMGThresholds),EMGCnt); hold on;
%     plot(sort(acrossMiceNR_LSP),100.*(1:numel(acrossMiceNR_LSP))./numel(acrossMiceNR_LSP),'k'); hold on;
%     plot(sort(acrossMiceNR_HSP),100.*(1:numel(acrossMiceNR_HSP))./numel(acrossMiceNR_HSP),'Color',npgCMap(1,:));
%     legend('low sleep pressure','high sleep pressure')
%     xlabel('arousal delay');ylabel('% stimulations');  
%     linkaxes
end

saveas(f1,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','tArousal_allDataPoints'))
saveas(f1,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','tArousal_allDataPoints.png'))
saveas(f2,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','tArousal_1DataPointPerMouse'))
saveas(f2,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','tArousal_1DataPointPerMouse.png'))
%% plot average opto-EMG variance 
figure;
iteration = 2; %which threshold to pklot
conditions = {opto_LSP,opto_HSP};
colors = {[0,0,0],HSPColor};
lHandles = {};
for lowHigh = 1:2
    for mouseCnt = 1:nMice
        currCond = conditions{lowHigh};
        subplot(4,2,mouseCnt)
        timeInS = -lsp.settings.tPre:0.25:lsp.settings.tPost-0.25;
        currEMGVar = currCond.optoEvokedEMG_NR{iteration,mouseCnt}; 
        currEMGVar = smoothdata(currEMGVar,2,'movmean',25);
        meanVar= median(currEMGVar,1);
        VarPtiles = prctile(currEMGVar,[25,75],1);
        meanTArous=mean(currCond.tArousal_NR{iteration,mouseCnt}); %mean time to crossing "arousal" threshold
        currThresh = currCond.thresholdsUsed{iteration,mouseCnt};
        lHandles{lowHigh}=plot(timeInS,meanVar,'Color',colors{lowHigh}); hold on; % plot mean evoked variance
        patch([timeInS,fliplr(timeInS)],[VarPtiles(1,:),fliplr(VarPtiles(2,:))],'k','FaceAlpha',0.3,'EdgeColor','none','FaceColor',colors{lowHigh}); % add 25/75%percentiles as shading
        xlabel('time to stim [s]'); ylabel('log EMG var'); box off; set(gca,'TickDir','out','YScale','log');xlim([-5,15]); %cosmetics
        
        if mouseCnt ==1 && lowHigh==2
%             pause
            legend([lHandles{1} lHandles{2}],{'median LSP','median HSP'},'AutoUpdate','off'); legend boxoff          
        end
%         line([meanTArous,meanTArous],get(gca,'YLim'),'Color',colors{lowHigh})
%         line(get(gca,'XLim'),[currThresh,currThresh],'Color',colors{lowHigh})
    end
end
% linkaxes
% saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','mean_evoked_EMG_Var'))
% saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','mean_evoked_EMG_Var.png'))
%% sanity check: plot swan at time of stim in the two conditions 
hsp.settings.EMGThresholds = [2]'; % EMG threshold
hsp.settings.nEpochsSWA =  ([1,2,4,8])';
hsp.settings.forceThreshold = higherThresholds;
lsp.settings.EMGThresholds = [2]'; % EMG threshold
lsp.settings.nEpochsSWA =  ([1,2,4,8])';
lsp.settings.forceThreshold = higherThresholds;

opto_HSP=extract_features_from_opto_evoked_arousals(hsp.EMGs,hsp.stims,hsp.scorings,hsp.EMGVarRatios,hsp.EMGVarBinning,hsp.settings); % runs the script
opto_LSP=extract_features_from_opto_evoked_arousals(lsp.EMGs,lsp.stims,lsp.scorings,lsp.EMGVarRatios,lsp.EMGVarBinning,lsp.settings); % runs the script
pValsPairedT=nan(numel(hsp.settings.nEpochsSWA ),1);
f1=figure;
f2=figure;
for nSWA = 1:numel(lsp.settings.nEpochsSWA)
    tStimLSP = opto_LSP.epIdxNRStim(nSWA,:);
    tStimHSP = opto_HSP.epIdxNRStim(nSWA,:);
    currSWALSP = opto_LSP.preStimSWA(nSWA,:);
    currSWAHSP = opto_HSP.preStimSWA(nSWA,:);
    currtArousLSP = opto_LSP.tArousal_NR(nSWA,:);
    currtArousHSP = opto_HSP.tArousal_NR(nSWA,:); 
    figure(f1);
    for mouseCnt = 1:numel(currSWALSP)
        subplot(numel(lsp.settings.nEpochsSWA),numel(currSWALSP),mouseCnt + (nSWA-1)*numel(currSWALSP))
%         subplot(1,numel(currSWALSP),mouseCnt)

        meanSWACurrMouse = nanmean([currSWALSP{mouseCnt};currSWAHSP{mouseCnt}]);
        meanTArousCurrMouse = nanmean([currtArousLSP{mouseCnt},currtArousHSP{mouseCnt}]);
        plot(tStimLSP{mouseCnt}.*4./3600,100.*currSWALSP{mouseCnt}./meanSWACurrMouse,'k-o','LineWidth',1.5,'MarkerFaceColor','k','MarkerSize',3); hold on;
%         plot(currtArousLSP{mouseCnt}./meanTArousCurrMouse,'k--','LineWidth',1.5); hold on;
%         plot(currtArousHSP{mouseCnt}./meanTArousCurrMouse,'Color',npgCMap(1,:),'LineWidth',1.5); 
        plot(tStimHSP{mouseCnt}.*4./3600,100.*currSWAHSP{mouseCnt}./meanSWACurrMouse,'o-','Color',npgCMap(1,:),'LineWidth',1.5,'MarkerFaceColor',npgCMap(1,:),'MarkerSize',3); 
        ylabel('SWA [% mean]'); xlabel('hours (ZT)'); box off; set(gca,'TickDir','out'); 
    end
    
    %plot mean SWA high vs low for each mouse
    figure(f2)
    subplot(1,numel(lsp.settings.nEpochsSWA),nSWA)
    meanByMouse_LSP = log10(cellfun(@nanmean,currSWALSP));
    meanByMouse_HSP = log10(cellfun(@nanmean,currSWAHSP));
    combined = [meanByMouse_LSP;meanByMouse_HSP];
    plot(1:2,combined,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3); hold on;
    errorbar(1:2,nanmean(combined,2),std(combined,[],2)./sqrt(size(combined,2)),'ok','MarkerFaceColor','k')
    title(['n Epochs SWA: ',num2str(lsp.settings.nEpochsSWA(nSWA))]);
    ylabel('SWA [log]'); ax = gca; ax.XTickLabel ={'LSP','HSP'};
    set(gca,'TickDir','out'); box off
    ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415]); 
    [h,pValsPairedT(nSWA)]=ttest(meanByMouse_LSP,meanByMouse_HSP,'tail','both');
    linkaxes
    text(1,max(get(gca,'YLim'))-0.1,['p= ',num2str(pValsPairedT(nSWA))])
end
saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','changeInSWA_sanityCheck'))
saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','changeInSWA_sanityCheck.png'))
%% compare tArousal at time of highest SWA in SD condition 
hsp.settings.EMGThresholds = [1,2,3,4,6,8]'; % EMG threshold
hsp.settings.nEpochsSWA =  ([4])';
hsp.settings.forceThreshold = higherThresholds;
lsp.settings.EMGThresholds = [1,2,3,4,6,8]'; % EMG threshold
lsp.settings.nEpochsSWA =  ([4])';
lsp.settings.forceThreshold = higherThresholds;


opto_HSP=extract_features_from_opto_evoked_arousals(hsp.EMGs,hsp.stims,hsp.scorings,hsp.EMGVarRatios,hsp.EMGVarBinning,hsp.settings); % runs the script
opto_LSP=extract_features_from_opto_evoked_arousals(lsp.EMGs,lsp.stims,lsp.scorings,lsp.EMGVarRatios,lsp.EMGVarBinning,lsp.settings); % runs the script

f1=figure('name','tArousal max SWA, 1 val per mouse');
f2=figure('name',' max SWA, 1 val per mouse');
tArousNR_lsp=opto_LSP.tArousal_NR(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
tArousNR_hsp=opto_HSP.tArousal_NR(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
maxSWA=nan(4,nMice); % 1-2 tArousal LSP HSP, 3-4 SWA LSP HSP 
pValsPairedT = nan(numel(EMGThresholds))
for EMGCnt = 1:numel(EMGThresholds) % cycle through various emg thresholds
    currTArousLSP = tArousNR_lsp(EMGCnt,:);
    currTArousHSP = tArousNR_hsp(EMGCnt,:);
    currSWALSP = opto_LSP.preStimSWA(EMGCnt,:);
    currSWAHSP = opto_HSP.preStimSWA(EMGCnt,:);
    [maxSWA(4,:),maxSWAIdxHSP] = cellfun(@max,currSWAHSP); % get the stim nr where the SWA was highest
    maxSWA(2,:)= arrayfun(@(x) currTArousHSP{x}(maxSWAIdxHSP(x)),1:nMice);
    tMaxStim=arrayfun(@(mCnt) opto_HSP.epIdxNRStim{EMGCnt,mCnt}(maxSWAIdxHSP(mCnt)),1:nMice); % find epoch nr of max SWA in HSP
    [diff,maxSWAIdxLSP]=arrayfun(@(mCnt) min(abs(opto_LSP.epIdxNRStim{EMGCnt,mCnt}-tMaxStim(mCnt))),1:nMice); % find closest NREM stim epoch in LSP
    maxSWA(1,:)= arrayfun(@(x) currTArousLSP{x}(maxSWAIdxLSP(x)),1:nMice);
    maxSWA(3,:)= arrayfun(@(x) currSWALSP{x}(maxSWAIdxLSP(x)),1:nMice);

    figure(f1)
    subplot(1,numel(EMGThresholds),EMGCnt)
    combined = [maxSWA(1,:);maxSWA(2,:)];
    plot(1:2,combined,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3); hold on;
    errorbar(1:2,nanmean(combined,2),std(combined,[],2)./sqrt(size(combined,2)),'ok','MarkerFaceColor','k')
    title(['EMG thresh: ',num2str(lsp.settings.EMGThresholds(EMGCnt))]);
    ylabel('tArousal'); ax = gca; ax.XTickLabel ={'LSP','HSP'};
    set(gca,'TickDir','out'); box off
    ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415]); 
    [h,pValsPairedT(EMGCnt)]=ttest(maxSWA(1,:),maxSWA(2,:),'tail','both');
    linkaxes
    text(1,max(get(gca,'YLim'))-0.1,['p= ',num2str(pValsPairedT(EMGCnt))])

    figure(f2)
    subplot(1,numel(EMGThresholds),EMGCnt)
    combined = [maxSWA(3,:);maxSWA(4,:)];
    plot(1:2,combined,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3); hold on;
    errorbar(1:2,nanmean(combined,2),std(combined,[],2)./sqrt(size(combined,2)),'ok','MarkerFaceColor','k')
    title(['EMG thresh: ',num2str(lsp.settings.EMGThresholds(EMGCnt))]);
    ylabel('SWA'); ax = gca; ax.XTickLabel ={'LSP','HSP'};
    set(gca,'TickDir','out'); box off
    ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415]); ax.YScale = 'log';
    [h,currP]=ttest(maxSWA(3,:),maxSWA(4,:),'tail','both');
    linkaxes
    text(1,max(get(gca,'YLim'))-0.1,['p= ',num2str(currP)])
    
end

saveas(f1,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','tArousal_maxSWAHSP'))
saveas(f1,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','tArousal_maxSWAHSP.png'))
saveas(f2,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','SWA_maxSWAHSP'))
saveas(f2,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','SWA_maxSWAHSP.png'))