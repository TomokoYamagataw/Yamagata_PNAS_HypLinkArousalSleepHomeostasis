%% compare arousal delay high and low sleep pressure conditions
clear all
%% load emg data and spectra etc 
paths.extractedSignals ='F:\Tomoko_OptoStim\EMGAnalysis_matfiles';
sedation=load(fullfile(paths.extractedSignals,'sedation'));
tenHz=load(fullfile(paths.extractedSignals,'tenHzSameAsSedation'));


% % bit ugly but this makes sure we are looking only at mice that had 10Hz stim as well as sedation
% missingMouseIdx = 2; % the second mouse has not been in the sedation dataset
% for currField = (fieldnames(tenHz))'
%    if iscell(tenHz.(currField{1}))
%        tenHz.(currField{1})(missingMouseIdx)=[];
%    end    
% end
tenHz.stims = cellfun(@(x) x((x(:,1)./(256*3600))<12,:),tenHz.stims,'un',0); % remove dark period

% colormap 
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 240 228 66;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue,yellow
npgCMap = [npgCMap;npgCMap];

%% extract features: arousal delay, preceeding power in delta, theta, sigma, time to last waking
EMGThresholds= [1,2,3,4,6,8]'; % EMG threshold
nEpochsSWA= ([1])'; % EMG threshold
sedation.settings.EMGThresholds = EMGThresholds;
sedation.settings.nEpochsSWA = nEpochsSWA;
sedation.settings.sedation = 1; % to use "sedation" as a scoring output
sedation.settings.minDur = 4; % to use "sedation" as a scoring output
tenHz.settings.EMGThresholds = EMGThresholds;
tenHz.settings.nEpochsSWA = nEpochsSWA;
tenHz.settings.minDur = 4;

% detect arousals based on thresholds as before

try sedation.settings=rmfield(sedation.settings,'forceThreshold'); catch; end
try tenHz.settings=rmfield(tenHz.settings,'forceThreshold'); catch; end
opto_sedation=extract_features_from_opto_evoked_arousals(sedation.EMGs,sedation.stims,sedation.scorings,sedation.EMGVarRatios,sedation.EMGVarBinning,sedation.settings); % runs the script
opto_tenHz=extract_features_from_opto_evoked_arousals(tenHz.EMGs,tenHz.stims,tenHz.scorings,tenHz.EMGVarRatios,tenHz.EMGVarBinning,tenHz.settings); % runs the script


%% comapre arousal delays sedation vs NREM
% Note: we use all NREM stims in the 24h 10 Hz but only the first 4 stims in the sedation group
% However, the effect is weaker but the same (and also significant) when using all stims, the anesthetic just runs out after 2h
f1=figure('name','ignoring mouse ID');
f2=figure('name','using mean of each mouse');
% f3=figure('name','distribution of arousal delay in each mouse');
tArousNR_tenHz=opto_tenHz.tArousal_NR(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
tArousNR_sedation=opto_sedation.tArousal_Sed(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
pValsPairedT = nan(numel(EMGThresholds),1);
pValsLME = nan(numel(EMGThresholds),1);
nMice = size(tArousNR_tenHz,2);
EMGCntForFigs = 2; % which one was used for figures 

for EMGCnt = 1:numel(EMGThresholds) % cycle through various emg thresholds
    %extract the arousalT based on the current EMG threshold
    acrossMiceNR_tenHz = cell2mat(tArousNR_tenHz(EMGCnt,:));
    withinMiceNR_tenHz = cellfun(@(x)mean(x),tArousNR_tenHz(EMGCnt,:));
    acrossMiceNR_Sed = cell2mat(tArousNR_sedation(EMGCnt,:));
    withinMiceNR_Sed = cellfun(@(x)nanmean(x(1:4)),tArousNR_sedation(EMGCnt,:));

%     plot difference HSP LSP across all samples, ignoring mouse ID
    figure(f1)
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    combined = [acrossMiceNR_tenHz,acrossMiceNR_Sed]';
    groups = [ones(size(acrossMiceNR_tenHz)),2.*ones(size(acrossMiceNR_Sed))]';
    boxplot(combined,groups,'notch','on');
    h = findobj(gca,'Tag','Box');
    h2= findobj(gca,'Tag','Outliers');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),npgCMap(j,:),'FaceAlpha',.4);
        h(j).Color = 'k'; % set box edge color
        h2(j).MarkerEdgeColor=[0.8 0.8 0.8];
    end
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]); 
    ylabel('arousal delay [s]'); ax = gca; ax.XTickLabel ={'10Hz 12h NREM','sedation'};
    set(gca,'TickDir','out'); box off;set(gcf,'Position',[205 484 1264 415]);
    linkaxes

    if EMGCnt == EMGCntForFigs
        display(['mean sedation: ',num2str(mean(withinMiceNR_Sed))])
        display(['SEM  sedation: ',num2str(std(withinMiceNR_Sed)./sqrt(size(withinMiceNR_Sed,2)))])
        display(['mean  NREM: ',num2str(mean(withinMiceNR_tenHz))])
        display(['SEM NREM: ',num2str(std(withinMiceNR_tenHz)./sqrt(size(withinMiceNR_tenHz,2)))])
    end
    
    %do stats
    nStimsLSP=arrayfun(@(x) x.*ones(numel(tArousNR_tenHz{EMGCnt,x}),1),1:nMice,'un',0);
    nStimsHSP=arrayfun(@(x) x.*ones(numel(tArousNR_sedation{EMGCnt,x}),1),1:nMice,'un',0);
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
    combined = [withinMiceNR_tenHz;withinMiceNR_Sed];
    boxplot([withinMiceNR_tenHz,withinMiceNR_Sed],[ones(size(withinMiceNR_tenHz)),2.*ones(size(withinMiceNR_Sed))]);hold on;
    plot(1:2,[withinMiceNR_tenHz;withinMiceNR_Sed],'-ok','MarkerFaceColor','k','MarkerSize',4); hold on;
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
%     
%     
%     plot(1:2,combined,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3)
%     errorbar(1:2,nanmean(combined,2),std(combined,[],2)./sqrt(size(combined,2)),'ok','MarkerFaceColor','k')
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('arousal delay [s]'); ax = gca; ax.XTickLabel ={'NREM','sedation'};
    set(gca,'TickDir','out'); box off
%     ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415])
    [h,pValsPairedT(EMGCnt)]=ttest(withinMiceNR_tenHz,withinMiceNR_Sed,'tail','both');
    linkaxes
    text(1,max(get(gca,'YLim'))-1,['p= ',num2str(pValsPairedT(EMGCnt))])

%        % plot distributions of arousal delay in HSP and LSP 
%     figure(f3)
%     subplot(1,numel(EMGThresholds),EMGCnt); hold on;
%     plot(sort(acrossMiceNR_LSP),100.*(1:numel(acrossMiceNR_LSP))./numel(acrossMiceNR_LSP),'k'); hold on;
%     plot(sort(acrossMiceNR_HSP),100.*(1:numel(acrossMiceNR_HSP))./numel(acrossMiceNR_HSP),'Color',npgCMap(1,:));
%     legend('low sleep pressure','high sleep pressure')
%     xlabel('arousal delay');ylabel('% stimulations');  
%     linkaxes
end

saveas(f1,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','tArousal_allDataPoints'))
saveas(f1,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','tArousal_allDataPoints.png'))
saveas(f2,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','tArousal_1DataPointPerMouse'))
saveas(f2,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','tArousal_1DataPointPerMouse.png'))
%% plot average opto-EMG variance 
figure;
iteration = 1; %which threshold to pklot
conditions = {opto_tenHz,opto_sedation};
colors = {[0,0,0],npgCMap(1,:)};
for lowHigh = 1:2
    for mouseCnt = 1:nMice
        currCond = conditions{lowHigh};
        subplot(4,2,mouseCnt)
        timeInS = -tenHz.settings.tPre:0.25:tenHz.settings.tPost-0.25;
        if lowHigh == 1
            meanVar= median(currCond.optoEvokedEMG_NR{iteration,mouseCnt},1);
            VarPtiles = prctile(currCond.optoEvokedEMG_NR{iteration,mouseCnt},[25,75],1);
            meanTArous=mean(currCond.tArousal_NR{iteration,mouseCnt}); %mean time to crossing "arousal" threshold
        else
            meanVar= median(currCond.optoEvokedEMG_Sed{iteration,mouseCnt},1);
            VarPtiles = prctile(currCond.optoEvokedEMG_Sed{iteration,mouseCnt},[25,75],1);
            meanTArous=mean(currCond.tArousal_Sed{iteration,mouseCnt}); %mean time to crossing "arousal" threshold
        end
        
        currThresh = currCond.thresholdsUsed{iteration,mouseCnt};
        plot(timeInS,meanVar,'Color',colors{lowHigh}); hold on; % plot mean evoked variance
        patch([timeInS,fliplr(timeInS)],[VarPtiles(1,:),fliplr(VarPtiles(2,:))],'k','FaceAlpha',0.3,'EdgeColor','none','FaceColor',colors{lowHigh}); % add 25/75%percentiles as shading
        xlabel('time to stim [s]'); ylabel('log EMG var'); box off; set(gca,'TickDir','out','YScale','log');xlim([0,20]); %cosmetics
        if mouseCnt ==1 && lowHigh==1
            legend('median LSP','25/75 %tiles','AutoUpdate','off'); legend boxoff          
        end
        line([meanTArous,meanTArous],get(gca,'YLim'),'Color',colors{lowHigh})
        line(get(gca,'XLim'),[currThresh,currThresh],'Color',colors{lowHigh})
    end
end
linkaxes
saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','mean_evoked_EMG_Var'))
saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','mean_evoked_EMG_Var.png'))

%% rest is just copied from high vs low arousal, dont think I'll need to look at this 
% %% sanity check: plot swan at time of stim in the two conditions 
% sedation.settings.EMGThresholds = [1]'; % EMG threshold
% sedation.settings.nEpochsSWA =  ([1,2,4,8])';
% sedation.settings.forceThreshold = higherThresholds;
% tenHz.settings.EMGThresholds = [1]'; % EMG threshold
% tenHz.settings.nEpochsSWA =  ([1,2,4,8])';
% tenHz.settings.forceThreshold = higherThresholds;
% 
% opto_sedation=extract_features_from_opto_evoked_arousals(sedation.EMGs,sedation.stims,sedation.scorings,sedation.EMGVarRatios,sedation.EMGVarBinning,sedation.settings); % runs the script
% opto_tenHz=extract_features_from_opto_evoked_arousals(tenHz.EMGs,tenHz.stims,tenHz.scorings,tenHz.EMGVarRatios,tenHz.EMGVarBinning,tenHz.settings); % runs the script
% pValsPairedT=nan(numel(sedation.settings.nEpochsSWA ),1);
% f1=figure;
% f2=figure;
% for nSWA = 1:numel(settings.nEpochsSWA)
%     tStimLSP = opto_tenHz.epIdxNRStim(nSWA,:);
%     tStimHSP = opto_sedation.epIdxNRStim(nSWA,:);
%     currSWALSP = opto_tenHz.preStimSWA(nSWA,:);
%     currSWAHSP = opto_sedation.preStimSWA(nSWA,:);
%     currtArousLSP = opto_tenHz.tArousal_NR(nSWA,:);
%     currtArousHSP = opto_sedation.tArousal_NR(nSWA,:); 
%     figure(f1);
%     for mouseCnt = 1:numel(currSWALSP)
%         subplot(numel(settings.nEpochsSWA),numel(currSWALSP),mouseCnt + (nSWA-1)*numel(currSWALSP))
% %         subplot(1,numel(currSWALSP),mouseCnt)
% 
%         meanSWACurrMouse = nanmean([currSWALSP{mouseCnt};currSWAHSP{mouseCnt}]);
%         meanTArousCurrMouse = nanmean([currtArousLSP{mouseCnt},currtArousHSP{mouseCnt}]);
%         plot(tStimLSP{mouseCnt}.*4./3600,100.*currSWALSP{mouseCnt}./meanSWACurrMouse,'k-o','LineWidth',1.5,'MarkerFaceColor','k','MarkerSize',3); hold on;
% %         plot(currtArousLSP{mouseCnt}./meanTArousCurrMouse,'k--','LineWidth',1.5); hold on;
% %         plot(currtArousHSP{mouseCnt}./meanTArousCurrMouse,'Color',npgCMap(1,:),'LineWidth',1.5); 
%         plot(tStimHSP{mouseCnt}.*4./3600,100.*currSWAHSP{mouseCnt}./meanSWACurrMouse,'o-','Color',npgCMap(1,:),'LineWidth',1.5,'MarkerFaceColor',npgCMap(1,:),'MarkerSize',3); 
%         ylabel('SWA [% mean]'); xlabel('hours (ZT)'); box off; set(gca,'TickDir','out'); 
%     end
%     
%     %plot mean SWA high vs low for each mouse
%     figure(f2)
%     subplot(1,numel(settings.nEpochsSWA),nSWA)
%     meanByMouse_LSP = log10(cellfun(@nanmean,currSWALSP));
%     meanByMouse_HSP = log10(cellfun(@nanmean,currSWAHSP));
%     combined = [meanByMouse_LSP;meanByMouse_HSP];
%     plot(1:2,combined,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3); hold on;
%     errorbar(1:2,nanmean(combined,2),std(combined,[],2)./sqrt(size(combined,2)),'ok','MarkerFaceColor','k')
%     title(['n Epochs SWA: ',num2str(settings.nEpochsSWA(nSWA))]);
%     ylabel('SWA [log]'); ax = gca; ax.XTickLabel ={'10Hz 12h NREM','sedation'};
%     set(gca,'TickDir','out'); box off
%     ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415]); 
%     [h,pValsPairedT(nSWA)]=ttest(meanByMouse_LSP,meanByMouse_HSP,'tail','both');
%     linkaxes
%     text(1,max(get(gca,'YLim'))-0.1,['p= ',num2str(pValsPairedT(nSWA))])
% end
% saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','changeInSWA_sanityCheck'))
% saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','changeInSWA_sanityCheck.png'))
% %% compare tArousal at time of highest SWA in SD condition 
% sedation.settings.EMGThresholds = [1,2,3,4,6,8]'; % EMG threshold
% sedation.settings.nEpochsSWA =  ([4])';
% sedation.settings.forceThreshold = higherThresholds;
% tenHz.settings.EMGThresholds = [1,2,3,4,6,8]'; % EMG threshold
% tenHz.settings.nEpochsSWA =  ([4])';
% tenHz.settings.forceThreshold = higherThresholds;
% 
% 
% EMGThresholds = settings.EMGThresholds;
% opto_sedation=extract_features_from_opto_evoked_arousals(sedation.EMGs,sedation.stims,sedation.scorings,sedation.EMGVarRatios,sedation.EMGVarBinning,sedation.settings); % runs the script
% opto_tenHz=extract_features_from_opto_evoked_arousals(tenHz.EMGs,tenHz.stims,tenHz.scorings,tenHz.EMGVarRatios,tenHz.EMGVarBinning,tenHz.settings); % runs the script
% 
% f1=figure('name','tArousal max SWA, 1 val per mouse');
% f2=figure('name',' max SWA, 1 val per mouse');
% tArousNR_tenHz=opto_tenHz.tArousal_NR(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
% tArousNR_sedation=opto_sedation.tArousal_NR(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
% maxSWA=nan(4,nMice); % 1-2 tArousal LSP HSP, 3-4 SWA LSP HSP 
% pValsPairedT = nan(numel(EMGThresholds))
% for EMGCnt = 1:numel(EMGThresholds) % cycle through various emg thresholds
%     currTArousLSP = tArousNR_tenHz(EMGCnt,:);
%     currTArousHSP = tArousNR_sedation(EMGCnt,:);
%     currSWALSP = opto_tenHz.preStimSWA(EMGCnt,:);
%     currSWAHSP = opto_sedation.preStimSWA(EMGCnt,:);
%     [maxSWA(4,:),maxSWAIdxHSP] = cellfun(@max,currSWAHSP); % get the stim nr where the SWA was highest
%     maxSWA(2,:)= arrayfun(@(x) currTArousHSP{x}(maxSWAIdxHSP(x)),1:nMice);
%     tMaxStim=arrayfun(@(mCnt) opto_sedation.epIdxNRStim{EMGCnt,mCnt}(maxSWAIdxHSP(mCnt)),1:nMice); % find epoch nr of max SWA in HSP
%     [diff,maxSWAIdxLSP]=arrayfun(@(mCnt) min(abs(opto_tenHz.epIdxNRStim{EMGCnt,mCnt}-tMaxStim(mCnt))),1:nMice); % find closest NREM stim epoch in LSP
%     maxSWA(1,:)= arrayfun(@(x) currTArousLSP{x}(maxSWAIdxLSP(x)),1:nMice);
%     maxSWA(3,:)= arrayfun(@(x) currSWALSP{x}(maxSWAIdxLSP(x)),1:nMice);
% 
%     figure(f1)
%     subplot(1,numel(EMGThresholds),EMGCnt)
%     combined = [maxSWA(1,:);maxSWA(2,:)];
%     plot(1:2,combined,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3); hold on;
%     errorbar(1:2,nanmean(combined,2),std(combined,[],2)./sqrt(size(combined,2)),'ok','MarkerFaceColor','k')
%     title(['EMG thresh: ',num2str(settings.EMGThresholds(EMGCnt))]);
%     ylabel('tArousal'); ax = gca; ax.XTickLabel ={'10Hz 12h NREM','sedation'};
%     set(gca,'TickDir','out'); box off
%     ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415]); 
%     [h,pValsPairedT(EMGCnt)]=ttest(maxSWA(1,:),maxSWA(2,:),'tail','both');
%     linkaxes
%     text(1,max(get(gca,'YLim'))-0.1,['p= ',num2str(pValsPairedT(EMGCnt))])
% 
%     figure(f2)
%     subplot(1,numel(EMGThresholds),EMGCnt)
%     combined = [maxSWA(3,:);maxSWA(4,:)];
%     plot(1:2,combined,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3); hold on;
%     errorbar(1:2,nanmean(combined,2),std(combined,[],2)./sqrt(size(combined,2)),'ok','MarkerFaceColor','k')
%     title(['EMG thresh: ',num2str(settings.EMGThresholds(EMGCnt))]);
%     ylabel('SWA'); ax = gca; ax.XTickLabel ={'10Hz 12h NREM','sedation'};
%     set(gca,'TickDir','out'); box off
%     ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415]); ax.YScale = 'log';
%     [h,currP]=ttest(maxSWA(3,:),maxSWA(4,:),'tail','both');
%     linkaxes
%     text(1,max(get(gca,'YLim'))-0.1,['p= ',num2str(currP)])
%     
% end
% 
% saveas(f1,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','tArousal_maxSWAHSP'))
% saveas(f1,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','tArousal_maxSWAHSP.png'))
% saveas(f2,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','SWA_maxSWAHSP'))
% saveas(f2,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\sedationVsTenHz','SWA_maxSWAHSP.png'))