%% compare arousal delay high and low sleep pressure conditions
clear all
% V3: include 20Hz and GFP, which don't have the same mice nr as the other conditions - check stuff with LMEs
% V 5: replace bar charts by box plots
% change stats to RM ANOVA 
%% load emg data and spectra etc
paths.extractedSignals ='F:\Tomoko_OptoStim\EMGAnalysis_matfiles';
fiveHz=load(fullfile(paths.extractedSignals,'fiveHz24h'));
oneHz=load(fullfile(paths.extractedSignals,'oneHz24h'));
tenHz=load(fullfile(paths.extractedSignals,'TenHz24hSameMiceAs5Hz'));
twoHz=load(fullfile(paths.extractedSignals,'twoHz24h'));
baseline=load(fullfile(paths.extractedSignals,'LPOBl'));
twentyHz=load(fullfile(paths.extractedSignals,'twentyHz24h'));



pathForFigs ='F:\Tomoko_OptoStim\figures_EMGVar\oneVs2Vs5Vs10Hz';
% colormap
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 240 228 66;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue,yellow
npgCMap = [npgCMap;npgCMap];
nMice = numel(twoHz.EMGs);
addpath(genpath('C:\Users\Martin Kahn\Desktop\Scripts_and_Code\Matlab\'))

%% extract features: arousal delay, preceeding power in delta, theta, sigma, time to last waking
EMGThresholds= [1,2,3,4,6]'; % EMG threshold
nEpochsSWA= ([1])'; % EMG threshold
conditions = {baseline,oneHz,twoHz,fiveHz,tenHz,twentyHz};
opto = cell(numel(conditions),1);
thresholds = nan(numel(conditions),numel(EMGThresholds),nMice);
for condCnt = 1:numel(conditions)
    currCond = conditions{condCnt};
    currCond.settings.EMGThresholds = EMGThresholds;
    currCond.settings.nEpochsSWA = nEpochsSWA; 
     currCond.settings.minDur = 4;
    % detect arousals based on thresholds as before
    
    try settings=rmfield(settings,'forceThreshold'); catch; end
    opto{condCnt}=extract_features_from_opto_evoked_arousals(currCond.EMGs,currCond.stims,currCond.scorings,currCond.EMGVarRatios,currCond.EMGVarBinning,currCond.settings); % runs the script
end
%% now redo this all with same threshold across conditions
% optoSameThresh = cell(numel(conditions),1);
%     highestThresholds = squeeze(max(thresholds,[],1)); % get highest threshold for each setting and mouse across the conditions
% for condCnt = 1:numel(conditions)
%     % for each EMG threshold, get the higher value between hsp and lsp
%     %rerun threshold based detection using the same thresholds for both conditions
%     % detect arousals based on thresholds as before
%     currCond = conditions{condCnt};
%     currCond.settings.minDur = 4;
% %     currCond.settings.forceThreshold =highestThresholds ;
%     currCond.settings.EMGThresholds = EMGThresholds;
%     currCond.settings.nEpochsSWA = nEpochsSWA; 
%     optoSameThresh{condCnt}=extract_features_from_opto_evoked_arousals(currCond.EMGs,currCond.stims,currCond.scorings,currCond.EMGVarRatios,currCond.EMGVarBinning,currCond.settings); % runs the script
% end

%% comapre arousal delays
f1=figure('name','ignoring mouse ID');
f2=figure('name','using mean of each mouse');
f3=figure('name','percentage arousals after stim, using mean of each mouse');
f4=figure('name','variability in arousal delay, using mean of each mouse');
f5=figure('name','ratio arousal delay NREM/REM');
f6= figure('name','overview figure with single threshold');
XTickLabelString = {'no stim','1Hz','2Hz','5Hz','10Hz','20Hz'};
% f3=figure('name','distribution of arousal delay in each mouse');
pValsPairedT = nan(numel(EMGThresholds),1);
pValsLME = nan(numel(EMGThresholds),1);
EMGCntForFigure = 2;
nConds = numel(conditions); 
sPlots =[];
for EMGCnt = 1:numel(EMGThresholds) % cycle through various emg thresholds
    %extract the arousalT based on the current EMG threshold
    [tArousWithinMiceNR,nArousalsWithinMiceNR,varTArousalWithinMiceNR,tArousWithinMiceR,nArousalsWithinMiceR] = deal(cell(1,nConds));
    tArousAcrossMiceNR = cell(4,1);
    for condCnt = 1:numel(conditions)
        currTArousal =opto{condCnt}.tArousal_NR;
        currTArousalR =opto{condCnt}.tArousal_R;
        tArousAcrossMiceNR{condCnt} = cell2mat(currTArousal(EMGCnt,:));
        tArousWithinMiceNR{condCnt}= cellfun(@nanmean,currTArousal(EMGCnt,:));
        tArousWithinMiceR{condCnt}= cellfun(@nanmean,currTArousalR(EMGCnt,:));

        nArousalsWithinMiceNR{condCnt}=cellfun(@(x) 100.*sum(isnan(x)==0)./numel(x),currTArousal(EMGCnt,:)); % percent stims crossing arousal limit within 120s
        nArousalsWithinMiceR{condCnt}=cellfun(@(x) 100.*sum(isnan(x)==0)./numel(x),currTArousalR(EMGCnt,:)); % percent stims crossing arousal limit within 120s
        varTArousalWithinMiceNR{condCnt}=cellfun(@nanstd,currTArousal(EMGCnt,:));
    end
    % plot difference HSP LSP across all samples, ignoring mouse ID
    figure(f1)
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    combined = cat(2,tArousAcrossMiceNR{:});
    groups = arrayfun(@(x)ones(numel(tArousAcrossMiceNR{x}),1).*x,1:numel(conditions),'un',0);
    groups =cat(1,groups{:}); 
    boxplot(combined,groups,'notch','on');
    h = findobj(gca,'Tag','Box');
    h2= findobj(gca,'Tag','Outliers');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),npgCMap(j,:),'FaceAlpha',.4);
        h(j).Color = 'k'; % set box edge color
        h2(j).MarkerEdgeColor=[0.8 0.8 0.8];
    end
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('time to awakening [s]'); ax = gca; ax.XTickLabel =XTickLabelString;
    set(gca,'TickDir','out'); box off;set(gcf,'Position',[205 484 1264 415]);
    linkaxes
    sPlots1=gca;
    
    % plot difference HSP LSP for each mouse sepparately
    figure(f2)
    sPlots2=subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    combined = tArousWithinMiceNR;
    combinedGroups = arrayfun(@(x) x.*ones(numel(combined{x}),1),1:numel(combined),'un',0);  % mean per condition across all mice
    boxplot(cat(2,combined{:}),cat(1,combinedGroups{:})); hold on;
    h = findobj(gca,'Tag','Box');
    m =  findobj(gca,'Tag','Median');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',.4);
        h(j).Color = 'k'; % set box edge color
        h(j).MarkerEdgeColor=[0.8 0.8 0.8];
        m(j).LineWidth = 1.7;
    end
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('time to awakening [s]'); ax = gca; ax.XTickLabel =XTickLabelString;
    set(gca,'TickDir','out'); box off
    ax.XTick = 1:nConds; ax.XLim =  [0,nConds+1]; set(gcf,'Position',[205 484 1264 415])
    linkaxes
    sPlots2=gca;
    
    stats = cat(1,combined{1:5})';
    combos =  nchoosek(1:5,2);
    [posthoc,~]=arrayfun(@(x)ttest(stats(:,combos(x,1)),stats(:,combos(x,2))),1:size(combos,1),'un',0);
    
    figure(f3)
    sPlots3=subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    combined = nArousalsWithinMiceNR;
     combinedGroups = arrayfun(@(x) x.*ones(numel(combined{x}),1),1:numel(combined),'un',0);  % mean per condition across all mice
    boxplot(cat(2,combined{:}),cat(1,combinedGroups{:})); hold on;
    h = findobj(gca,'Tag','Box');
    m =  findobj(gca,'Tag','Median');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',.4);
        h(j).Color = 'k'; % set box edge color
        h(j).MarkerEdgeColor=[0.8 0.8 0.8];
        m(j).LineWidth = 1.7;
    end
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('trials with arousal [% total]'); ax = gca; ax.XTickLabel =XTickLabelString;
    set(gca,'TickDir','out'); box off
    ax.XTick = 1:nConds; ax.XLim = [0,nConds+1]; set(gcf,'Position',[205 484 1264 415])
    linkaxes;sPlots3=gca;
   
    stats = cat(1,combined{1:5})';
    combos =  nchoosek(1:5,2);
    [posthoc,~]=arrayfun(@(x)ttest(stats(:,combos(x,1)),stats(:,combos(x,2))),1:size(combos,1),'un',0);
    sidaksAlpha = 1- (1-0.05).^(1/10)
%     statsStruct = struct();
%     statsStruct.delay = cat(2,combined{:})';
%     conditions  =cellfun(@numel,combined);
%     mouseID = arrayfun(@(x) 1:conditions(x),1:nConds,'un',0);
%     mouseID{6}= mouseID{6}+3;
%     conditions=arrayfun(@(x) x.*ones(conditions(x),1),1:nConds,'un',0);
%     statsStruct.conditionID = cat(1,conditions{:});
%     statsStruct.mouseID = cat(2,mouseID{:})';
%     statsStruct = struct2table(statsStruct);
%     statsStruct(statsStruct.conditionID>4,:)=[];
%     statsStruct.mouseID = nominal(statsStruct.mouseID);
%     statsStruct.conditionID = nominal(statsStruct.conditionID);
%     currLME = fitlme(statsStruct,'delay~1+conditionID+(1|mouseID)','Fitmethod','ml');
%     currNullLME = fitlme(statsStruct,'delay~1+(1|mouseID)','Fitmethod','ml');
%     compare(currNullLME,currLME)
%     checkLMEAssumptions(currLME,statsStruct)
    
    
    figure(f4)
    sPlots4=subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    combined = varTArousalWithinMiceNR;
    combinedMean = cellfun(@mean,combined);  % mean per condition across all mice 
    combinedSEM = cellfun(@std,combined)./sqrt(cellfun(@numel,combined)); 
    bar(1:numel(conditions),combinedMean,'k','FaceAlpha',0.3,'EdgeColor','none'); hold on; 
    errorbar(1:nConds,combinedMean,combinedSEM,'ok','MarkerFaceColor','none','MarkerEdgeColor','none')
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('std of arousal delay [s]'); ax = gca; ax.XTickLabel =XTickLabelString;
    set(gca,'TickDir','out'); box off
    ax.XTick = 1:nConds; ax.XLim =  [0,nConds+1]; set(gcf,'Position',[205 484 1264 415])
    linkaxes;sPlots4=gca;
    
    
    figure(f5)
    sPlots5=subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    combined = arrayfun(@(x)100.*tArousWithinMiceR{x}./tArousWithinMiceNR{x},1:nConds,'un',0);
    combinedMean = cellfun(@mean,combined);  % mean per condition across all mice 
    combinedSEM = cellfun(@std,combined)./sqrt(cellfun(@numel,combined)); 
    errorbar(1:nConds,combinedMean,combinedSEM,'ok','MarkerFaceColor','k')
    
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('ratio NREM/REM [s]'); ax = gca; ax.XTickLabel =XTickLabelString(1:4);
    set(gca,'TickDir','out'); box off
    ax.XTick = 1:4; ax.XLim =  [0,4+1]; set(gcf,'Position',[205 484 1264 415])
    linkaxes;sPlots5=gca;
    

    if EMGCnt == EMGCntForFigure
      figure(f6); 
      new=subplot(2,3,1); new.Visible='off';s=copyobj(sPlots1,f6);s.Position =new.Position; s.Title.String='arousal delays across all trials and mice';
      new=subplot(2,3,2); new.Visible='off';s=copyobj(sPlots2,f6);s.Position =new.Position; s.Title.String='arousal delay, mean across mice';
      new=subplot(2,3,3); new.Visible='off';s=copyobj(sPlots3,f6);s.Position =new.Position; s.Title.String='percentage of trials with arousal during stim';
      new=subplot(2,3,4); new.Visible='off';s=copyobj(sPlots4,f6);s.Position =new.Position; s.Title.String='variability of arousal delay';
      new=subplot(2,3,5); new.Visible='off';s=copyobj(sPlots5,f6);s.Position =new.Position; s.Title.String='ratio NREM/REM';

    end
    
end
%%
saveas(f1,fullfile(pathForFigs,'tArousal_allDataPoints'))
saveas(f1,fullfile(pathForFigs,'tArousal_allDataPoints.png'))
saveas(f2,fullfile(pathForFigs,'tArousal_1DataPointPerMouse'))
saveas(f2,fullfile(pathForFigs,'tArousal_1DataPointPerMouse.png'))
saveas(f3,fullfile(pathForFigs,'nArousals_1DataPointPerMouse'))
saveas(f3,fullfile(pathForFigs,'nArousals_1DataPointPerMouse.png'))
saveas(f4,fullfile(pathForFigs,'stdTArousal_1DataPPerMouse'))
saveas(f4,fullfile(pathForFigs,'stdTArousal_1DataPPerMouse.png'))
saveas(f6,fullfile(pathForFigs,'overview'))
saveas(f6,fullfile(pathForFigs,'overview.png'))


%% plot average opto-EMG variance
figure;
iteration =2; %which threshold to pklot
colors = {[0,0,0],npgCMap(1,:),npgCMap(2,:),npgCMap(3,:),npgCMap(4,:),npgCMap(5,:)};
for condCnt = 1:nConds
    for mouseCnt = 1:nMice
        currCond =optoSameThresh{condCnt};
        s=subplot(4,2,mouseCnt);
        timeInS = -oneHz.settings.tPre:0.25:oneHz.settings.tPost-0.25;
        meanVar= median(currCond.optoEvokedEMG_NR{iteration,mouseCnt},1);
        VarPtiles = prctile(currCond.optoEvokedEMG_NR{iteration,mouseCnt},[25,75],1);
        meanTArous=nanmean(currCond.tArousal_NR{iteration,mouseCnt}); %mean time to crossing "arousal" threshold
        currThresh = currCond.thresholdsUsed{iteration,mouseCnt};
        p1=plot(timeInS,meanVar,'Color',colors{condCnt}); hold on; % plot mean evoked variance
        patch([timeInS,fliplr(timeInS)],[VarPtiles(1,:),fliplr(VarPtiles(2,:))],'k','FaceAlpha',0.3,'EdgeColor','none','FaceColor',colors{condCnt}); % add 25/75%percentiles as shading
        xlabel('time to stim [s]'); ylabel('log EMG var'); box off; set(gca,'TickDir','out','YScale','log');%xlim([0,20]); %cosmetics
        if mouseCnt ==1 && condCnt==4
            h = findobj(s,'Type','patch');
            legend(h,{'LSP','HSP','10Hz','5Hz','2Hz','1Hz'},'AutoUpdate','off'); legend boxoff
        end
        line([meanTArous,meanTArous],get(gca,'YLim'),'Color',colors{condCnt})
        line(get(gca,'XLim'),[currThresh,currThresh],'Color',colors{condCnt})
        title(['mouse Nr',num2str(mouseCnt)])
    end
end
linkaxes
% saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','mean_evoked_EMG_Var'))
% saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\HSPvsLSP','mean_evoked_EMG_Var.png'))
%% do correlation between SWA and tArousal for each condition and mouse 
%% extract features: arousal delay, preceeding power in delta, theta, sigma, time to last waking
freqs = 0:0.25:30;
EMGThresholds= [2]'; % EMG threshold
nEpochsSWA= ([1,2,4,6,8])'; % EMG threshold
conditions = {oneHz,twoHz,fiveHz,tenHz,baseline,twentyHz};
conditionNames = {'1Hz','2Hz','5Hz','10Hz','HSP','LSP'}; 
opto = cell(numel(conditions),1);
for condCnt = 1:nConds
    currCond = conditions{condCnt};
    currCond.settings.EMGThresholds = EMGThresholds;
    currCond.settings.nEpochsSWA = nEpochsSWA; 
    currCond.settings.minDur = 4;
    % detect arousals based on thresholds as before
    
    try settings=rmfield(settings,'forceThreshold'); catch; end
    opto{condCnt}=extract_features_from_opto_evoked_arousals(currCond.EMGs,currCond.stims,currCond.scorings,currCond.EMGVarRatios,currCond.EMGVarBinning,currCond.settings); % runs the script
end

%% now do correlations on looking 4 epochs back 
nEpochs = 5; 
figure;
% corrC = nan(4,nMice,121);
corrC = nan(nConds,nMice,59);
% SWARange = 3:9;
newFreqs = freqs;
SWARange = find(and(newFreqs>=0.5,newFreqs<=4));
spindleRange = find(and(newFreqs>=10,newFreqs<=15));
thetaRange = find(and(newFreqs>=5,newFreqs<=8));
[diffSWA,diffSpindles,diffTheta,diffTLastW,diffNRDur] = deal(nan(4,nMice));
sPlots = cell(3,nConds);
for condCnt = 1:nConds
    currOpto = opto{condCnt}; 
    currSpectra = currOpto.prestimSpectra(nEpochs,:);
    currMeanNrSpectra = currOpto.nrSpectra;
    currTArousal = currOpto.tArousal_NR(nEpochs,:);
    currTLastW = currOpto.tSinceWake(nEpochs,:);
    currNRDur = currOpto.NRBoutDur(nEpochs,:);

    nanEntries = arrayfun(@(x)isnan(currSpectra{x}(:,1)),1:nMice,'un',0);% some spectra will be nan as there weere insufficent preceeding epochs
%     currSpectra = arrayfun(@(x)currSpectra{x}(nanEntries{x}==0,:)./currMeanNrSpectra{nEpochs,x},1:nMice,'un',0);
    currSpectra = arrayfun(@(x)currSpectra{x}(nanEntries{x}==0,:),1:nMice,'un',0);
%     currSpectra = cellfun(@(x) mean(reshape(x(:,3:end-1),size(x,1),59,2),3),currSpectra,'un',0);
    nanArousals = cellfun(@isnan,currTArousal,'un',0);
    
%     for mouseCnt = 1:nMice
%         currTArousal{mouseCnt}(nanArousals{mouseCnt})=120;
%     end 
    currTArousal = arrayfun(@(x)currTArousal{x}(and(nanEntries{x}'==0,nanArousals{x}==0)),1:nMice,'un',0); % get rif of arousals without spectra
    currSpectra = arrayfun(@(x)currSpectra{x}(and(nanEntries{x}'==0,nanArousals{x}==0),:),1:nMice,'un',0); % get rif of arousals without spectra
    currTLastW = arrayfun(@(x)currTLastW{x}(and(nanEntries{x}'==0,nanArousals{x}==0),:),1:nMice,'un',0); % get rif of arousals without spectra
    currNRDur= arrayfun(@(x)currNRDur{x}(and(nanEntries{x}'==0,nanArousals{x}==0),:),1:nMice,'un',0); % get r
    
    % correlate spectra with tArousal
    for iFreq = 1:size(currSpectra{1},2)
        tmp=arrayfun(@(x)corrcoef(currTArousal{x},currSpectra{x}(:,iFreq)'),1:nMice,'un',0);
        corrC(condCnt,:,iFreq)= cellfun(@(x)x(2),tmp);
    end
    
    % compare and plot spectra of biggest with smallest 25% of arousal times 
    [~,I]=cellfun(@sort,currTArousal,'un',0);
    smallest20percent = cellfun(@(x)x(1:round(numel(x)/5)),I,'un',0);
    biggest20percent = cellfun(@(x)x(round(3*numel(x)/5):end),I,'un',0);
    tLastWSmallest  = arrayfun(@(x)nanmean(currTLastW{x}(smallest20percent{x})),1:nMice);
    tLastWLargest  = arrayfun(@(x)nanmean(currTLastW{x}(biggest20percent{x})),1:nMice);
    NRDurSmallest  = arrayfun(@(x)nanmean(currNRDur{x}(smallest20percent{x})),1:nMice);
    NRDurWLargest  = arrayfun(@(x)nanmean(currNRDur{x}(biggest20percent{x})),1:nMice);    
%     tArousSmallest  = arrayfun(@(x)nanmean(currTLastW{x}(smallest20percent{x})),1:nMice);
%     tArousLargest  = arrayfun(@(x)nanmean(currTLastW{x}(biggest20percent{x})),1:nMice);
    spectraSmallest  = cell2mat(arrayfun(@(x)nanmean(currSpectra{x}(smallest20percent{x},:),1),1:nMice,'un',0)');
    spectraLargest  = cell2mat(arrayfun(@(x)nanmean(currSpectra{x}(biggest20percent{x},:),1),1:nMice,'un',0)');
    spectraDiff = 100.*spectraSmallest./spectraLargest;
    diffSWA(condCnt,:) = mean(spectraDiff(:,SWARange),2);
    diffSpindles(condCnt,:) = mean(spectraDiff(:,spindleRange),2);
    diffTheta(condCnt,:) = mean(spectraDiff(:,thetaRange),2);
    diffTLastW(condCnt,:) = 100.*tLastWSmallest./tLastWLargest;
    diffNRDur(condCnt,:) = 100.*NRDurSmallest./NRDurWLargest;
    
    sPlots{1,condCnt}=subplot(3,nConds,condCnt);
    p1=patchMeUp(newFreqs,mean(spectraLargest,1),std(spectraLargest,[],1)./sqrt(nMice),'k',0.3); hold on;
    p2=patchMeUp(newFreqs,mean(spectraSmallest,1),std(spectraSmallest,[],1)./sqrt(nMice),npgCMap(1,:),0.3); hold on;
    ax = gca; ax.YScale = 'log';
    xlim([0,20]); title(conditionNames{condCnt}); legend([p1,p2],'25% longest','25% shortest'); legend('boxoff')
    ylabel('power [uV^2]');xlabel('frequency [Hz]')
    
    sPlots{2,condCnt}=subplot(3,nConds,condCnt+nConds);
    SEM = (std(spectraDiff,[],1)./sqrt(nMice));
    patchMeUp(newFreqs,mean(spectraDiff,1),SEM,'k',0.3); hold on;
    xlim([0,20]);  ylabel('difference in power [uV^2]');xlabel('frequency [Hz]')

    
    sPlots{3,condCnt}=subplot(3,nConds,condCnt+nConds*2);
    SEM = squeeze((std(corrC(condCnt,:,:),[],2)./sqrt(nMice)));
    patchMeUp(newFreqs,squeeze(mean(corrC(condCnt,:,:),2)),SEM,'k',0.3); hold on;
    ylabel('perasons R'); xlabel('frequency [Hz]') ;xlim([0,20]); 

end
linkaxes(cat(1,sPlots{1,:}));linkaxes(cat(1,sPlots{2,:}));linkaxes(cat(1,sPlots{3,:}));

% make boxplots of binned frequencies 
figure; 
groups = {'1Hz','2Hz','5Hz','10Hz','HSP','LSP'};
subplot(5,1,1);
boxplot(diffSWA',groups); hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),npgCMap(j,:),'FaceAlpha',.4);
    h(j).Color = 'k'; % set box edge color
    h(j).MarkerEdgeColor=[0.8 0.8 0.8];
end
line([0 nConds+1],[100 100],'LineStyle','--','Color','k'); box off; set(gca,'TickDir','out');
title('SWA'); ylabel('power ratio smallest/largest [%]');

subplot(5,1,2);
boxplot(diffTheta',groups);hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),npgCMap(j,:),'FaceAlpha',.4);
    h(j).Color = 'k'; % set box edge color
    h(j).MarkerEdgeColor=[0.8 0.8 0.8];
end
line([0 nConds+1],[100 100],'LineStyle','--','Color','k'); box off; set(gca,'TickDir','out');
title('theta'); ylabel('power ratio smallest/largest [%]');


subplot(5,1,3);
boxplot(diffSpindles',groups);hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),npgCMap(j,:),'FaceAlpha',.4);
    h(j).Color = 'k'; % set box edge color
    h(j).MarkerEdgeColor=[0.8 0.8 0.8];
end
line([0 nConds+1],[100 100],'LineStyle','--','Color','k'); box off; set(gca,'TickDir','out');
title('sigma'); ylabel('power ratio smallest/largest [%]');

subplot(5,1,4);
boxplot(diffTLastW',groups);hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),npgCMap(j,:),'FaceAlpha',.4);
    h(j).Color = 'k'; % set box edge color
    h(j).MarkerEdgeColor=[0.8 0.8 0.8];
end
line([0 nConds+1],[100 100],'LineStyle','--','Color','k'); box off; set(gca,'TickDir','out');
title('t since last wake'); ylabel('ratio smallest/largest [%]');

subplot(5,1,5);
boxplot(diffNRDur',groups);hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),npgCMap(j,:),'FaceAlpha',.4);
    h(j).Color = 'k'; % set box edge color
    h(j).MarkerEdgeColor=[0.8 0.8 0.8];
end
line([0 nConds+1],[100 100],'LineStyle','--','Color','k'); box off; set(gca,'TickDir','out');
title('preceeding NREM duration'); ylabel('ratio smallest/largest [%]');


% [coeff,score,latent,tsquared,explained,mu] = pca(currSpectra{5});
% % try multilinear regression 
% spectra =cat(1,currSpectra{:}); 
% forModel = array2table(spectra);
% forModel.tArousal = currTArousal{1}';
% test = fitlm(forModel);