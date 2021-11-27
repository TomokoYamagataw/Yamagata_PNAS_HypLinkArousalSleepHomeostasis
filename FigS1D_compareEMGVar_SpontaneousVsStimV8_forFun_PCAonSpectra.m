clear all
pathToData='F:\Tomoko_OptoStim\EMGAnalysis_matfiles';
pathFigs = 'F:\Tomoko_OptoStim\figs_10Hz\spontWakeVsStimWake'; % path to save figure
conditions={'TenHz24hWithExclusion','TenHz24hGFP','TenHzNonLPO','twentyHz24h'};
pathToHelpfulFunctions = 'C:\Users\Martin Kahn\Desktop\Scripts_and_Code\Matlab';
%colormap
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 240 228 66;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue,yellow
npgCMap = [npgCMap;npgCMap];
tomokoColors =[0 0 255;0 192 0; 249 64 64]./255;

%% add paths 
addpath(genpath(pathToHelpfulFunctions));
%% extract EMG variance and spectra in spontaneous wake compared to wake with stim


%%
nConds = numel(conditions);
[wakeSpectra_fro,wakeSpectra_occ] = deal(cell(nConds,1));
for condCnt = 1:nConds
    load(fullfile(pathToData,conditions{condCnt})) % load current condition   
    sR = settings.sR;
    nMice = numel(EMGs);
    for mouseCnt = 1:nMice
        mouseCnt
        % extract the relevant data from the loaded mat file
        currEMG= EMGs{mouseCnt} ;
        currCond=stims{mouseCnt} ;
        stimEpochs = floor(currCond/(sR*4)); % get epoch idx for each stimulation onset
        allStimEpochs = cell2mat(arrayfun(@(x)stimEpochs(x,1):stimEpochs(x,2)',1:size(stimEpochs,1),'un',0));
        currScoring=scorings{mouseCnt} ;
        currScoring2=scorings2{mouseCnt} ; % occipital derivation
        
        % convert wake spectra into PCA space 
        wakeOrNREM = [currScoring.w;currScoring.nr];
        [~,wakeOrNREMStim,~] = intersect(wakeOrNREM,allStimEpochs);
        currFroSpect = log10(cat(1,currScoring.spectr(currScoring.w,:),currScoring.spectr(currScoring.nr,:)));
        currOccSpect = log10(cat(1,currScoring2.spectr(currScoring.w,:),currScoring2.spectr(currScoring.nr,:)));
        currFroSpect = cat(2,currFroSpect,currOccSpect);
        currFroSpect(sum(isinf(currFroSpect),2)>0,:)=0;

         nWake = length(currScoring.w);
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(currFroSpect);
        figure;
        scatter3(SCORE(1:nWake,1),SCORE(1:nWake,2),SCORE(1:nWake,3),4,'k','MarkerFaceAlpha',0.3,'MarkerFaceColor','k','MarkerEdgeColor','none'); hold on;
        scatter3(SCORE(nWake+1:end,1),SCORE(nWake+1:end,2),SCORE(nWake+1:end,3),4,'r','MarkerFaceAlpha',0.3,'MarkerFaceColor',npgCMap(2,:),'MarkerEdgeColor','none')
        scatter3(SCORE(wakeOrNREMStim,1),SCORE(wakeOrNREMStim,2),SCORE(wakeOrNREMStim,3),4,'b','MarkerFaceAlpha',0.3,'MarkerFaceColor',npgCMap(3,:),'MarkerEdgeColor','none')
    end
        
        % get wake epochs and epochs with stimulation
        currStates = rmfield(currScoring,{'bastend','spectr','derivation'}); % remove non-scoring realted fields
        hypno = struct2cell(structfun(@max,currStates,'un',0)); % make empty array according to number of epochs
        hypno(cellfun(@isempty,hypno))=[];
        hypno= nan(max(cell2mat(hypno)),1);
        hypno(currStates.w)=1;
        hypno(currStates.w1)=1.1;
        hypno(currStates.nr)=2;
        hypno(currStates.nr2)=2.1;
        hypno(currStates.r)=3;
        hypno(currStates.r)=3.1;
        hypno(currStates.mt)=4;
        hypno(currStates.ma)=4.1;
        
        if stimEpochs(end,1)+44>length(hypno)
            stimEpochs(end,:)=[];
        end
        % now find stimulations which are surrounded by 2 minutes of wakefulness on each side
        if stimEpochs(1,1)<16
            stimEpochs(1,:)=[];
        end
        stimHypno = arrayfun(@(x) hypno(x-15:x+44),stimEpochs(:,1),'un',0); % get hypnogram surrounding each stim
        stimSpectra_fro = arrayfun(@(x) currScoring.spectr(x-15:x+44,:),stimEpochs(:,1),'un',0); % get hypnogram surrounding each stim
        stimSpectra_occ = arrayfun(@(x) currScoring2.spectr(x-15:x+44,:),stimEpochs(:,1),'un',0); % get hypnogram surrounding each stim
        
        stimOnlyWake = cellfun(@(x)sum(x>1.5)==0,stimHypno); % stimulations entirely surrounded by waking
        stimDuringL = stimEpochs(:,1)<(12*3600/4);
        stimsToUse=find(stimOnlyWake);
        
        
        
        [wakeSpectra_fro{condCnt}{mouseCnt},wakeSpectra_occ{condCnt}{mouseCnt}] = deal(nan(numel(stimsToUse),size(stimSpectra_fro{1},1),size(stimSpectra_fro{1},2)));
        for stimCnt = 1 :numel(stimsToUse)
            currHypno = stimHypno{stimsToUse(stimCnt)};
            currSpectra_fro = stimSpectra_fro{stimsToUse(stimCnt)};
            currSpectra_fro(currHypno~=1,:)= nan;
            wakeSpectra_fro{condCnt}{mouseCnt}(stimCnt,:,:)=currSpectra_fro;
            %same for occipital
            currSpectra_occ = stimSpectra_occ{stimsToUse(stimCnt)};
            currSpectra_occ(currHypno~=1,:)= nan;
            wakeSpectra_occ{condCnt}{mouseCnt}(stimCnt,:,:)=currSpectra_occ;
        end
        
    end
end
%% plot sepctrogram and average spectra
tInMin = (-15:44).*4-4; % from epochs to seconds
tInMin = tInMin./60;

[meanWSpectra_fro,meanWSpectra_occ]=deal(cell(nConds,1));
for condCnt = 1:nConds
    meanWSpectra_fro{condCnt} = cellfun(@(x)squeeze(nanmean(x,1)),wakeSpectra_fro{condCnt},'un',0);
    meanWSpectra_fro{condCnt} = cat(3,meanWSpectra_fro{condCnt}{:});
    meanWSpectra_occ{condCnt} = cellfun(@(x)squeeze(nanmean(x,1)),wakeSpectra_occ{condCnt},'un',0);
    meanWSpectra_occ{condCnt} = cat(3,meanWSpectra_occ{condCnt}{:});
end

freqs = 0:0.25:30;

epochsNoStim =  or(tInMin<0,tInMin>2);

%% plot individual mice 
occSpectraGFP = wakeSpectra_fro{1}; 
occSpectraGFP = cellfun(@(x) squeeze(nanmean(x,1)),occSpectraGFP,'un',0);
figure;
stimEpochs = 17:44; 
postStimEpochs = 46:size(occSpectraGFP{1},1);
preStimEpochs = 1:15;

for mouseCnt = 1:numel(occSpectraGFP)
    subplot(1,7,mouseCnt);
    currPreStim = squeeze(nanmean(occSpectraGFP{mouseCnt}(1:15,:),1)); 
    currStim = squeeze(nanmean(occSpectraGFP{mouseCnt}(postStimEpochs,:),1)); 
    plot(freqs,currPreStim,'Color',npgCMap(1,:)); hold on
    plot(freqs,currStim,'Color',npgCMap(2,:)); 
    legend('prestim','stim'); xlabel('freq [Hz]'); ylabel('power [uV^2/Hz]')
    set(gca,'Yscale','log')
    xlim([0,20])
end

%% plot absolute average spectra across mice
% define qwhich epochs to use 
stimEpochs = 17:44; 
postStimEpochs = 48:48+7;
% postStimEpochs = 48+7:size(meanWSpectra_occ{1},1);
preStimEpochs = 8:15;
% postStimEpochs = stimEpochs;

toPlot={meanWSpectra_fro{1},meanWSpectra_occ{1},meanWSpectra_fro{3},meanWSpectra_occ{3},meanWSpectra_fro{4},meanWSpectra_occ{4},meanWSpectra_fro{2},meanWSpectra_occ{2}};
titleStrings = {'fro_LPO','occ_LPO','fro_nLPO','occ_nLPO','fro_20Hz','occ_20Hz','fro_GFP','occ_GFP'};
nPlots = numel(toPlot); 

figure;
for plotCnt = 1:nPlots
    
    subplot(nPlots/2,2,plotCnt)
    currPlt = toPlot{plotCnt};
    grandAvgPreStim = squeeze(mean(currPlt(preStimEpochs,:,:),1));
%     grandSEMPreStim = 1.57*(iqr(grandAvgPreStim,2)./sqrt(nMice));% use CI for median:.
    grandSEMPreStim = std(grandAvgPreStim,[],2)./sqrt(nMice);
    
    grandAvgPreStim = mean(grandAvgPreStim,2);
    
    grandAvgStim = squeeze(median(currPlt(postStimEpochs,:,:),1)); % being conservative with the edges here
%     grandSEMStim = 1.57*(iqr(grandAvgStim,2)./sqrt(nMice));% use CI for median:.
    grandSEMStim = std(grandAvgStim,[],2)./sqrt(nMice);
    grandAvgStim = mean(grandAvgStim,2);
    
    p1 = patchMeUp(freqs,grandAvgPreStim,grandSEMPreStim,npgCMap(1,:),0.3); hold on;
    p2 = patchMeUp(freqs,grandAvgStim,grandSEMStim,npgCMap(2,:),0.3);
    xlabel('frequency [hz]'); ylabel('power [uV^2/Hz]'); legend([p1,p2],{'pre-stim','postStim'});legend('boxoff')
    set(gca,'TickDir','out','YScale','log'); box off; title(titleStrings{plotCnt},'interpreter','none');
    
end
%% plot proportional spectra (i.e. % of total power) across mice
% define qwhich epochs to use 
stimEpochs = 17:44; 
postStimEpochs = 48:47+7;
preStimEpochs = 8:15;

toPlot={meanWSpectra_fro{1},meanWSpectra_occ{1},meanWSpectra_fro{3},meanWSpectra_occ{3},meanWSpectra_fro{2},meanWSpectra_occ{2}};
titleStrings = {'fro_LPO','occ_LPO','fro_nLPO','occ_nLPO','fro_GFP','occ_GFP'};
nPlots = numel(toPlot); 

figure;
for plotCnt = 1:nPlots
    
    subplot(nPlots/2,2,plotCnt)
    currPlt = toPlot{plotCnt};
    currPlt = 100.*currPlt./sum(currPlt,2);
    grandAvgPreStim = squeeze(mean(currPlt(preStimEpochs,:,:),1));
%     grandSEMPreStim = 1.57*(iqr(grandAvgPreStim,2)./sqrt(nMice));% use CI for median:.
    grandSEMPreStim = std(grandAvgPreStim,[],2)./sqrt(nMice);
    
    grandAvgPreStim = mean(grandAvgPreStim,2);
    
    grandAvgStim = squeeze(mean(currPlt(postStimEpochs,:,:),1)); % being conservative with the edges here
%     grandSEMStim = 1.57*(iqr(grandAvgStim,2)./sqrt(nMice));% use CI for median:.
    grandSEMStim = std(grandAvgStim,[],2)./sqrt(nMice);
    grandAvgStim = mean(grandAvgStim,2);
    
    p1 = patchMeUp(freqs,grandAvgPreStim,grandSEMPreStim,npgCMap(1,:),0.3); hold on;
    p2 = patchMeUp(freqs,grandAvgStim,grandSEMStim,npgCMap(2,:),0.3);
    xlabel('frequency [hz]'); ylabel('power [uV^2/Hz]'); legend([p1,p2],{'pre-stim','postStim'});legend('boxoff')
    set(gca,'TickDir','out','YScale','log'); box off; title(titleStrings{plotCnt},'interpreter','none');
    
end



%% plot relative average spectra across mice
% put GFP and Chr2 on the same graph 

titleStrings = {'frontal EEG','occipital EEG'};
currCond = 1; 
toPlot={meanWSpectra_fro,meanWSpectra_occ};

figure;
p = {[],[]};
for plotCnt = 1:2
    for condCnt = [1,3,2]
        subplot(2,1,plotCnt)
        currCond = toPlot{plotCnt}{condCnt};
        grandAvg = mean(currCond(postStimEpochs,:,:),1)./mean(currCond(preStimEpochs,:,:),1);
        grandAvg = squeeze(10.*log10(grandAvg)); % convert to decibel
        %     grandSEMPreStim = 1.57*(iqr(grandAvgPreStim,2)./sqrt(nMice));% use CI for median:.
        grandSEM= std(grandAvg,[],2)./sqrt(size(grandAvg,2));
        grandAvg = mean(grandAvg,2);
        p{condCnt} = patchMeUp(freqs,grandAvg,grandSEM,npgCMap(condCnt,:),0.3); hold on;
        line(get(gca,'XLim'),[0 0 ],'Color','k','linestyle','--');
    end
    xlabel('frequency [hz]'); ylabel(' power stim/prestim [dB]'); legend([p{1},p{2},p{3}],{'LPO','GFP','nonLPO'});legend('boxoff')
    set(gca,'TickDir','out'); box off; title(titleStrings{plotCnt},'interpreter','none');
end
    
%% make boxplots for delta and theta power 
% total of 8 conditions: fro vs occ (2) SWA and theta (2) for gfp and and Chr2 (2) - 2x 2 x2 = 8 
thetaRange = and(freqs>=6,freqs<=9);
SWARange = and(freqs>0.5,freqs<=4);

toPlot{1}  = cellfun(@(x)squeeze(mean(x(:,SWARange,:),2)),meanWSpectra_fro,'un',0);
toPlot{2}  = cellfun(@(x)squeeze(mean(x(:,SWARange,:),2)),meanWSpectra_occ,'un',0);
toPlot{3}  = cellfun(@(x)squeeze(mean(x(:,thetaRange,:),2)),meanWSpectra_fro,'un',0);
toPlot{4}  = cellfun(@(x)squeeze(mean(x(:,thetaRange,:),2)),meanWSpectra_occ,'un',0);
% toPlot = cellfun(@(x) cat(3,x{:}),toPlot,'un',0);
titleStrings = {'SWA_fro','SWA_occ','theta_fro','theta_occ'};
figure; 
for plotCnt = 1:4
    % get mean power pre-stim and during stim 
    currPre = cellfun(@(x) squeeze(mean(x(preStimEpochs,:,:),1)),toPlot{plotCnt},'un',0); 
    currStim = cellfun(@(x) squeeze(mean(x(postStimEpochs,:,:),1)),toPlot{plotCnt},'un',0); 

    currChR2 = [currPre{1};currStim{1}]';
    currGFP = [currPre{2};currStim{2}]';

    subplot(2,2,plotCnt)
    groups = [ones(size(currChR2,1),1);2.*ones(size(currChR2,1),1);3.*ones(size(currGFP,1),1);4.*ones(size(currGFP,1),1)];
    boxplot([currChR2(:);currGFP(:)],groups); hold on;
    h = findobj(gca,'Tag','Box');
    h2= findobj(gca,'Tag','Outliers');
    colors = {npgCMap(1,:),npgCMap(1,:),npgCMap(2,:),npgCMap(2,:)};
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors{j},'FaceAlpha',.4);
        h(j).Color = colors{j}; % set box edge color
        h2(j).MarkerEdgeColor=[0.8 0.8 0.8];
    end
    p1= plot(1:2,currChR2,'o-','Color',npgCMap(2,:),'MarkerFaceColor',npgCMap(2,:),'MarkerSize',2);
    p2 = plot(3:4,currGFP,'o-','Color',npgCMap(1,:),'MarkerFaceColor',npgCMap(1,:),'MarkerSize',2);

    title(titleStrings{plotCnt},'interpreter','none'); box off; legend([p1(1) p2(1)],{'Chr2','GFP'}); legend('boxoff')
    ylabel('power density (uV^2/Hz)'); 
    
    YPosTxt = max(get(gca,'YLim'))+ 1;
    pChR2=signrank(currChR2(:,1),currChR2(:,2));
    pGFP=signrank(currGFP(:,1),currGFP(:,2));
    text(1,YPosTxt,['p=',num2str(round(pChR2,3))]);
    text(3,YPosTxt,['p=',num2str(round(pGFP,3))]);

end

