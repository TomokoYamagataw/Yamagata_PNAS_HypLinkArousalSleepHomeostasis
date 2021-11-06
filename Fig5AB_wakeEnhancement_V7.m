% In V4, I tried to explain the differences in NREM time after SD by wake/sleep bout duration. But it's really hard because most
% bouts don't nicely restrict themselves to the 30 minutes of interest
% in V5 I'm trying to ask whether thwe differences in SWA can simply be attributed to the difference in time asleep afterr SD.
% I.e. my theory is that what happens is that opto-SD just accumulates more awake time and thus sleep pressure than non-opoto SD

clear all
pathToData='F:\Tomoko_OptoStim\EMGAnalysis_matfiles';
pathFigs = 'F:\Tomoko_OptoStim\figsWakeEnhancement'; % path to save figure
pathToHelpfulFunctions = 'C:\Users\Martin Kahn\Desktop\Scripts_and_Code\Matlab';
%colormap
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 240 228 66;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue,yellow
npgCMap = [npgCMap;npgCMap];
tomokoColors =[0 0 255;0 192 0; 249 64 64]./255;

% freq range for theta and SWA according to 0.25 Hz bins
thetaRange = 25:37;  % 6 - 9 Hz
SWARange = 3:17; % 0.5 - 4 Hz
clrnLPO = [83,182,200]./255;  %define colors according to tomookos plot
clrGFP = [85,160,251]./255;
colorsToUse = {clrGFP,'b',clrnLPO};
% V6: compare LPO and nLPO and particularly compare stim-waking vs nostim waking
% V7: adjust figures so only a few selected plots are made and directly show comparison between nLPO and LPO 
%% load EMG and sleep scoring
addpath(genpath(pathToHelpfulFunctions));
pathFigs = 'F:\Tomoko_OptoStim\figsWakeEnhancement';
paths.extractedSignals ='F:\Tomoko_OptoStim\EMGAnalysis_matfiles';
SDOptoLPO=load(fullfile(paths.extractedSignals,'SDOpto_LPO'));
SDCtrlLPO=load(fullfile(paths.extractedSignals,'SDctrl_LPO'));
SDOptonLPO=load(fullfile(paths.extractedSignals,'SDOpto_nLPO'));
SDCtrlnLPO=load(fullfile(paths.extractedSignals,'SDctrl_nLPO'));
mNoFro = 3;  % GDCh 17 has no frontal derivation so we'll not use it for spectra but for the sleep dur etc
%% get wake bouts and NREM bouts and then get the corresponding spectra
[wSpectra_occ,wSpectra_fro,NREMSpectra_occ,NREMSpectra_fro,wakeBoutOnOffsets,NREMboutOnOffsets,wBoutDuration,....
    nRBoutDuration,nremEpochs,wakeEpochs,spectraFro,spectraOcc,stimSpectraFro,stimSpectraOcc] = deal(cell(2,1));


for condCnt = 1:2
    if condCnt == 1
        currCond = {SDOptoLPO,SDCtrlLPO};
    else
        currCond = {SDOptonLPO,SDCtrlnLPO};
    end
    
    nMice = numel(currCond{1}.EMGs);
    
    [wSpectra_occ{condCnt},wSpectra_fro{condCnt},NREMSpectra_occ{condCnt},NREMSpectra_fro{condCnt},wakeBoutOnOffsets{condCnt}...
        ,NREMboutOnOffsets{condCnt},wBoutDuration{condCnt},nRBoutDuration{condCnt},nremEpochs{condCnt},...
        wakeEpochs{condCnt},spectraFro{condCnt},spectraOcc{condCnt},stimSpectraFro{condCnt},stimSpectraOcc{condCnt}] = deal(cell(nMice,2));
    
    
    
    for mouseCnt = 1:nMice
        for optoOrCtrl = 1:2
            currExp = currCond{optoOrCtrl};
            % extract the relevant data from the loaded mat file
            currScoring=currExp.scorings{mouseCnt} ; % frontal
            currScoring2=currExp.scorings2{mouseCnt} ; % occipital derivation
            spectra_occ = currScoring2.spectr;
            spectra_fro= currScoring.spectr;
            spectraFro{condCnt}{mouseCnt,optoOrCtrl}=spectra_fro;
            spectraOcc{condCnt}{mouseCnt,optoOrCtrl}=spectra_occ;
            sR = currExp.settings.sR;
            
            % use stimTimings from opto experiment for both days
            if optoOrCtrl ==1
                currStims = currExp.stims{mouseCnt};
                stimEpochs = floor(currStims/(sR*4)); % get epoch idx for each stimulation onset
                allStimEpochs = cell2mat(arrayfun(@(x)stimEpochs(x,1):stimEpochs(x,2)',1:size(stimEpochs,1),'un',0));
            end
            epochsToUse = arrayfun(@(x) intersect(x:x+15,currScoring.w),stimEpochs(:,1),'un',0);
            currStimSpectra_fro = cellfun(@(x) nanmean(currScoring.spectr(x,:),1),epochsToUse,'un',0); % get hypnogram surrounding each stim
            currStimSpectra_fro=cat(1,currStimSpectra_fro{:});
            currStimSpectra_occ = cellfun(@(x) nanmean(currScoring2.spectr(x,:),1),epochsToUse,'un',0); % get hypnogram surrounding each stim
            currStimSpectra_occ=cat(1,currStimSpectra_occ{:});
            
            
            stimSpectraFro{condCnt}{mouseCnt,optoOrCtrl} = currStimSpectra_fro;
            stimSpectraOcc{condCnt}{mouseCnt,optoOrCtrl} = currStimSpectra_occ;
            
            %unscored epochs
            allScored=sort([currScoring.ma;currScoring.mt;currScoring.nr;currScoring.nr2;currScoring.r;currScoring.r3;currScoring.w;currScoring.w1]);
            unscored=find(arrayfun(@(x)sum(allScored==x)==0,1:numel(allScored)));
            
            % get wake epochs and epochs with stimulation
            wEpochs = sort([currScoring.w;currScoring.w1;currScoring.mt;currScoring.ma]);
            wakeBouts = [0;find(abs(diff(wEpochs))~=1);numel(wEpochs)]; % finds the end ioxs of each wake bout
            wakeBouts = [wakeBouts(1:end-1)+1,wakeBouts(2:end)];
            
            
            
            % big monkeypatch: there are some unscored epochs, we must make sure that wake and sleep bouts are not artificially
            % interrupted by them
            for boutCnt = 1:size(wakeBouts,1)-1
                if boutCnt > size(wakeBouts,1)-1
                    break
                end
                inBetweenEpochs = wEpochs(wakeBouts(boutCnt,2))+1:wEpochs(wakeBouts(boutCnt+1,1))-1; % epoch IDxs between curretn wake bouts
                unscoredEps = intersect(unscored,inBetweenEpochs); % find number of in between epochs that are unscored
                if numel(unscoredEps)==numel(inBetweenEpochs) % if all epochs between wake bouts are unscored, assume continuous waking
                    wakeBouts(boutCnt,2) =  wakeBouts(boutCnt+1,2);
                    wakeBouts(boutCnt+1,:)=[];
                    display(['excluded',num2str(boutCnt),'mouse',num2str(mouseCnt)])
                end
            end
            wakeEpochs{condCnt}{mouseCnt,optoOrCtrl}=currScoring.w;
            %         wakeBouts(boutDurs<15,:)=[]; % exclude wakefulness bouts shorter than 1 min
            
            nrEpochs = sort([currScoring.nr;currScoring.nr2;currScoring.r;currScoring.r3;]);
            nrBouts= [0;find(abs(diff(nrEpochs))~=1);numel(nrEpochs)]; % finds the end ioxs of each wake bout
            nrBouts = [nrBouts(1:end-1)+1,nrBouts(2:end)];
            nremEpochs{condCnt}{mouseCnt,optoOrCtrl} = currScoring.nr;
            
            
            % current problem: NREM bouts can either start after rem or after wake, so just define them as starting after wake endss
            % get spectra during artefact free wake
            [wSpectra_occ{condCnt}{mouseCnt,optoOrCtrl},wSpectra_fro{condCnt}{mouseCnt,optoOrCtrl},NREMSpectra_occ{condCnt}{mouseCnt,optoOrCtrl},NREMSpectra_fro{mouseCnt,condCnt}] = deal(cell(size(wakeBouts,1),1));
            [wBoutDuration{condCnt}{mouseCnt,optoOrCtrl},nRBoutDuration{condCnt}{mouseCnt,optoOrCtrl}] = deal(nan(size(wakeBouts,1),1));
            for boutCnt = 1:size(wakeBouts,1)
                currBout = wEpochs(wakeBouts(boutCnt,1):wakeBouts(boutCnt,2)); % epochs of current  wake bout
                nrBoutIDx = find(nrEpochs(nrBouts(:,1))>currBout(end),1); % which NREM bout follows this wake bout
                if isempty(nrBoutIDx) % if the last wake bout is not followed by a NREM bout
                    break
                end
                currNRBout = nrEpochs(nrBouts(nrBoutIDx,1):nrBouts(nrBoutIDx,2));
                
                
                wBoutDuration{condCnt}{mouseCnt,optoOrCtrl}(boutCnt)= (size(currBout,1)*4)./60;
                nRBoutDuration{condCnt}{mouseCnt,optoOrCtrl}(boutCnt)= (size(currNRBout,1)*4)./60; % in minutes, get duration of bout
                
                % store bout onset time
                NREMboutOnOffsets{condCnt}{mouseCnt,optoOrCtrl}(boutCnt,:)=(currNRBout([1,end])*4)/60; % in minutes
                wakeBoutOnOffsets{condCnt}{mouseCnt,optoOrCtrl}(boutCnt,:)=(currBout([1,end])*4)/60; % in minutes
            end
        end
    end
end

%% plot wake spectra during stim 
%%compare wake spectra stim day vs nonstim day
AvgFroLPO_bl = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraFro{1}(:,2),'un',0));
AvgFroLPO_stim = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraFro{1}(:,1),'un',0));
AvgFronLPO_bl = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraFro{2}(:,2),'un',0));
AvgFronLPO_stim = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraFro{2}(:,1),'un',0));

AvgOccLPO_bl = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraOcc{1}(:,2),'un',0));
AvgOccLPO_stim = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraOcc{1}(:,1),'un',0));
AvgOccnLPO_bl = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraOcc{2}(:,2),'un',0));
AvgOccnLPO_stim = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraOcc{2}(:,1),'un',0));

figure;
freq = 0:0.25:30;
for condCnt = 1:2
    AvgFroLPO_bl = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraFro{condCnt}(:,2),'un',0));
    AvgFroLPO_stim = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraFro{condCnt}(:,1),'un',0));
    meanFro_ratio = 10.*log10(AvgFroLPO_stim./AvgFroLPO_bl);
    
    AvgOccLPO_bl = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraOcc{condCnt}(:,2),'un',0));
    AvgOccLPO_stim = cell2mat(cellfun(@(x)nanmean(x,1),stimSpectraOcc{condCnt}(:,1),'un',0));
    meanOcc_ratio = 10.*log10(AvgOccLPO_stim./AvgOccLPO_bl);
    
    
    subplot(2,1,1)
    plots{1,condCnt}=patchMeUp(freq,mean(meanFro_ratio),std(meanFro_ratio)./sqrt(size(meanFro_ratio,1)),colorsToUse{condCnt},0.3); hold on;
    p1=plot(freq,meanFro_ratio,'Color',colorsToUse{condCnt}); hold on;for i = 1:numel(p1); p1(i).Color(4)=0.3;end
    
    box off; xlabel('frequency [Hz]'); ylabel('power relative to baseline [dB]');title('frontal');
    xlim([0.5,25]); line(get(gca,'XLim'),[0 0 ],'Color','k','LineStyle','--'); set(gca,'TickDir','out')
    
    subplot(2,1,2)
    plots{2,condCnt}=patchMeUp(freq,mean(meanOcc_ratio),std(meanOcc_ratio)./sqrt(size(meanOcc_ratio,1)),colorsToUse{condCnt},0.3); hold on;
    p1=plot(freq,meanOcc_ratio,'Color',colorsToUse{condCnt}); hold on;for i = 1:numel(p1); p1(i).Color(4)=0.3;end
    
    
    box off; xlabel('frequency [Hz]'); ylabel('power relative to baseline [dB]');title('occipital');
    xlim([0.5,25]);line(get(gca,'XLim'),[0 0 ],'Color','k','LineStyle','--'); set(gca,'TickDir','out');
    
    
end

%% look at average spectra
fRawSpectra=figure('name','raw spectra');
labels={'SDctrl','SD_opto'};
freqs = 0.25:0.25:30; % let's ignore the 0 freq as it can have 0 power due to filtering
[epochs,normNremSpectra,meanNormNremSpectra,meanSpectraFro,meanSpectraOcc,allNremSpectraFro,...
    allNremSpectraOcc] = deal(cell(2,1));

for condCnt = 1:2

    if condCnt == 1
        currCond = {SDOptoLPO,SDCtrlLPO};
    else
        currCond = {SDOptonLPO,SDCtrlnLPO};
    end
        nMice = numel(currCond{1}.EMGs);
                    [epochs{condCnt},normNremSpectra{condCnt}] = deal(struct());
        [meanNormNremSpectra{condCnt},meanSpectraFro{condCnt},meanSpectraOcc{condCnt}] = deal(nan(nMice,2,9,numel(freqs)));
        [allNremSpectraFro{condCnt},allNremSpectraOcc{condCnt}] = deal(cell(nMice,2));

    for optoOrCtrl = 1:2

        for mouseCnt = 1:nMice
            
            if and(mouseCnt == mNoFro,condCnt==1) % skip the mouse without frontal derivation
                continue
            end
            

            % get NREM spectra before stim
            currNrEpochs =  nremEpochs{condCnt}{mouseCnt,optoOrCtrl};
            currWEpochs =  wakeEpochs{condCnt}{mouseCnt,optoOrCtrl};
            
            % get epoch-wise EMG variance           
            currEMG = currCond{optoOrCtrl}.EMGs{mouseCnt};
            binSize = sR/(1/currCond{optoOrCtrl}.settings.binSize_EMGVar); % binSize for EMG variance in samples            
            % calculate variance for each EMG epoch and also calculate mean amplitude of positive half-width (Vyazovskiy 2005)
            w_EMG = cell2mat(arrayfun(@(x) currEMG(x*sR+1:x*sR+sR*4),currWEpochs*4 -4,'un',0)); % trial-wise opto-evoked EMG
            w_EMGVar = squeeze(var(reshape(w_EMG,size(w_EMG,1),binSize,size(w_EMG,2)/binSize),[],2));
            w_EMG(w_EMG<0)=0;
            [PKS,LOCS] = arrayfun(@(x)findpeaks(w_EMG(x,:)),1:size(w_EMG,1),'un',0); % find peak idxs
            meanEMGPeaks = log(cellfun(@mean,PKS));
            
            w_EMGVar = log(mean(w_EMGVar,2));
            w_EMGVarPre= w_EMGVar(currWEpochs<(7200/4));
            w_EMGVarDuring= w_EMGVar(and(currWEpochs>((2*3600)/4),currWEpochs<((4*3600)/4)));
            
            % subdivide EMG into high and low epochs aguely according to Vyazovskiy 2005
            %         lowEMGIdxs= w_EMGVar<(median(w_EMGVar)-std(w_EMGVar)); %based on EMG variance
            lowEMGIdxs= transpose(meanEMGPeaks<(nanmean(meanEMGPeaks)-nanstd(meanEMGPeaks))); %based on EMG amplitude (Vyazovskiy 2005)
            
            % subdivide epochs
            epochs{condCnt}.pre = currNrEpochs(currNrEpochs<((2*3600)/4));
            epochs{condCnt}.post4h = currNrEpochs(and(currNrEpochs>((4*3600)/4),currNrEpochs<((8*3600)/4)));
            epochs{condCnt}.post2h = currNrEpochs(and(currNrEpochs>((4*3600)/4),currNrEpochs<((6*3600)/4)));
            epochs{condCnt}.Wpre = currWEpochs(currWEpochs<(7200/4));
            epochs{condCnt}.WDuring = currWEpochs(and(currWEpochs>((2*3600)/4),currWEpochs<((4*3600)/4)));
            epochs{condCnt}.Wpre_lowEMG = currWEpochs(and(currWEpochs<(7200/4),lowEMGIdxs));
            epochs{condCnt}.WDuring_lowEMG = currWEpochs(and(and(currWEpochs>((2*3600)/4),currWEpochs<((4*3600)/4)),lowEMGIdxs));
            epochs{condCnt}.Wpre_highEMG = currWEpochs(and(currWEpochs<(7200/4),lowEMGIdxs==0));
            epochs{condCnt}.WDuring_highEMG = currWEpochs(and(and(currWEpochs>((2*3600)/4),currWEpochs<((4*3600)/4)),lowEMGIdxs==0));
            
            % now get corresponding spectra
            currSpectraFro = spectraFro{condCnt}{mouseCnt,optoOrCtrl}(:,2:end);
            currSpectraFro = structfun(@(x) currSpectraFro(x,:),epochs{condCnt},'un',0);
            currMeanSpectraFro = structfun(@(x) nanmean(x,1),currSpectraFro,'un',0);
            
            currSpectraOcc = spectraOcc{condCnt}{mouseCnt,optoOrCtrl}(:,2:end);
            currSpectraOcc = structfun(@(x) currSpectraOcc(x,:),epochs{condCnt},'un',0);
            currMeanSpectraOcc = structfun(@(x) nanmean(x,1),currSpectraOcc,'un',0);
            
            % normalise spectra to mean pre stim baseline
            if optoOrCtrl >0
                meanPreBaseline =  nanmean(currMeanSpectraFro.pre(2:end));
            end
            allNremSpectraFro{condCnt}{mouseCnt,optoOrCtrl}=currSpectraFro;
            tmp=struct2cell(currMeanSpectraFro);
            meanSpectraFro{condCnt}(mouseCnt,optoOrCtrl,:,:) =cell2mat(tmp);
            
            allNremSpectraOcc{condCnt}{mouseCnt,optoOrCtrl}=currSpectraOcc;
            tmp=struct2cell(currMeanSpectraOcc);
            meanSpectraOcc{condCnt}(mouseCnt,optoOrCtrl,:,:) =cell2mat(tmp);
        end
    end
end
% linkaxes()

%% NREM spectra plot grand average (i.e. across mice) of spectra normalised in different ways
% Q: Does opto stim in addition to SD change spectra, specifically, does it change SWA?
% dimensions of meanSpectraFro: (mice,noStimStim,timing(see above),frequency)  
    
fDifferenceInSpectra=figure('name','difference in spectra');
hp2 = uipanel('Title','difference in spectra','FontSize',12);
conds = {'LPO','nLPO'};
for condCnt = 1:2
    pVals = nan(3,3);
    nMiceForThis = size(meanSpectraFro{condCnt},1);
    
    % normalise each condition to the mean of the 2h before SD
%     toPlot = 10.*log10(meanSpectraFro{condCnt}./nanmean(meanSpectraFro{condCnt}(:,:,1,:),4));
%     toPlot = 100.*(meanSpectraFro{condCnt}./nanmean(meanSpectraFro{condCnt}(:,2,1,:),4));
    toPlot = 100.*(meanSpectraFro{condCnt}./meanSpectraFro{condCnt}(:,:,1,:));

    yAxLabel = {'power [% pre-SD mean]'};
    currTitle=[conds{condCnt},'sepparately normalised to avg pwr 2h pre SD'];

    grandAvg = squeeze(nanmean(toPlot,1));
    grandSEM = squeeze(nanstd(toPlot,[],1)./sqrt(nMiceForThis));
    % difference between spectra before/after sleep dep/ sleep dep & opto
    diffMeanSpectra =  squeeze(toPlot(:,1,:,:) - toPlot(:,2,:,:));
    grandAvg_diff = squeeze(nanmean(diffMeanSpectra,1));
    grandSEM_diff = squeeze(nanstd(diffMeanSpectra,[],1)./sqrt(nMiceForThis));
    
    fNormalisedSpectra=figure('name',currTitle);
    hp1 = uipanel('Title',currTitle,'FontSize',12);
    titles={'2h before sleep dep','first 2h after SD','first 4h after SD'};
    splots = []; splots2=[];
    for sPlotCnt = 1:3
        figure(fNormalisedSpectra)
        sPlots(sPlotCnt)=subplot(1,3,sPlotCnt,'Parent',hp1); hold on;
        p1 =patchMeUp(freqs,squeeze(grandAvg(2,sPlotCnt,:)),squeeze(grandSEM(1,sPlotCnt,:)),npgCMap(1,:),0.3); hold on;
        p2=patchMeUp(freqs,squeeze(grandAvg(1,sPlotCnt,:)),squeeze(grandSEM(2,sPlotCnt,:)),npgCMap(3,:),0.3); hold on;
        title(titles{sPlotCnt}); xlabel('frequency [hz]'); ylabel(yAxLabel);
        legend([p1,p2],{'SD only','opto + SD'},'autoupdate','off')
        
%         % add the indlividual lines for each mouse
%         p1=plot(repmat(freqs,size(toPlot,1),1)', squeeze(toPlot(:,2,sPlotCnt,:))','Color',npgCMap(1,:));
%         p2=plot(repmat(freqs,size(toPlot,1),1)',squeeze(toPlot(:,1,sPlotCnt,:))','Color',npgCMap(3,:));
%         for i = 1:numel(p1); p1(i).Color(4)=0.3; p2(i).Color(4)=0.3;end
        xlim([0 15]); ax = gca; %ax.YScale = 'log';
        
        
        figure(fDifferenceInSpectra)
        sPlots2(sPlotCnt)=subplot(1,3,sPlotCnt,'Parent',hp2); hold on;
        p1 =patchMeUp(freqs,squeeze(grandAvg_diff(sPlotCnt,:)),squeeze(grandSEM_diff(sPlotCnt,:)),colorsToUse{condCnt+1},0.3); hold on;
        xlabel('frequency [hz]'); ylabel(['difference in ',yAxLabel]);
        legend(p1,'difference SD vs SD + opto','autoupdate','off')
        
        % add the indlividual lines for each mouse
%         p1=plot(repmat(freqs,size(toPlot,1),1)', squeeze(diffMeanSpectra(:,sPlotCnt,:))','Color',npgCMap(1,:));
%         for i = 1:numel(p1); p1(i).Color(4)=0.3;end
        xlim([0 15]); ax = gca; %ax.YScale = 'log';
    end
    linkaxes(sPlots2);
    linkaxes(sPlots);
    set(gcf,'Position',[1939 360  1100  767]);
    %     saveas(gcf,fullfile(pathFigs,['avgSpectr',currTitle,'.fig']));
    %     saveas(gcf,fullfile(pathFigs,['avgSpectr',currTitle,'.png']));
    
    % plot overall effects on SWA
    figure();
    meanSWA = squeeze(nanmean(toPlot(:,:,:,6:16),4));
    grandMeanSWA = squeeze(nanmean(meanSWA,1));
    grandSEMSWA =  squeeze(nanstd(meanSWA,[],1)./sqrt(nMiceForThis));
    for sPlotCnt = 1:3
        errorbar([sPlotCnt-0.2,sPlotCnt+0.2],grandMeanSWA([2,1],sPlotCnt),grandSEMSWA([2,1],sPlotCnt),'ok','Color',npgCMap(7,:),'MarkerFaceColor',npgCMap(7,:)); hold on;
        plot([sPlotCnt-0.2,sPlotCnt+0.2],squeeze(meanSWA(:,[2,1],sPlotCnt)),'Color',[0.8 0.8 0.8])
        [~,pVals(plotCnt,sPlotCnt),~,STATS]=ttest(squeeze(meanSWA(:,2,sPlotCnt)),squeeze(meanSWA(:,1,sPlotCnt)));
    end
    xlim([0.2, 3.8]); ax = gca; ax.XTick = 1:3; ax.XTickLabel = {'2h before','2h after','4h after'};ylabel(yAxLabel);
    title(currTitle)
    
end
%% wake spectra plot grand average (i.e. across mice) of spectra normalised in different ways
% Q: Does opto stim in addition to SD change spectra, specifically, does it change SWA?
fSWA=figure();
pVals = nan(3,3);
for plotCnt = 3:3
    if plotCnt==1
        toPlot = meanSpectraOcc(:,:,4:end,:);%for the wake condition, just use the wake preSD and wake during SD ignore NREM
        yAxLabel = {'raw power [uV^2]'};
        currTitle='unnormalised';
    elseif plotCnt==2 % normalise each condition to the mean of the SD only condition 2h before SD, so both conditions normalised to same number
        toPlot =meanSpectraOcc(:,:,4:end,:);%for the wake condition, just use the wake preSD and wake during SD ignore NREM
        toPlot(:,:,1:2,:)=  100.*toPlot(:,:,1:2,:)./mean(toPlot(:,1,1,:),4); % when not sorting acc to EMG
        toPlot(:,:,3:4,:)=  100.*toPlot(:,:,3:4,:)./mean(toPlot(:,1,3,:),4); % normalise low EMG to low EMG SD only
        toPlot(:,:,5:6,:)=  100.*toPlot(:,:,5:6,:)./mean(toPlot(:,1,5,:),4); % normalise low EMG to low EMG SD only
        yAxLabel = {'power [% SD-only pre-SD mean]'};
        currTitle='both conditions normalised to mean power in 2h pre SD without opto';
    elseif plotCnt==3 % normalise each condition to the mean of the 2h before SD
        toPlot =meanSpectraOcc(:,:,4:end,:);%for the wake condition, just use the wake preSD and wake during SD ignore NREM
        toPlot(:,:,1:2,:)=  100.*toPlot(:,:,1:2,:)./mean(toPlot(:,1:2,1,:),4); % when not sorting acc to EMG
        toPlot(:,:,3:4,:)=  100.*toPlot(:,:,3:4,:)./mean(toPlot(:,1:2,3,:),4); % normalise low EMG to low EMG SD only
        toPlot(:,:,5:6,:)=  100.*toPlot(:,:,5:6,:)./mean(toPlot(:,1:2,5,:),4); % normalise low EMG to low EMG SD only
        yAxLabel = {'power [% pre-SD mean]'};
        currTitle='each condition sepparately normalised to mean power in 2h pre SD';
    end
    %     toPlot = toPlot(:,:,4:end,:); %for the wake condition, just use the wake pre and wake during
    fNormalisedSpectra=figure('name',currTitle);
    hp = uipanel('Title',currTitle,'FontSize',12);
    
    grandAvg = squeeze(mean(toPlot,1));
    grandSEM = squeeze(nanstd(toPlot,[],1)./sqrt(nMice));
    
    % difference between spectra before/after sleep dep/ sleep dep & opto
    diffMeanSpectra =  squeeze(toPlot(:,2,:,:) - toPlot(:,1,:,:));
    grandAvg_diff = squeeze(mean(diffMeanSpectra,1));
    grandSEM_diff = squeeze(nanstd(diffMeanSpectra,[],1)./sqrt(nMice));
    
    titles={'before','during','before-low','during-low','before-high','during-high'};
    splots = []; splots2=[];
    for condCnt = 1:6
        sPlots(condCnt)=subplot(2,6,condCnt,'Parent',hp); hold on;
        p1 =patchMeUp(freqs,squeeze(grandAvg(1,condCnt,:)),squeeze(grandSEM(1,condCnt,:)),npgCMap(1,:),0.3); hold on;
        p2=patchMeUp(freqs,squeeze(grandAvg(2,condCnt,:)),squeeze(grandSEM(2,condCnt,:)),npgCMap(3,:),0.3); hold on;
        title(titles{condCnt}); xlabel('frequency [hz]'); ylabel(yAxLabel);
        legend([p1,p2],{'SD only','opto + SD'},'autoupdate','off')
        
        % add the indlividual lines for each mouse
        p1=plot(repmat(freqs,nMice,1)', squeeze(toPlot(:,1,condCnt,:))','Color',npgCMap(1,:));
        p2=plot(repmat(freqs,nMice,1)',squeeze(toPlot(:,2,condCnt,:))','Color',npgCMap(3,:));
        for i = 1:numel(p1); p1(i).Color(4)=0.3; p2(i).Color(4)=0.3;end
        xlim([0 15]); ax = gca; %ax.YScale = 'log';
        
        sPlots2(condCnt)=subplot(2,6,condCnt+6,'Parent',hp); hold on;
        p1 =patchMeUp(freqs,squeeze(grandAvg_diff(condCnt,:)),squeeze(grandSEM_diff(condCnt,:)),npgCMap(1,:),0.3); hold on;
        xlabel('frequency [hz]'); ylabel(['difference in ',yAxLabel]);
        legend(p1,'difference SD vs SD + opto','autoupdate','off')
        
        % add the indlividual lines for each mouse
        p1=plot(repmat(freqs,nMice,1)', squeeze(diffMeanSpectra(:,condCnt,:))','Color',npgCMap(1,:));
        for i = 1:numel(p1); p1(i).Color(4)=0.3;end
        xlim([0 15]); ax = gca; %ax.YScale = 'log';
    end
    linkaxes(sPlots2);
    linkaxes(sPlots);
    set(gcf,'Position',[1939 360  1100  767]);
    %     saveas(gcf,fullfile(pathFigs,['avgSpectr',currTitle,'.fig']));
    %     saveas(gcf,fullfile(pathFigs,['avgSpectr',currTitle,'.png']));
    
    % plot overall effects on SWA
    figure(fSWA);
    subplot(1,3,plotCnt)
    meanSWA = squeeze(mean(toPlot(:,:,:,thetaRange),4));
    grandMeanSWA = squeeze(mean(meanSWA,1));
    grandSEMSWA =  squeeze(std(meanSWA,[],1)./sqrt(nMice));
    for condCnt = 1:2
        errorbar([condCnt-0.2,condCnt+0.2],grandMeanSWA(:,condCnt),grandSEMSWA(:,condCnt),'ok','Color',npgCMap(7,:),'MarkerFaceColor',npgCMap(7,:)); hold on;
        plot([condCnt-0.2,condCnt+0.2],squeeze(meanSWA(:,:,condCnt)),'Color',[0.8 0.8 0.8])
        [~,pVals(plotCnt,condCnt),~,STATS]=ttest(squeeze(meanSWA(:,1,condCnt)),squeeze(meanSWA(:,2,condCnt)));
    end
    xlim([0.2, 3.8]); ax = gca; ax.XTick = 1:3; ax.XTickLabel = {'2h before','2h after','4h after'};ylabel(yAxLabel);
    title(currTitle)
end

%% check why logging the power makes a difference
figure;
for mouseCnt = 1:nMice
    freq = 13; % so 4 Hz
    currLogged_Ctrl = allNremSpectraFro{mouseCnt,1};
    currLogged_Opto = allNremSpectraFro{mouseCnt,2};
    
    curUnLogged_Ctrl = exp(currLogged_Ctrl.post4h(:,freq))./mean(exp(currLogged_Ctrl.pre(:,freq)));
    curUnLogged_Opto =exp(currLogged_Opto.post4h(:,freq))./mean(exp(currLogged_Opto.pre(:,freq)));
    
    currLogged_Ctrl = currLogged_Ctrl.post4h(:,freq)./mean(currLogged_Ctrl.pre(:,freq));
    currLogged_Opto =currLogged_Opto.post4h(:,freq)./mean(currLogged_Opto.pre(:,freq));
    
    subplot(2,nMice,mouseCnt)
    [histcounts_Ctrl,edges]=histcounts(curUnLogged_Ctrl,100,'Normalization','probability');
    histcounts_Opto=histcounts(curUnLogged_Opto,edges,'Normalization','probability');
    plot(edges(2:end),histcounts_Ctrl,'k'); hold on;
    plot(edges(2:end),histcounts_Opto,'r'); hold on;
    line([mean(curUnLogged_Ctrl),mean(curUnLogged_Ctrl)],get(gca,'YLim'),'Color','k','LineStyle','--');
    line([mean(curUnLogged_Opto),mean(curUnLogged_Opto)],get(gca,'YLim'),'Color','r','LineStyle','--');
    line([median(curUnLogged_Ctrl),median(curUnLogged_Ctrl)],get(gca,'YLim'),'Color','k','LineStyle','-');
    line([median(curUnLogged_Opto),median(curUnLogged_Opto)],get(gca,'YLim'),'Color','r','LineStyle','-');
    box off; title('raw data'); xlabel('power [uV^2]'); ylabel('probability');
    ax = gca; ax.XScale = 'log';xlim([10^-1. 10^2]);
    
    subplot(2,nMice,mouseCnt+nMice)
    [histcounts_Ctrl,edges]=histcounts(currLogged_Ctrl,100,'Normalization','probability');
    histcounts_Opto=histcounts(currLogged_Opto,edges,'Normalization','probability');
    plot(edges(2:end),histcounts_Ctrl,'k'); hold on;
    plot(edges(2:end),histcounts_Opto,'r'); hold on;
    line([mean(currLogged_Ctrl),mean(currLogged_Ctrl)],get(gca,'YLim'),'Color','k','LineStyle','--');
    line([mean(currLogged_Opto),mean(currLogged_Opto)],get(gca,'YLim'),'Color','r','LineStyle','--');
    line([median(currLogged_Ctrl),median(currLogged_Ctrl)],get(gca,'YLim'),'Color','k','LineStyle','-');
    line([median(currLogged_Opto),median(currLogged_Opto)],get(gca,'YLim'),'Color','r','LineStyle','-');
    box off; title('logtransformed data'); xlabel('power [log]'); ylabel('probability');xlim([-1 2]);
end
%% do mice have more sleep after opto-SD? They have less! Does SWA mirror this trend?
nHours = 8;
binSize = 30; % in minutes
nBins = nHours*60/binSize;
[NRAmount_Ctrl,NRAmount_Opto,RAmount_Opto,RAmount_Ctrl,meanSWA_Ctrl,meanSWA_Opto,sumSWACtrl,sumSWAOpto] = deal(nan(nMice,nBins));
freqs = 0:0.25:30; % let's ignore the 0 freq as it can have 0 power due to filtering
SWARange = find(and(freqs>=0.5,freqs<=4));

for mouseCnt = 1:nMice
    
    currScoringCtrl=SDCtrl.scorings{mouseCnt} ; % frontal
    currScoringOpto=SDOpto.scorings{mouseCnt} ; % frontal
    NRCtrl = sort([currScoringCtrl.nr;currScoringCtrl.nr2]);
    RCtrl = sort([currScoringCtrl.r;currScoringCtrl.r3]);
    NROpto  = sort([currScoringOpto.nr;currScoringOpto.nr2]);
    ROpto = sort([currScoringOpto.r;currScoringOpto.r3]);
    
    %calculate amount of NREM after sleep dep in 20 minute bins
    for binCnt = 1:nBins
        currEdges = ceil([4*60+(binCnt-1)*binSize,4*60+binCnt*binSize].*60./4); % edges in epoch numbers (i.e. 4s bins)
        
        currNREp = NRCtrl(and(NRCtrl>=currEdges(1),NRCtrl<=currEdges(2))); %NR epochs in the current bin
        NRAmount_Ctrl(mouseCnt,binCnt)= (numel(currNREp).*4)./60; % count epochs and convert to seconds
        RAmount_Ctrl(mouseCnt,binCnt)= sum(and(RCtrl>=currEdges(1),RCtrl<=currEdges(2)));
        % now get SWA timecourse, look at average SWA per bin as well as summed SWA (i.e either we dont or we do factor in how much total NR sleep)
        meanSWA_Ctrl(mouseCnt,binCnt) = nanmean(nanmean(currScoringCtrl.spectr(currNREp,SWARange),2),1);
        sumSWACtrl(mouseCnt,binCnt) = nansum(nanmean(currScoringCtrl.spectr(currNREp,SWARange),2),1);
        
        
        % same as above for opto
        currNREp = NROpto(and(NROpto>=currEdges(1),NROpto<=currEdges(2))); %NR epochs in the current bin
        NRAmount_Opto(mouseCnt,binCnt)= (numel(currNREp).*4)./60; % count epochs and convert to minutes
        RAmount_Opto(mouseCnt,binCnt)= sum(and(ROpto>=currEdges(1),ROpto<=currEdges(2)));
        meanSWA_Opto(mouseCnt,binCnt) = nanmean(nanmean(currScoringOpto.spectr(currNREp,SWARange),2),1);
        sumSWAOpto(mouseCnt,binCnt) = nansum(nanmean(currScoringOpto.spectr(currNREp,SWARange),2),1);
        
        if sum(isnan(currScoringCtrl.spectr(:)))>0
            pause
        end
    end
    
end
% TTEST FIRST BIN
[p,h]=ttest(NRAmount_Ctrl(:,1),NRAmount_Opto(:,1))


% remove spectra for faulty frontal EEG in GDCh 17
[sumSWAOpto(mNoFro,:),sumSWACtrl(mNoFro,:),sumSWAOpto(mNoFro,:),meanSWA_Ctrl(mNoFro,:),meanSWA_Ctrl(mNoFro,:)]=deal(nan); % skip the mouse without frontal derivation

optoCumSumNR=cumsum(NRAmount_Opto,2);
ctrlCumSumNR=cumsum(NRAmount_Ctrl,2);
% optoCumSumR=cumsum(RAmount_Opto,2,'omitnan')./60;
% ctrlCumSumR=cumsum(RAmount_Ctrl,2,'omitnan' )./60;
optoSWACumSum=cumsum(sumSWAOpto,2);
ctrlSWACumSum=cumsum(sumSWACtrl,2);

figure;
toPlot = {{ctrlCumSumNR,optoCumSumNR},{NRAmount_Ctrl,NRAmount_Opto},{ctrlSWACumSum,optoSWACumSum},{meanSWA_Ctrl,meanSWA_Opto}};
labels={'cumulative time in NREM','mean time in NREM','cumulativer SWA','meaN SWA'};
tInMins = binSize:binSize:nHours*60;
for plotCnt = 1:4
    subplot(2,2,plotCnt)
    p1=patchMeUp(tInMins,nanmean(toPlot{plotCnt}{1},1),nanstd(toPlot{plotCnt}{1},1)./sqrt(nMice),'k',0.3); hold on;
    p2=patchMeUp(tInMins,nanmean(toPlot{plotCnt}{2},1),nanstd(toPlot{plotCnt}{2},1)./sqrt(nMice),'r',0.3); hold on;
    xlabel('minutes after SD end'); ylabel(labels{plotCnt}); l=legend([p1 p2],{'ctrl-SD','opto-SD'});legend('boxoff');
    set(gca,'TickDir','out');box off; title(labels{plotCnt}); l.Location = 'northwest';
end

saveas(gcf,fullfile(pathFigs,'timeCourses.fig'));
saveas(gcf,fullfile(pathFigs,'timeCourses.png'));

% plot same figure bot normalise each moiuse to itself to remove inter-mouse variability
figure;
toPlot{1} = {100.*ctrlCumSumNR./ctrlCumSumNR(:,end),100.*optoCumSumNR./ctrlCumSumNR(:,end)};
toPlot{3} = {100.*ctrlSWACumSum./ctrlSWACumSum(:,end),100.*optoSWACumSum./ctrlSWACumSum(:,end)};
toPlot{2} = {100.*NRAmount_Ctrl./mean(NRAmount_Ctrl,2),100.*NRAmount_Opto./mean(NRAmount_Ctrl,2)}; % normalise to mean of ctrl
toPlot{4} = {100.*meanSWA_Ctrl./mean(meanSWA_Ctrl,2),100.*meanSWA_Opto./mean(meanSWA_Ctrl,2)}; % normalise to mean of ctrl
labels={'cumulative time in NREM','mean time in NREM','cumulativer SWA','meaN SWA'};
yLabels = {'cum t asleep [% total ctrl]','mean t asleep [% mean ctrl]','cum SWA [% total ctrl]','mean SWA [% mean ctrl]'};
tInMins = binSize:binSize:nHours*60;
for plotCnt = 1:4
    subplot(2,2,plotCnt)
    p1=patchMeUp(tInMins,nanmean(toPlot{plotCnt}{1},1),nanstd(toPlot{plotCnt}{1},1)./sqrt(nMice),'k',0.3); hold on;
    p11=plot(tInMins,toPlot{plotCnt}{1},'k'); hold on; for i = 1:numel(p11); p11(i).Color(4)=0.2; end
    p2=patchMeUp(tInMins,nanmean(toPlot{plotCnt}{2},1),nanstd(toPlot{plotCnt}{2},1)./sqrt(nMice),'r',0.3); hold on;
    p22=plot(tInMins,toPlot{plotCnt}{2},'r'); hold on;for i = 1:numel(p22); p22(i).Color(4)=0.2; end
    
    xlabel('minutes after SD end'); ylabel(yLabels{plotCnt}); l=legend([p1 p2],{'ctrl-SD','opto-SD'});legend('boxoff');
    set(gca,'TickDir','out');box off; title(labels{plotCnt}); l.Location = 'northwest';
end

saveas(gcf,fullfile(pathFigs,'timeCoursesDifferencefig'));
saveas(gcf,fullfile(pathFigs,'timeCourseDifference.png'));

%% now just plot difference
figure;
toPlot{1} = 100.*optoCumSumNR./ctrlCumSumNR(:,end)-100.*ctrlCumSumNR./ctrlCumSumNR(:,end);
toPlot{3} = 100.*optoSWACumSum./optoSWACumSum(:,end)-100.*ctrlSWACumSum./ctrlSWACumSum(:,end);

toPlot{2} = 100.*NRAmount_Opto./mean(NRAmount_Ctrl,2)-100.*NRAmount_Ctrl./mean(NRAmount_Ctrl,2); % normalise to mean of ctrl
toPlot{4} = 100.*meanSWA_Opto./mean(meanSWA_Ctrl,2)-100.*meanSWA_Ctrl./mean(meanSWA_Ctrl,2); % normalise to mean of ctrl

labels={'cumulative time in NREM','mean time in NREM','cumulative SWA','mean SWA'};
yLabels = {'cum t asleep [% total ctrl]','mean t asleep [% mean ctrl]','cum SWA [% total ctrl]','mean SWA [% mean ctrl]'};
tInMins = binSize:binSize:nHours*60;
for plotCnt = 1:4
    subplot(2,2,plotCnt)
    p1=patchMeUp(tInMins,nanmean(toPlot{plotCnt},1),nanstd(toPlot{plotCnt},1)./sqrt(nMice),'k',0.3); hold on;
    p11 = plot(tInMins,toPlot{plotCnt},'k'); for i = 1:numel(p11); p11(i).Color(4)=0.2; end
    xlabel('minutes after SD end'); ylabel(yLabels{plotCnt});
    set(gca,'TickDir','out');box off; title(labels{plotCnt});
    line(get(gca,'XLim'),[0 0],'Color','k','LineStyle','--');
end

saveas(gcf,fullfile(pathFigs,'timeCoursesRelative.fig'));
saveas(gcf,fullfile(pathFigs,'timeCoursesRelative.png'));



%% Is there a difference in time to first NREM sleep? And can it explain the difference in SWA?
[tFirstNROpto,tFirstNRCtrl,] = deal(nan(nMice,1));
[spectraOptoPre,spectraOptoPost,spectraCtrlPre,spectraCtrlPost]= deal(nan(nMice,121));

maxTime = 100;
for mouseCnt = 1:nMice
    %calculate amount of NREM after sleep dep in 20 minute bins
    currEdges = [4*60,4*60+maxTime]; % edges in minutes
    currNRBouts = and(NREMboutOnOffsets{mouseCnt,2}>= currEdges(1),NREMboutOnOffsets{mouseCnt,2}<= currEdges(2));
    currNRBouts = find(sum(currNRBouts,2)>=2); % only analyse bouts that start and end in the time we look at
    tFirstNROpto(mouseCnt) = NREMboutOnOffsets{mouseCnt,2}(currNRBouts(1)) - 4*60;
    currNrEp =  nremEpochs{mouseCnt,2};
    currNrEpPre =  currNrEp(currNrEp<((120*60)/4));
    currNrEpPost=  currNrEp(and(currNrEp>((4*60*60)/4),currNrEp<((5*60*60)/4)));
    spectraOptoPre(mouseCnt,:) = mean(SDOpto.scorings{mouseCnt}.spectr(currNrEpPre,:),1);
    spectraOptoPost(mouseCnt,:) = mean(SDOpto.scorings{mouseCnt}.spectr(currNrEpPost,:),1);
    
    currNRBouts = and(NREMboutOnOffsets{mouseCnt,1}>= currEdges(1),NREMboutOnOffsets{mouseCnt,1}<= currEdges(2));
    currNRBouts = sum(currNRBouts,2)>=2; % only analyse bouts that start and end in the time we look at
    tFirstNRCtrl(mouseCnt) = NREMboutOnOffsets{mouseCnt,1}(find(currNRBouts(:,1),1)) - 4*60;
    currNrEp =  nremEpochs{mouseCnt,1};
    currNrEpPre =  currNrEp(currNrEp<((120*60)/4));
    currNrEpPost=  currNrEp(and(currNrEp>((4*60*60)/4),currNrEp<((5*60*60)/4)));
    spectraCtrlPre(mouseCnt,:) = mean(SDCtrl.scorings{mouseCnt}.spectr(currNrEpPre,:),1);
    spectraCtrlPost(mouseCnt,:) = mean(SDCtrl.scorings{mouseCnt}.spectr(currNrEpPost,:),1);
end

SWACtrlPre = mean(spectraCtrlPre(:,SWARange),2);
SWACtrlPost = mean(spectraCtrlPost(:,SWARange),2);
SWAOptoPre = mean(spectraOptoPre(:,SWARange),2);
SWAOptoPost = mean(spectraOptoPost(:,SWARange),2);

figure;
toPlot = {[tFirstNROpto,tFirstNRCtrl]};
pVals =nan(size(toPlot));ax = cell(size(toPlot));
labels={'time to first NREM'};
for plotCnt = 1:numel(toPlot)
    p1=plot(1:2,toPlot{plotCnt},'o-k','MarkerFaceColor',[0.8,0.8,0.8]); hold on;
    for i = 1:numel(p1); p1(i).Color(4)=0.3; end
    boxplot(toPlot{plotCnt},'labels',{'opto-SD','ctrl-SD'},'notch','off');hold on;
    [~,pVals(plotCnt),CI,stats] = ttest(toPlot{plotCnt}(:,1),toPlot{plotCnt}(:,2));
    %     [pVals(plotCnt),~,stats] = signrank(toPlot{plotCnt}(:,1),toPlot{plotCnt}(:,2));
    
    % cosmetics
    h = findobj(gca,'Tag','Box');
    h2= findobj(gca,'Tag','Outliers');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',.4);
        h(j).Color = 'k'; % set box edge color
        h2(j).MarkerEdgeColor=[0.8 0.8 0.8];
    end
    ylabel(labels{plotCnt});
    set(gca,'TickDir','out');box off; title(labels{plotCnt});
    ax{plotCnt}=gca;
    line(get(gca,'XLim'),[0 0],'Color','k','LineStyle','--')
end
saveas(gcf,fullfile(pathFigs,'tAwake.fig'));
saveas(gcf,fullfile(pathFigs,'tAwake.png'));



% now correlate change in SWA with time awake
miceToUse = 1:nMice;
miceToUse(mNoFro)=[];
figure;
changeInSWACtrl = SWACtrlPost-SWACtrlPre;
changeInSWAOpto = SWAOptoPost-SWAOptoPre;
scatter(tFirstNROpto(miceToUse),changeInSWAOpto(miceToUse),'r','MarkerFaceColor','r'); hold on;
fitlineOpto = polyfit(tFirstNROpto(miceToUse),changeInSWAOpto(miceToUse),1);
xlabel('time awake'); ylabel('increase in SWA [uV^2]');

scatter(tFirstNRCtrl(miceToUse),changeInSWACtrl(miceToUse),'k','MarkerFaceColor','k'); hold on;
fitlineCtrl = polyfit(tFirstNRCtrl(miceToUse),changeInSWACtrl(miceToUse),1);
% plot(get(gca,'XTick'),fitlineCtrl(2)+get(gca,'XTick').*fitlineCtrl(1),'k');
plot(get(gca,'XTick'),fitlineOpto(2)+get(gca,'XTick').*fitlineOpto(1),'r');
xlabel('time awake'); ylabel('increase in SWA [uV^2]');

saveas(gcf,fullfile(pathFigs,'corrSWATAwake.fig'));
saveas(gcf,fullfile(pathFigs,'corrSWATAwake.png'));

%% now check whether sleep architecture is somehow altered after the animals fell asleep
% hypothesis: perhaps increased arousal leads to more disrupted sleep

% only look at bouts within 4h after sleep dep
relevantWBoutIdxs=cellfun(@(x)sum(and(x>240,x<480),2)==2,wakeBoutOnOffsets,'un',0);
wBoutDurCtrl=arrayfun(@(x) wBoutDuration{x,1}(relevantWBoutIdxs{x,1}),1:nMice,'un',0);
wBoutDurOpto=arrayfun(@(x) wBoutDuration{x,2}(relevantWBoutIdxs{x,2}),1:nMice,'un',0);
wBoutDurCtrl=cellfun(@log,wBoutDurCtrl,'un',0);
wBoutDurOpto=cellfun(@log,wBoutDurOpto,'un',0);

extremaWakeDur = [cell2mat(wBoutDurCtrl');cell2mat(wBoutDurOpto')];
extremaWakeDur = [min(extremaWakeDur),max(extremaWakeDur)];

relevantNRBoutIdxs=cellfun(@(x)sum(and(x>240,x<480),2)==2,NREMboutOnOffsets,'un',0);
NRBoutDurCtrl=arrayfun(@(x) nRBoutDuration{x,1}(relevantNRBoutIdxs{x,1}),1:nMice,'un',0);
NRBoutDurOpto=arrayfun(@(x) nRBoutDuration{x,2}(relevantNRBoutIdxs{x,2}),1:nMice,'un',0);
NRBoutDurCtrl=cellfun(@log,NRBoutDurCtrl,'un',0);
NRBoutDurOpto=cellfun(@log,NRBoutDurOpto,'un',0);

% histNRdurOpto =
extremaNRDur = [cell2mat(NRBoutDurCtrl');cell2mat(NRBoutDurOpto')];
extremaNRDur = [min(extremaNRDur),max(extremaNRDur)];


figure;
for mouseCnt = 1:nMice
    subplot(2,nMice,mouseCnt)
    edges = linspace(extremaWakeDur(1),extremaWakeDur(2),10);
    ctrl=histcounts((wBoutDurCtrl{mouseCnt}),edges,'Normalization','probability');
    opto=histcounts(((wBoutDurOpto{mouseCnt})),edges,'Normalization','probability');
    plot(edges(2:end),ctrl,'k');hold on;
    plot(edges(2:end),opto,'r');hold on;
    legend('Ctrl-SD','opto-SD'); xlabel('bout duration'); ylabel('probability'); title('w epochs')
    
    subplot(2,nMice,mouseCnt+nMice)
    edges = linspace(extremaNRDur(1),extremaNRDur(2),10);
    ctrl=histcounts((NRBoutDurCtrl{mouseCnt}),edges,'Normalization','probability');
    opto=histcounts(((NRBoutDurOpto{mouseCnt})),edges,'Normalization','probability');
    plot(edges(2:end),ctrl,'k');hold on;
    plot(edges(2:end),opto,'r');hold on;
    legend('Ctrl-SD','opto-SD'); xlabel('bout duration'); ylabel('probability'); title('nr epochs')
end

% get mean duration and number of bouts for wake and NREM
meanDurWCtrl=  cellfun(@nanmean,wBoutDurCtrl); % get mean duration
meanDurWOpto=  cellfun(@nanmean,wBoutDurOpto);% get mean duration
meanDurWakeCombined= [meanDurWCtrl;meanDurWOpto]';% combine vectors for plotting
nWakeEpCtrl = cellfun(@numel,wBoutDurCtrl); % number of bouts
nBriefWakeCtrl= cellfun(@(x)sum(x<0.25),wBoutDurCtrl);
nBriefWakeOpto= cellfun(@(x)sum(x<0.25),wBoutDurOpto);
nWakeEpOpto= cellfun(@numel,wBoutDurOpto);
nWakeCombined = [nWakeEpCtrl;nWakeEpOpto]';

% same for NR
meanDurNRCtrl=  cellfun(@nanmean,NRBoutDurCtrl);
meanDurNROpto=  cellfun(@nanmean,NRBoutDurOpto);
meanDurNRCombined= [meanDurNRCtrl;meanDurNROpto]';
nNREpCtrl = cellfun(@numel,NRBoutDurCtrl);
nNREpOpto= cellfun(@numel,NRBoutDurOpto);
nNRCombined = [nNREpCtrl;nNREpOpto]';


figure;
subplot(2,2,1)
plot(1:2,nWakeCombined,'o-','Color',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
xlim([0,3]); box off; hold on; ylabel('mean number of wake epochs');
errorbar(1:2,mean(nWakeCombined,1),std(nWakeCombined,[],1)./sqrt(nMice),'Color','k');

subplot(2,2,2)
plot(1:2,meanDurWakeCombined,'o-','Color',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
xlim([0,3]); box off; hold on; ylabel('mean duration wake epochs');
errorbar(1:2,mean(meanDurWakeCombined,1),std(meanDurWakeCombined,[],1)./sqrt(nMice),'Color','k');

subplot(2,2,3)
plot(1:2,nNRCombined,'o-','Color',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
xlim([0,3]); box off; hold on; ylabel('mean number of NREM epochs');
errorbar(1:2,mean(nNRCombined,1),std(nNRCombined,[],1)./sqrt(nMice),'Color','k');

subplot(2,2,4)
plot(1:2,meanDurNRCombined,'o-','Color',[0.8 0.8 0.8],'MarkerFaceColor',[0.8 0.8 0.8]);
xlim([0,3]); box off; hold on; ylabel('mean duration NREM epochs');
errorbar(1:2,mean(meanDurNRCombined,1),std(meanDurNRCombined,[],1)./sqrt(nMice),'Color','k');