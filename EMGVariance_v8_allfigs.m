%V3 sepparating the file loading from extracting relevant stuff  - this should allow rapid iterations across parameters
%V5: calculate emg threshold based on the distribution of emg variance between sleep and wake
%v6: get rid of arousal "intensity" and plot spindle correlation instead
%v7: sepparating figures into different functions: (i) correlation of spectra with arousal delay (ii) high vs low sleep pressure experiments
%%
% path to files etc
clear all

paths.stimFiles = 'F:\Tomoko_OptoStim\STIMs';
paths.sleepScoring = 'F:\Tomoko_OptoStim\outputVSchr';
paths.rawSignals = 'F:\Tomoko_OptoStim\OutputSignals';
paths.extractedSignals ='F:\Tomoko_OptoStim\EMGAnalysis_matfiles';
paths.code = 'C:\Users\Martin Kahn\Desktop\Scripts_and_Code\Matlab';

% set parameters
settings.sR = 256;
settings.anmPrefix = 'GDCh';
% settings.anmPrefix = 'GFP';
settings.tPre = 4; % extract EMG x seconds before stim
settings.tPost = 120; % extract EMG x seconds after stim
settings.EMGthresh = 5; % how many times SD is an arousal
settings.minDur = 0.25; % how many secs for an arousal
settings.minDur = settings.minDur*4; % convert to samples (i.e. 2 seconds min dur with 250ms per variance bin = 8 samples)
settings.nEpochsForSWA = 2;
settings.binSize_EMGVar = 0.25; % over how many s will emg variance be calculated
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 240 228 66;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue,yellow
npgCMap = [npgCMap;npgCMap];

% animal numbers
mouseNr.TenHz24h = [5,15:21];
mouseNr.LPOBl = [15:21];
mouseNr.TenHz24hGFP = [1:4,6:8];
mouseNr.TenHz24hWithExclusion = [5,15:21];
mouseNr.TenHzNonLPO = [1,6,7,8,9,12];
mouseNr.TenHz24hSameMiceAs5Hz = 15:21;
mouseNr.TenHz24hBaseline = [5,15:19,21];
mouseNr.twentyHz24h = 18:21; 
mouseNr.fiveHz24h = 15:21;
mouseNr.fiveHz24hWithExclusion = [15:19,21];
mouseNr.twoHz24h = 15:21;
mouseNr.oneHz24h = 15:21;
mouseNr.HSP = [5,15:21]; % high sleep pressure
mouseNr.LSP = [5,15:21]; % low sleep pressure
mouseNr.sedation = [5,16:19,21]; % low sleep pressure
mouseNr.tenHzSameAsSedation = [5,16:19,21];
mouseNr.SDOpto_LPO=[5,16,17,18,19,20,21];
mouseNr.SDctrl_LPO = [5,16,17,18,19,20,21];
mouseNr.SDOpto_nLPO=[1,6,7,8,9,12];
mouseNr.SDctrl_nLPO = [1,6,7,8,9,12];


dates.TenHz24h = [030418 191218 191218 191218 160819 160819 180919 160819];
dates.LPOBl = [181218 181218 181218 100819 100819 070819 100819];

dates.TenHz24hGFP = [260218,120418,260618,040718,061018,061018,070119];
dates.oneHz24h = [241218 241218 241218 190819 190819 190819 190819];
dates.twoHz24h = [ 251218 251218 251218 200819 200819 200819 200819];
dates.fiveHz24h = [261218 261218 261218	210819 210819 210819 210819];
dates.twentyHz24h = [140919 140919 140919 140919];

dates.fiveHz24hWithExclusion = [261218 261218 261218 200819 200819 200819];
dates.TenHz24hWithExclusion = [030418 191218 191218 191218 160819 160819 180919 160819];
dates.TenHz24hSameMiceAs5Hz = [191218 191218 191218 160819 160819 180919 160819];
dates.TenHzNonLPO = [140218,120418,130518,150618,150618,090618];
dates.twentyHz24h = [140919,140919,140919,140919];

dates.TenHz24hBaseline = [020418 181218 181218 181218 100819 100819 100819];
dates.HSP = [310318 211218 211218 211218 060919 060919 060919 060919];
dates.LSP = [300318 201218 201218 201218 040919 040919 040919 040919];
dates.sedation = [010518 150119 150119 220819 220819 220819];
dates.tenHzSameAsSedation =[030418 191218 191218 160819 160819 160819];
dates.SDOpto_LPO = [040518 090119 090119 080819 110819 080819 110819];  % anm 17 may not have frontal EEG 
dates.SDctrl_LPO = [060518 110119 110119 110819 080819 110819 080819]; 
dates.SDOpto_nLPO = [040518 300418 030718 130618 130618 130618];  % anm 17 may not have frontal EEG 
dates.SDctrl_nLPO = [060518 020518 010718 110618 110618 110618]; 

% choose which condition to use
currCond = 'fiveHz24h';
%% load data and filter EMG
addpath(genpath(paths.code));
nMice = numel(mouseNr.(currCond));
[EMGs, scorings,scorings2, stims,EMGVarRatios,EMGVarBinning,fros,occs] = deal(cell(nMice,1));
for mouseCnt = 1:nMice
    mouseCnt
    % get all filenames matching the current date and mouse in all folders
    currMouse = [settings.anmPrefix,num2str(mouseNr.(currCond)(mouseCnt))];
    currDate = dates.(currCond)(mouseCnt);
    if isnan(currDate); error('trying to access nonexistent recording'); end
    fileNames = structfun(@(x)dir(fullfile(x,[currMouse,'*',num2str(currDate,'%06d'),'*'])),paths,'un',0); % search all relevant folders for mousename & date
    
    % stim (sR: 256) - if this is a baseline recording skip this
    if or(contains(currCond,'Baseline'),contains(currCond,'SDctrl'))
        currStim = [];
    elseif contains(currCond,'LPOBl') % if this is LPO baseline, just use stim times from stim day 
        currStimDate = dates.TenHz24hSameMiceAs5Hz(mouseCnt);
        stimFile = dir(fullfile(paths.stimFiles,[currMouse,'*',num2str(currStimDate,'%06d'),'*']));
        currStim = load(fullfile(paths.stimFiles,stimFile.name),'startend'); % import start and end times of stim
        currStim = currStim.startend;
    else
        currStim = load(fullfile(paths.stimFiles,fileNames.stimFiles.name),'startend'); % import start and end times of stim
        currStim = currStim.startend;
    end
    %sleepScoring
%     if contains(currCond,'SD')
%         epochSO = load(fullfile(paths.sleepScoring,fileNames.sleepScoring(1).name));
%         tmpIdx = 2; % monkeypatch: in the sleep enhancement condition there is a file called epochSO which stores when sleep onset is
%     else
%         tmpIdx = 1; % monkeypatch: imn all other conditions, this doesn't exists and file 1 is simply derivation 1 and so on
%     end
    tmpIdx = 1;   
    currScoring = load(fullfile(paths.sleepScoring,fileNames.sleepScoring(tmpIdx).name));
    currScoring.derivation = fileNames.sleepScoring(tmpIdx).name(end-6:end-4);
    tmp = load(fullfile(paths.sleepScoring,fileNames.sleepScoring(tmpIdx+1).name));
    currScoring2.spectr = tmp.spectr;
    currScoring2.derivation = fileNames.sleepScoring(tmpIdx+1).name(end-6:end-4);
    
    % load emg signal (sR: 256)
    currEMGIdx = arrayfun(@ (x)contains(fileNames.rawSignals(x).name,'emg'),1:numel(fileNames.rawSignals)); % find emg (as opposed to EEG etc)
    currEMGIdx = fileNames.rawSignals(find(currEMGIdx,1)).name;
    currEMG = load(fullfile(paths.rawSignals,currEMGIdx)); % load emg
    currEMG = currEMG.resampled_sig;
    currFroIdx = arrayfun(@ (x)contains(fileNames.rawSignals(x).name,'fr'),1:numel(fileNames.rawSignals)); % find emg (as opposed to EEG etc)
    currFroIdx = fileNames.rawSignals(find(currFroIdx,1)).name;
    currFro = load(fullfile(paths.rawSignals,currFroIdx)); % load emg
    currFro = currFro.resampled_sig;
    currOccIdx = arrayfun(@ (x)contains(fileNames.rawSignals(x).name,'oc'),1:numel(fileNames.rawSignals)); % find emg (as opposed to EEG etc)
    currOccIdx = fileNames.rawSignals(find(currOccIdx,1)).name;
    currOcc = load(fullfile(paths.rawSignals,currOccIdx)); % load emg
    currOcc = currOcc.resampled_sig;

    
    % filter EMG using chebychev bandpass filter
    sR =settings.sR;
    p1=5; p2=100; s1=4; s2=120;
    Wp=[p1 p2]/(sR/2); Ws=[s1 s2]/(sR/2); Rp=3; Rs=30; [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
    [bb1,aa1]=cheby2(n,Rs,Wn);
    currEMG=filtfilt(bb1,aa1,currEMG);
    
    %get distribution of EMG variance across sleep epochs and wake epochs
    binSize = sR/(1/settings.binSize_EMGVar);
    nREpochs = currScoring.nr;
    nR_EMG = cell2mat(arrayfun(@(x) currEMG(x*sR+1:x*sR+sR*4),nREpochs*4 -4,'un',0)); % trial-wise opto-evoked EMG
    nR_EMGVar = squeeze(var(reshape(nR_EMG,size(nR_EMG,1),binSize,size(nR_EMG,2)/binSize),[],2));
    nR_EMGVar = mean(nR_EMGVar,2);
    wEpochs = currScoring.w(1:end-1); % to avoid running over the edge at the end
    w_EMG = cell2mat(arrayfun(@(x) currEMG(x*sR+1:x*sR+sR*4),wEpochs*4 -4,'un',0)); % trial-wise opto-evoked EMG
    w_EMGVar = squeeze(var(reshape(w_EMG,size(w_EMG,1),binSize,size(w_EMG,2)/binSize),[],2));
    w_EMGVar = mean(w_EMGVar,2);

    %calculate probability (relative to all w/nrem trials) to display a given EMG variance
    [wCounts,EMGVarBins]=histcounts(log10(w_EMGVar),500);
    [nrCounts,~]=histcounts(log10(nR_EMGVar),EMGVarBins);
    wCounts = 100.*wCounts./(numel(w_EMGVar)+numel(nR_EMGVar)); % normalise to total amount of trials
    nrCounts= 100.*nrCounts./(numel(w_EMGVar)+numel(nR_EMGVar));
    wCounts=smooth(wCounts,5);nrCounts=smooth(nrCounts,5);
    EMGVarRatio=wCounts./nrCounts;
    EMGVarRatio =smooth(EMGVarRatio,10);
    
    % show above-calculated distribution of EMG variance across vigilance states
    figure;
    plot(EMGVarBins(1:end-1),wCounts,'b'); hold on; plot(EMGVarBins(1:end-1),nrCounts,'Color','g');%plot(edges(1:end-1),rCounts,'r');
    xlabel('log EMG variance');ylabel('% total trials');  
    %define threshold as ratio = 2
    %     line([ratioBasedThresh ratioBasedThresh],get(gca,'YLim'),'Color','r','LineStyle','--')
    yyaxis right
    plot(EMGVarBins(1:end-1),EMGVarRatio,'-','Color','k')
    title(['GDCh',num2str(mouseNr.(currCond)(mouseCnt))]); % title
    box off; set(gca,'TickDir','out'); ylabel('ratio NREM/wake'); ax = gca; ax.YColor = npgCMap(3,:);
    
    % store data from this mouse
    EMGs{mouseCnt} = currEMG;
    fros{mouseCnt} = currFro;
    occs{mouseCnt} = currOcc;

    stims{mouseCnt} = currStim;
    scorings{mouseCnt} = currScoring;
    scorings2{mouseCnt} = currScoring2;
    EMGVarRatios{mouseCnt} = EMGVarRatio;
    EMGVarBinning{mouseCnt} = EMGVarBins;
end
%% save result
save(fullfile(paths.extractedSignals,currCond),'EMGs','fros','occs','stims','scorings','EMGVarRatios','EMGVarBinning','settings','scorings2');
%% get wavelets or load them 
pathFig = 'F:\Tomoko_OptoStim\Fgiures_wavelets';
settings.EMGThresholds= [2]'; % EMG threshold
settings.nEpochsSWA= (1)'; % EMG threshold
settings.minDur= 4; % EMG threshold
settings.tPost = 60;  % this will determine extracted time (sec) after arousal, also influences min arousal duration criterion for plotted episodes
settings.tPre = 20;  % this will determine extracted time (sec) after arousal, also influences min arousal duration criterion for plotted episodes

output=extract_features_from_opto_evoked_arousals_wavelets(EMGs,stims,scorings,EMGVarRatios,EMGVarBinning,settings,fros,occs); % runs the script

% save(fullfile(paths.extractedSignals,[currCond,'wavelets']),'output','-v7.3');
load(fullfile(paths.extractedSignals,[currCond,'wavelets']),'output');

%% plot mean wavelets for each mouse
tInSecs = -output.nEpsForWavelets*4+1/256:1/256:settings.tPost; % time in seconds
freqs = output.waveletFrqs{1};
cMap = 'parula';
labelling = {'light','dark'};
currWavs = struct();
for mouseCnt = 1:nMice
    spontStimTInH = output.stimTimesInS_spont{mouseCnt}./(settings.sR*3600);% get stimTimes in hours
    evokedStimTInH = output.stimTimesInS_evoked{mouseCnt}./(settings.sR*3600);% get stimTimes in hours
    
    % get mean wavelets sepparately for light and dark periods
    currSpontStims_light = spontStimTInH<12; % select light or dark phase stims
    currEvokedStims_light =evokedStimTInH<12; % select light or dark phase stims
    currSpontStims_dark = spontStimTInH<12; % select light or dark phase stims
    currEvokedStims_dark =evokedStimTInH<12; % select light or dark phase stims
    
    currWavs.OptoFro_L = mean(output.optoFroNR_wav{mouseCnt}(:,:,evokedStimTInH<12),3);
    currWavs.OptoOcc_L = mean(output.optoOccNR_wav{mouseCnt}(:,:,evokedStimTInH<12),3);
    currWavs.SpontFro_L = mean(output.spontFroNR_wav{mouseCnt}(:,:,spontStimTInH<12),3);
    currWavs.SpontOcc_L = mean(output.spontOccNR_wav{mouseCnt}(:,:,spontStimTInH<12),3);
    
    currWavs.OptoFro_D = mean(output.optoFroNR_wav{mouseCnt}(:,:,evokedStimTInH>12),3);
    currWavs.OptoOcc_D = mean(output.optoOccNR_wav{mouseCnt}(:,:,evokedStimTInH>12),3);
    currWavs.SpontFro_D = mean(output.spontFroNR_wav{mouseCnt}(:,:,spontStimTInH>12),3);
    currWavs.SpontOcc_D = mean(output.spontOccNR_wav{mouseCnt}(:,:,spontStimTInH>12),3);
    
    toPlot = {'SpontFro_L','SpontFro_D','SpontOcc_L','SpontOcc_D','OptoFro_L','OptoFro_D','OptoOcc_L','OptoOcc_D'};
    % determine plot limits to have the same color ranges in different plots
    minVals = structfun(@(x) min(x(:)),currWavs);
    minVals = [min(minVals([1,3,5,7])),min(minVals(2:2:8))];
    maxVals = structfun(@(x) max(x(:)),currWavs);
    maxVals = [max(maxVals([1,3,5,7])),max(maxVals(2:2:8))];
    
    
    figure;
    for pCnt = 1:numel(toPlot)
        subplot(2,4,pCnt);
        surface(tInSecs,freqs',currWavs.(toPlot{pCnt})); shading flat
        %         surface(currWavs.(toPlot{pCnt})); shading flat
        
        ax = gca; ax.YDir = 'normal';title(toPlot{pCnt},'Interpreter','none'); xlabel('time [sec]'); ylabel('frequency [Hz]');
        colormap(cMap);colorbar();
        if contains(toPlot{pCnt},'Fro')
            minIdx =1;
        else
            minIdx =2;
        end
        caxis([minVals(minIdx),maxVals(minIdx)]);colormap(cMap);colorbar()
    end
    set(gcf,'Position',[1933 296 1080 790]);title(labelling{lightDark});
    
    pause(1)
    saveas(gcf,fullfile(pathFig,[currCond,' mouse',num2str(mouseCnt),labelling{lightDark},'.png']));
    saveas(gcf,fullfile(pathFig,[currCond,' mouse',num2str(mouseCnt),'.fig']));
    
    close all
    pause(1)
end


%% plot mean opto-evoked EMGs in all states
settings.EMGThresholds= [2]'; % EMG threshold
settings.nEpochsSWA= (1)'; % EMG threshold
tomokoColors =[0 114 178;230 159 0; 213 94 0 ]./255;
tPrePost = {[8,15],[60,120]};
tLineWidth =[1.2,0.9];
settings.minDur = 4;
mouseForFig = 8;
for plotMeTwice = 1:2
    settings.tPre= tPrePost{plotMeTwice}(1)'; % EMG threshold
    settings.tPost= tPrePost{plotMeTwice}(2)'; % EMG threshold
    
    output=extract_features_from_opto_evoked_arousals(EMGs,stims,scorings,EMGVarRatios,EMGVarBinning,settings); % runs the script
    [dEMGVar_W,dEMGVar_R,dEMGVar_NR]=deal(cell(nMice,1));
    figure(plotMeTwice);
    for mouseCnt = 1:nMice
        currMouse = [settings.anmPrefix,num2str(mouseNr.(currCond)(mouseCnt))];
        
        
        %wake
        subplot(4,2,mouseCnt)
        timeInS = -settings.tPre:0.25:settings.tPost-0.25;
        %wake
        if plotMeTwice==2
            currEMG = output.optoEvokedEMG_W_corrected{mouseCnt};
            currEMG = cell2mat(arrayfun(@(x)smooth(currEMG(x,:),12),1:size(currEMG,1),'un',0))';
        else
           currEMG = output.optoEvokedEMG_W{mouseCnt};
           currEMG = smoothdata(currEMG,2,'movmean',25);
        end
        meanVar= median(currEMG,1);
        VarPtiles = prctile(currEMG,[25,75],1);
        p1=plot(timeInS,meanVar,'Color',tomokoColors(1,:),'lineWidth',tLineWidth(plotMeTwice)); hold on; % plot mean evoked variance
        pa1=patch([timeInS,fliplr(timeInS)],[VarPtiles(1,:),fliplr(VarPtiles(2,:))],tomokoColors(1,:),'FaceAlpha',0.3,'EdgeColor','none','FaceColor',tomokoColors(1,:)); % add 25/75%percentiles as shading
        
        %NREM
        currEMG = output.optoEvokedEMG_NR{mouseCnt};
        if plotMeTwice==2
            currEMG = cell2mat(arrayfun(@(x)smooth(currEMG(x,:),12),1:size(currEMG,1),'un',0))';
        end
        meanVar= median(currEMG,1);
        VarPtiles = prctile(currEMG,[25,75],1);
        p2=plot(timeInS,meanVar,'Color',tomokoColors(2,:),'lineWidth',tLineWidth(plotMeTwice)); hold on; % plot mean evoked variance
        patch([timeInS,fliplr(timeInS)],[VarPtiles(1,:),fliplr(VarPtiles(2,:))],tomokoColors(2,:),'FaceAlpha',0.3,'EdgeColor','none','FaceColor',tomokoColors(2,:)); % add 25/75%percentiles as shading
        dEMGVar_NR{mouseCnt} = mean(currEMG(:,settings.tPre/0.25:settings.tPre/0.25+16),2)./mean(currEMG(:,1:settings.tPre/0.25),2);
        
        %REM
        currEMG = output.optoEvokedEMG_R{mouseCnt};
        if plotMeTwice==2
            currEMG = cell2mat(arrayfun(@(x)smooth(currEMG(x,:),12),1:size(currEMG,1),'un',0))';
        end
        meanVar= median(currEMG,1);
        VarPtiles = prctile(currEMG,[25,75],1);
        p3=plot(timeInS,meanVar,'Color',tomokoColors(3,:),'lineWidth',tLineWidth(plotMeTwice)); hold on; % plot mean evoked variance
        patch([timeInS,fliplr(timeInS)],[VarPtiles(1,:),fliplr(VarPtiles(2,:))],tomokoColors(3,:),'FaceAlpha',0.3,'EdgeColor','none','FaceColor',tomokoColors(3,:)); % add 25/75%percentiles as shading
        dEMGVar_R{mouseCnt} = mean(currEMG(:,settings.tPre/0.25:settings.tPre/0.25+16),2)./mean(currEMG(:,1:settings.tPre/0.25),2);
        

        dEMGVar_W{mouseCnt} = mean(currEMG(:,settings.tPre/0.25:settings.tPre/0.25+16),2)./mean(currEMG(:,1:settings.tPre/0.25),2);
        xlabel('time to stim [s]');ax = gca; ylabel('log EMG var'); box off; set(gca,'TickDir','out','YScale','log');%xlim([-settings.tPre,15]); %cosmetics
        if mouseCnt ==1
            legend([p1,p2,p3,pa1],'median w','median NR','median_REM','25/75 %tiles','AutoUpdate','off','Interpreter','none'); legend boxoff
        end
        line([0,0],[ax.YLim(1),6*10^4],'Color','k','lineStyle','--'); title(currMouse);
    end
    linkaxes
end
figure(1);
saveas(gcf,fullfile('F:\Tomoko_OptoStim\figs_10Hz','opto_evoked_EMG_var'))
saveas(gcf,fullfile('F:\Tomoko_OptoStim\figs_10Hz','opto_evoked_EMG_var.png'))
figure(2);
saveas(gcf,fullfile('F:\Tomoko_OptoStim\figs_10Hz','opto_triggered_EMG_var_zoomOut'))
saveas(gcf,fullfile('F:\Tomoko_OptoStim\figs_10Hz','opto_triggered_EMG_var_zoomOut.png'))

%plot change in EMG variance 
%% calculate correlations without plotting
% for each set of parameters (i.e. 1 emg threshold and one SWA duration, we'll do 4 corelations
% extract features: arousal delay, preceeding power in delta, theta, sigma, time to last waking
EMGThresholds= [1,2,3,4,6,8]'; % EMG threshold
nEpochsSWA= (1:2:10)'; % EMG threshold
settings.EMGThresholds = EMGThresholds;
settings.nEpochsSWA = nEpochsSWA;
settings.minDur = 4;
settings.tPre= 4; 
settings.tPost= 15; % EMG threshold
output=extract_features_from_opto_evoked_arousals(EMGs,stims,scorings,EMGVarRatios,EMGVarBinning,settings); % runs the script

[corrCoefs,corrPVals] = deal(nan(4,numel(EMGThresholds),numel(nEpochsSWA)));
clearvars xData yData
iteration = 0;
for EMGCnt = 1:numel(EMGThresholds)
    EMGthresh = EMGThresholds(EMGCnt);
    for SWACnt = 1:numel(nEpochsSWA)
        nEpochsForSWA = nEpochsSWA(SWACnt);
        iteration = iteration + 1;
        for mouseCnt = 1:nMice
            %load variable and get rid of nans (no arousal detected)
            tArous_NR = output.tArousal_NR{iteration,mouseCnt}';
            preSWA=output.preStimSWA{iteration,mouseCnt};
            preTheta=output.prestimTheta{iteration,mouseCnt};
            tToLastWake = output.tSinceWake{iteration,mouseCnt};
            preSigma=output.prestimSigma{iteration,mouseCnt} ;
            tToNROnset = output.NRBoutDur{iteration,mouseCnt};

            
            % flag nan values
            noArousals = isnan(tArous_NR);
            noSWA =  isnan(preSWA);
            noPriorWake = isnan(tToLastWake);
            noSigma = isnan(preSigma);
            noTheta = isnan(preTheta);
            
            % plot1: SWA v delay to arousal
            xData{1}{mouseCnt}=preSWA((noArousals+noSWA)==0); % zscore SWA and remove trials without arousal
            yData{1}{mouseCnt}=tArous_NR((noArousals+noSWA)==0);% zscore arousal time and remove trials without arousal
            labels{1}={'SWA vs arousal delay'}; % store X and Y label string
            
            % plot sigma v delay to arousal
            xData{2}{mouseCnt}=preSigma((noArousals+noSigma)==0);
            yData{2}{mouseCnt}=tArous_NR((noArousals+noSigma)==0);
            labels{2}={'sigma vs arousal intensity'}; % store X and Y label string
            
            % plot theta v dfelay to arousal
            xData{3}{mouseCnt}=preTheta((noArousals+noTheta)==0);
            yData{3}{mouseCnt}=tArous_NR((noArousals+noTheta)==0);
            labels{3}={'theta vs arousal delay '}; % store X and Y label string
            
            % plot tSinceWaking v dfelay to arousal
            xData{4}{mouseCnt}=tToLastWake((noArousals+noPriorWake)==0);
            yData{4}{mouseCnt}=tArous_NR((noArousals+noPriorWake)==0);
            labels{4}={'t to last wake vs arousal delay '}; % store X and Y label string
            
            % plot tSinceWaking v dfelay to arousal
            xData{5}{mouseCnt}=tToNROnset((noArousals+noPriorWake)==0);
            yData{5}{mouseCnt}=tArous_NR((noArousals+noPriorWake)==0);
            labels{5}={'NREM bout dur (excluding mt,REM,w)'}; % store X and Y label string  
            
        end
        
        % remove outliers ugly written,should move to loop below sorry 
        getNotOutliers = @(x) ~or(x>median(x)+1.5*iqr(x),x<median(x)-1.5*iqr(x)); % function to flag non-outlier data points
        notOutX= cellfun(@(x) cellfun(@(x)getNotOutliers(x),x,'un',0),xData,'un',0); % apply the above function to every cell
        notOutY= cellfun(@(x) cellfun(@(x)getNotOutliers(x),x,'un',0),yData,'un',0);
        xData = arrayfun(@(y) arrayfun(@(x) xData{y}{x}(and(notOutX{y}{x},notOutY{y}{x})),1:numel(xData{y}),'un',0),1:numel(xData),'un',0);
        yData = arrayfun(@(y) arrayfun(@(x) yData{y}{x}(and(notOutX{y}{x},notOutY{y}{x})),1:numel(yData{y}),'un',0),1:numel(yData),'un',0);
        
        % zscore all the above-created x and y data
        xData_Zsc = cellfun(@(x) cellfun(@zscore,x,'un',0),xData,'un',0);
        yData_Zsc = cellfun(@(x) cellfun(@zscore,x,'un',0),yData,'un',0);
        
        % now do overall correlations for these parameters
        for corrCnt = 1:numel(yData)
            overallX = cell2mat(xData_Zsc{corrCnt}');
            overallY = cell2mat(yData_Zsc{corrCnt}');
            [pearson,pVal] = corrcoef(overallX,overallY);
            corrCoefs(corrCnt,EMGCnt,SWACnt)=pearson(2);
            corrPVals(corrCnt,EMGCnt,SWACnt)=pVal(2);  
        end
    end
end

%% plot dependency of correlation coefficient on parameters
f1=figure;
f2 = figure;

% do correlations
for corrCnt = 1:numel(yData)
    
    % plot correlations
    figure(f1)
    subplot(2,3,corrCnt)
    for SWACnt = 1:numel(nEpochsSWA)
        plot(EMGThresholds,squeeze(corrCoefs(corrCnt,:,SWACnt)),'Color',npgCMap(SWACnt,:)); hold on;
    end
    % cosmetics
    if corrCnt <3
        legend([num2str(nEpochsSWA),repmat(' epochs',numel(nEpochsSWA),1)]); legend boxoff;
    end
    title(labels{corrCnt}); ylabel('pearsons R'); xlabel('EMG detections threshold [*SD]');
    box off; set(gca,'TickDir','out');linkaxes;
    
    % plot pvalue
    figure(f2)
    subplot(2,3,corrCnt)
    for SWACnt = 1:numel(nEpochsSWA)
        plot(EMGThresholds,squeeze(corrPVals(corrCnt,:,SWACnt)),'-o','Color',npgCMap(SWACnt,:),...
            'MarkerFaceColor',npgCMap(SWACnt,:),'MarkerSize',3); hold on;
    end
    %cosmetics
    if corrCnt <3
        legend([num2str(nEpochsSWA),repmat(' epochs',numel(nEpochsSWA),1)]); legend boxoff;
    end
    title(labels{corrCnt}); ylabel('p-value'); xlabel('EMG detections threshold [*SD]')
    line(get(gca,'XLim'),[0.05 0.05],'Color','k','LineStyle','--');
    box off; set(gca,'TickDir','out');linkaxes;
    
end

%% divide SWA into percentiles and plot relationship of arousal delay
% for each set of parameters (i.e. 1 emg threshold and one SWA duration, we'll do 4 corelations
[ttestPVals] = deal(nan(4,numel(EMGThresholds),numel(nEpochsSWA)));
clearvars xData yData
figure
iteration = 0;
for EMGCnt = 1:numel(EMGThresholds)
    EMGthresh = EMGThresholds(EMGCnt);
    for SWACnt = 1:numel(nEpochsSWA)
        nEpochsForSWA = nEpochsSWA(SWACnt);
        iteration = iteration + 1;
        for mouseCnt = 1:nMice
            %load variable and get rid of nans (no arousal detected)
            tArous_NR = output.tArousal_NR{iteration,mouseCnt}';
            preSWA=output.preStimSWA{iteration,mouseCnt};
            tToLastWake = output.tSinceWake{iteration,mouseCnt};
            preSigma=output.prestimSigma{iteration,mouseCnt} ;
            preTheta=output.prestimTheta{iteration,mouseCnt} ;
            
            measures = {preSWA,preSigma,preTheta,tToLastWake}; 
            for measureCnt = 1:4
                
                currMeasure = measures{measureCnt};
                currTArousal = tArous_NR;
                
                % flag&eliminate nan values
                nanVals = or(isnan(currMeasure),isnan(currTArousal));
                currMeasure(nanVals)=[];currTArousal(nanVals)=[];
                
                % get percentiles
                currPrctIdxs = prctile(currMeasure,0:20:100);
                currPrctIdxs = discretize(currMeasure,currPrctIdxs);
                
                % plot1: SWA v delay to arousal
                xData(measureCnt,mouseCnt,:)= arrayfun(@(x) mean(currMeasure(currPrctIdxs==x)),1:max(currPrctIdxs)); % zscore SWA and remove trials without arousal
                yData(measureCnt,mouseCnt,:)=arrayfun(@(x) mean(currTArousal(currPrctIdxs==x)),1:max(currPrctIdxs));% zscore arousal time and remove trials without arousal
                labels{1}='SWA'; % store X and Y label string
                labels{2}='sigma'; % store X and Y label string
                labels{3}='theta'; % store X and Y label string
                labels{4}='last'; % store X and Y label string
            end
        end
        % now do overall correlations for these parameters
        
        [~,pVal] = ttest(squeeze(yData(:,:,1))',squeeze(yData(:,:,size(yData,3)))');
        ttestPVals(:,EMGCnt,SWACnt)=pVal;
        mean1stPct = squeeze(mean(yData,2));
        %             meanlastPct = squeeze(mean(yData(:,:,end),2));
        %             mean1stPct = squeeze(mean(yData(:,:,1),2));
        semPrctile = squeeze(std(yData,[],2))./sqrt(size(yData,2));
        subplot(numel(EMGThresholds),numel(nEpochsSWA),iteration)
        b=bar(mean1stPct(:,[1,end])); hold on;
        e=errorbar([(1:numel(measures))-0.125,(1:numel(measures))+0.125],[mean1stPct(:,1);mean1stPct(:,end)],[semPrctile(:,1);semPrctile(:,end)],'k');
        e.LineStyle = 'none'; ax = gca; ax.XTickLabel = labels;
        if iteration == 1; legend('lowest %tile', 'highest %tile');end
    end
end

% plot dependency of ttest results on parameters 
figure; 
for comparison = 1:4
    subplot(2,2,comparison)
    plot(EMGThresholds,squeeze(ttestPVals(comparison,:,:))); 
    xlabel('EMG threshold'); ylabel('p value'); 
    legend([num2str(nEpochsSWA),repmat(' epochs',numel(nEpochsSWA),1)]); legend boxoff;
    line([EMGThresholds(1),EMGThresholds(end)],[0.05, 0.05],'LineStyle','--','Color','k')
end
linkaxes
%% calculate and plot correlations for each mouse sepparately and across mice (WITHOUT ITERATIONS)
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 240 228 66;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue,yellow
npgCMap = vertcat(npgCMap,npgCMap);
EMGSelect=4; % EMG threshold
SWASelect= 3; % SWA threshold
selectedIteration = (EMGSelect-1)*numel(nEpochsSWA) + SWASelect;
figs{1}=figure; figs{2}=figure;
clearvars xData yData
corrByMouse= nan(3,nMice,2);
for mouseCnt = 1:nMice
    
    %load variable and get rid of nans (no arousal detected)
    tArous_NR = output.tArousal_NR{selectedIteration,mouseCnt}';
%     preSWA=prestimTheta{selectedIteration,mouseCnt};
    preSWA=output.preStimSWA{selectedIteration,mouseCnt};
    preTheta=output.prestimTheta{selectedIteration,mouseCnt};

    tToLastWake = output.tSinceWake{selectedIteration,mouseCnt};
    preSigma=output.prestimSigma{selectedIteration,mouseCnt} ;
    
    % flag nan values
    noArousals = isnan(tArous_NR);
    noSWA =  isnan(preSWA);
    noPriorWake = isnan(tToLastWake);
    noSigma = isnan(preSigma);
    
    % plot1: SWA v delay to arousal
    xData{1}{mouseCnt}=preSWA((noArousals+noSWA)==0); % zscore SWA and remove trials without arousal
    yData{1}{mouseCnt}=tArous_NR((noArousals+noSWA)==0);% zscore arousal time and remove trials without arousal
    xLab{1}={'SWA'}; yLab{1}={'arousal delay'}; % store X and Y label string
    
    % plot SWA v intensity of arousal
    xData{2}{mouseCnt}=preSigma((noArousals+noSigma)==0);
    yData{2}{mouseCnt}=tArous_NR((noArousals+noSigma)==0);
    xLab{2}={'sigma'}; yLab{2}={'arousal delay'}; % store X and Y label string
    
    % plot tSinceWaking v dfelay to arousal
    xData{3}{mouseCnt}=tToLastWake((noArousals+noPriorWake)==0);
    yData{3}{mouseCnt}=tArous_NR((noArousals+noPriorWake)==0);
    xLab{3}={'t to last wake'}; yLab{3}={'arousal delay '}; % store X and Y label string
    
    % remove outliers ugly written,should move to loop below sorry
    getNotOutliers = @(x) ~or(x>median(x)+1.5*iqr(x),x<median(x)-1.5*iqr(x)); % function to flag non-outlier data points
    notOutX= cellfun(@(x) cellfun(@(x)getNotOutliers(x),x,'un',0),xData,'un',0); % apply the above function to every cell
    notOutY= cellfun(@(x) cellfun(@(x)getNotOutliers(x),x,'un',0),yData,'un',0);
    xData = arrayfun(@(y) arrayfun(@(x) xData{y}{x}(and(notOutX{y}{x},notOutY{y}{x})),1:numel(xData{y}),'un',0),1:numel(xData),'un',0);
    yData = arrayfun(@(y) arrayfun(@(x) yData{y}{x}(and(notOutX{y}{x},notOutY{y}{x})),1:numel(yData{y}),'un',0),1:numel(yData),'un',0);
    
    % zscore all the above-created x and y data
    xData_Zsc = cellfun(@(x) cellfun(@zscore,x,'un',0),xData,'un',0);
    yData_Zsc = cellfun(@(x) cellfun(@zscore,x,'un',0),yData,'un',0);
    
    % make the plots, specifically, make 2 figures, one with zscoring and one without
    for zscoring = 1:2
        % use either zscored or raw data
        if zscoring==1
            currXData = xData; currYData = yData;
        else
            currXData = xData_Zsc; currYData = yData_Zsc;
        end
        
        figure(figs{zscoring});
        for plotCnt = 1:numel(xData)
            subplot(2,2,plotCnt); hold on;
            p1=plot(currXData{plotCnt}{mouseCnt},currYData{plotCnt}{mouseCnt},'o','MarkerFaceColor',npgCMap(mouseCnt,:),'MarkerEdgeColor',npgCMap(mouseCnt,:),'MarkerSize',3); hold on;
            currCorr = corrcoef(currXData{plotCnt}{mouseCnt},currYData{plotCnt}{mouseCnt});
            corrByMouse(plotCnt,mouseCnt,zscoring) = currCorr(2);
            if mouseCnt > 7; p1.MarkerFaceColor = 'none';end % the colormap just has 7 entries, recycle if neccessary
            xlabel(xLab{plotCnt}); ylabel(yLab{plotCnt}); set(gca,'TickDir','out');
            % once last mouse was plotted, fir line across all points and for each mouse sepparately
            if mouseCnt == nMice
                overallX = cell2mat(currXData{plotCnt}'); overallY = cell2mat(currYData{plotCnt}');
                overallFit = polyfit(overallX,overallY,1);
                overallFit = polyval(overallFit,sort(overallX));
                p2=plot(sort(overallX),overallFit,'k-','LineWidth',2);
                h=lsline; % add regerssion line for each mouse sepparately
                for i = 1:numel(h); h(i).LineWidth=0.6;  h(i).LineStyle='-'; h(i).Color(4)=0.3;  end
                legend([p2,h(1)],{'across all mice','individual mice'}); legend boxoff
                [pearson,pVal] = corrcoef(overallX,overallY);
                title([num2str('Pearsons R: '),num2str(pearson(2),'%4.2f'),num2str(' pVal: '),num2str(pVal(2),'%4.2f')])
            
            end
        end
    end
end
%% What drives this clear correlation between time asleep and arousal threshold? 
%It looks like sleep durations <60 s are much more sensitive to stimulation 
% is this true and is this independent of SWA? 
EMGSelect= 4; % EMG threshold
SWASelect= 3; % SWA threshold
selectedIteration = (EMGSelect-1)*numel(nEpochsSWA) + SWASelect;
edges = [0,120,240,360,5000];
[bnd_tSinceW_byMouse,bnd_SWA_byMouse,bnd_sigma_byMouse,bnd_nRDur_byMouse,bnd_stimT_byMouse,bnd_tA_byMouse]= deal(nan(nMice,numel(edges)-1));
[bnd_tSinceW_all,bnd_SWA_all,bnd_sigma_all,bnd_nRDur_all,bnd_stimT_all,bnd_tA_all] = deal(cell(nMice,1));

for mouseCnt = 1:nMice
    tArous_NR = output.tArousal_NR{selectedIteration,mouseCnt};
    preSWA=output.preStimSWA{selectedIteration,mouseCnt};
    preSigma=output.prestimSigma{selectedIteration,mouseCnt};
    nRDur = 4*output.NRBoutDur{selectedIteration,mouseCnt}./60; % convert from epochs to s
    tToLastWake = output.tSinceWake{selectedIteration,mouseCnt}./60;
    stimTime = 4*output.epIdxNRStim{selectedIteration,mouseCnt}./60;

    %normalise to mean
    tArous_NR = 100.*tArous_NR./nanmean(tArous_NR);
    preSWA = 100.*preSWA./nanmean(preSWA);
    preSigma = 100.*preSigma./nanmean(preSigma);
    nRDur = 100.*nRDur./nanmean(nRDur);
    tToLastWake = 100.*tToLastWake./nanmean(tToLastWake);
    stimTime = 100.*stimTime./nanmean(stimTime);
    
    % bin according to arousal delay
    [binning,edges]=discretize(tArous_NR,numel(edges)-1); % sort observations into bins depending on time since last wake epoch
    bnd_tSinceW_all{mouseCnt}= arrayfun(@(x) tToLastWake(binning==x), 1:max(binning),'un',0); % apply binning to arousal delay
    bnd_SWA_all{mouseCnt}= arrayfun(@(x) preSWA(binning==x), 1:max(binning),'un',0);
    bnd_sigma_all{mouseCnt}=arrayfun(@(x) preSigma(binning==x), 1:max(binning),'un',0);
    bnd_nRDur_all{mouseCnt}=arrayfun(@(x) nRDur(binning==x), 1:max(binning),'un',0);
    bnd_stimT_all{mouseCnt}=arrayfun(@(x) stimTime(binning==x), 1:max(binning),'un',0);
    bnd_tA_all{mouseCnt}=arrayfun(@(x) tArous_NR(binning==x), 1:max(binning),'un',0);
    
    % get mean for each bin and mouse
    bnd_tSinceW_byMouse(mouseCnt,:)=cellfun(@nanmean,bnd_tSinceW_all{mouseCnt});  % get mean per bin and mouse
    bnd_SWA_byMouse(mouseCnt,:)=cellfun(@nanmean,bnd_SWA_all{mouseCnt});
    bnd_sigma_byMouse(mouseCnt,:)=cellfun(@nanmean,bnd_sigma_all{mouseCnt});
    bnd_nRDur_byMouse(mouseCnt,:)=cellfun(@nanmean,bnd_nRDur_all{mouseCnt});
    bnd_stimT_byMouse(mouseCnt,:)=cellfun(@nanmean,bnd_stimT_all{mouseCnt});
    bnd_tA_byMouse(mouseCnt,:)=cellfun(@nanmean,bnd_tA_all{mouseCnt});
end
% bnd_tSinceW_byMouse=100.*bnd_tSinceW_byMouse./nanmean(bnd_tSinceW_byMouse,2);
% bnd_SWA_byMouse=100.*bnd_SWA_byMouse./nanmean(bnd_SWA_byMouse,2);
% bnd_sigma_byMouse=100.*bnd_sigma_byMouse./nanmean(bnd_sigma_byMouse,2);
% bnd_nRDur_byMouse=100.*bnd_nRDur_byMouse./nanmean(bnd_nRDur_byMouse,2);
% bnd_stimT_byMouse=100.*bnd_stimT_byMouse./nanmean(bnd_stimT_byMouse,2);
% bnd_tA_byMouse = 100.*bnd_tA_byMouse./nanmean(bnd_tA_byMouse,2);

xLabelString = 'perecentile of arousal delays';
figure
subplot(2,3,1)
boxplot(bnd_tSinceW_byMouse); hold on;
plot(bnd_tSinceW_byMouse','o')
xlim([0 numel(edges)+1]); 
ylabel('time since wake [% mean]'); set(gca,'TickDir','out'); xlabel(xLabelString)
title('time since wake'); box off;

subplot(2,3,2)
boxplot(bnd_SWA_byMouse); hold on;
plot(bnd_SWA_byMouse','o')
xlim([0 numel(edges)+1])
ylabel('SWA [% mean]');  set(gca,'TickDir','out'); xlabel(xLabelString)
title('slow wave activity');box off;

subplot(2,3,3)
boxplot(bnd_sigma_byMouse); hold on;
plot(bnd_sigma_byMouse','o')
xlim([0 numel(edges)+1])
ylabel('sigma power [% mean]');  set(gca,'TickDir','out'); xlabel(xLabelString)
title('sigma power');box off;

subplot(2,3,4)
boxplot(bnd_nRDur_byMouse); hold on;
plot(bnd_nRDur_byMouse','o')
xlim([0 numel(edges)+1])
ylabel('NREM dur [% mean]');  set(gca,'TickDir','out'); xlabel(xLabelString)
title('NREM dur');box off;

subplot(2,3,5)
boxplot(bnd_stimT_byMouse); hold on;
plot(bnd_stimT_byMouse','o')
xlim([0 numel(edges)+1])
ylabel('ZT time  [% mean]');  set(gca,'TickDir','out'); xlabel(xLabelString)
title('ZT time  ');box off;

subplot(2,3,6)
boxplot(bnd_tA_byMouse); hold on;
plot(bnd_tA_byMouse','o')
xlim([0 numel(edges)+1])
ylabel('time to arousal [% mean]');  set(gca,'TickDir','out'); xlabel(xLabelString)
title('time to arousal  ');box off;

% linkaxes();

% nEntries = cellfun(@(x) cellfun(@numel,x),bnd_tSinceW_all,'un',0); % get number of entries per bin for each mouse
% nEntries=cat(1,nEntries{:});
% nEntries = min(nEntries,[],1);

saveas(gcf,fullfile('F:\Tomoko_OptoStim\tArousalVsSWAVsTAsleep',[currCond,'.fig']))
saveas(gcf,fullfile('F:\Tomoko_OptoStim\tArousalVsSWAVsTAsleep',[currCond,'.png']))

%% visualize trial-wise thresholds for one mouse
timeInS = -tPre:0.25:tPost-0.25;
figure;
for stimCnt = 1:size(opto_EMGVar_NR,1)
    plot(timeInS,opto_EMGVar_NR(stimCnt,:)); hold on; % plot variance
    plot(timeInS(threshCrossings_NR(stimCnt,:)),opto_EMGVar_NR(stimCnt,threshCrossings_NR(stimCnt,:)),'ro') % plot variance
    line([-tPre,tPost],[thresholds(stimCnt),thresholds(stimCnt)],'Color','r')
    line([tArous_NR(stimCnt),tArous_NR(stimCnt)],[get(gca,'YLim')],'Color','k','LineWidth',1.4)
    xlabel('time [s]'); ylabel('EMG variance')
    ylim([0 10000])
    pause; hold off;
end
%% visualize mean opto evoked EMG variance 
figure;
iteration = 30; %which threshold to pklot
for mouseCnt = 1:nMice
    subplot(4,2,mouseCnt)
    timeInS = -tPre:0.25:tPost-0.25;
    meanVar= median(optoEvokedEMG_NR{iteration,mouseCnt},1);
    VarPtiles = prctile(optoEvokedEMG_NR{iteration,mouseCnt},[25,75],1);
    meanTArous=mean(tArousal_NR{iteration,mouseCnt}); %mean time to crossing "arousal" threshold
    currThresh = thresholdsUsed{iteration,mouseCnt};
    plot(timeInS,meanVar,'k'); hold on; % plot mean evoked variance
    patch([timeInS,fliplr(timeInS)],[VarPtiles(1,:),fliplr(VarPtiles(2,:))],'k','FaceAlpha',0.3,'EdgeColor','none'); % add 25/75%percentiles as shading
    xlabel('time to stim [s]'); ylabel('log EMG var'); box off; set(gca,'TickDir','out','YScale','log');xlim([0,20]); %cosmetics
    title(['GDCh',num2str(mouseNr.(currCond)(mouseCnt)),'ntrials: ',num2str(size(optoEvokedEMG_NR{iteration,mouseCnt},1))]); % title
    if mouseCnt ==1
        legend('median','25/75 %tiles'); legend boxoff
    end
    line([meanTArous,meanTArous],get(gca,'YLim'),'Color','r')
    line(get(gca,'XLim'),[currThresh,currThresh],'Color','r')

end
linkaxes
%% plot latency to arousal depending on states 
f1=figure('name','ignoring mouse ID');
f2=figure('name','using mean of each mouse');
EMGThresholds= [1,2,3,4,6,8]';
settings.EMGThresholds=EMGThresholds ; % EMG threshold
settings.nEpochsSWA= [1]; % EMG threshold
output=extract_features_from_opto_evoked_arousals(EMGs,stims,scorings,EMGVarRatios,EMGVarBinning,settings); % runs the script
for EMGCnt = 1:numel(EMGThresholds) % cycle through various emg thresholds
    %extract the arousalT based on the current EMG threshold
    tArousNR=output.tArousal_NR(1+((EMGCnt-1)*numel(settings.nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
    tArousR=output.tArousal_R(1+((EMGCnt-1)*numel(settings.nEpochsSWA)),:);
    acrossMiceNR = cell2mat(tArousNR);
    acrossMiceR= cell2mat(tArousR);
    withinMiceNR = cellfun(@mean,tArousNR);
    withinMiceR = cellfun(@mean,tArousR);

    % plot difference NREM REM across all samples, ignoring mouse ID
    figure(f1)
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    %     bar(1,mean(acrossMiceNR),'FaceColor',npgCMap(1,:));bar(2,mean(acrossMiceR),'FaceColor',npgCMap(2,:));
    %     errorbar(1,mean(acrossMiceNR),std(acrossMiceNR)./sqrt(numel(acrossMiceNR)),'Color',npgCMap(1,:));
    %     errorbar(2,mean(acrossMiceR),std(acrossMiceR)./sqrt(numel(acrossMiceR)),'Color',npgCMap(2,:));
    %     xlim([0,3]);
    boxplot([acrossMiceNR,acrossMiceR],[ones(size(acrossMiceNR)),2.*ones(size(acrossMiceR))],'notch','on');
    h = findobj(gca,'Tag','Box');
    h2= findobj(gca,'Tag','Outliers');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),npgCMap(j,:),'FaceAlpha',.4);
        h(j).Color = 'k'; % set box edge color
        h2(j).MarkerEdgeColor=[0.8 0.8 0.8];
    end
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]); 
    ylabel('arousal delay [s]'); ax = gca; ax.XTickLabel ={'NREM','REM'};
    set(gca,'TickDir','out'); box off;set(gcf,'Position',[205 484 1264 415]);
    
    % plot difference NREM REM for each mouse sepparately
    figure(f2)
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    combined = [withinMiceNR;withinMiceR];
    plot(1:2,combined,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3)
    meanVals = mean(combined,2);
    SEMs = std(combined,[],2)./sqrt(size(combined,2));
    display(['mean NR/REM: ', num2str(meanVals'),' SEM: ',num2str(SEMs')])
    errorbar(1:2,meanVals,SEMs,'ok','MarkerFaceColor','k')
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('arousal delay [s]'); ax = gca; ax.XTickLabel ={'NREM','REM'};
    set(gca,'TickDir','out'); box off
    ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415])
    
    p=signrank(withinMiceNR,withinMiceR,'tail','both');
    stats.p=p;
    text(1,max(get(gca,'YLim'))-1,['p= ',num2str(p)])
    linkaxes()
end

%% plot time course of arousal delays
figure;
% get arousal delays for each mouse dependent on what threshold was used (i.e. ignore how many SWA episodes we look back)
tArousNR=output.tArousal_NR(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
tStimNR=output.epIdxNRStim(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal

for EMGCnt = 1:numel(EMGThresholds) % cycle through various emg thresholds
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    currTArous = tArousNR(EMGCnt,:);
    currStimTimes = tStimNR(EMGCnt,:); % get time of stim epoch in  minutes
    for mouseCnt = 1:size(tArousNR,2)
       plot(currStimTimes{mouseCnt}.*4./60,currTArous{mouseCnt},'Color',npgCMap(mouseCnt,:)); hold on; 
       
    end

end
