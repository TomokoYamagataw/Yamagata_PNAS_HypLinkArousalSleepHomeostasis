function [output]=extract_features_from_opto_evoked_arousals(EMGs,stims,scorings,EMGVarRatios,EMGVarBinning,settings,varargin)

% extract settings 
sR = settings.sR;
EMGThresholds = settings.EMGThresholds;
nEpochsSWA = settings.nEpochsSWA;
tPre = settings.tPre;
tPost = settings.tPost;
binSize_EMGVar = settings.binSize_EMGVar; % in seconds
binSize = sR/(1/binSize_EMGVar);%in samples
minDur = settings.minDur;

% innitialise variables
nIter = numel(nEpochsSWA)*numel(EMGThresholds);
nMice = numel(EMGs);
[thresholdsUsed,tArousal_NR,tArousal_R,preStimSWA,arousInt,tSinceWake,optoEvokedEMG_NR,optoEvokedEMG_R,optoEvokedEMG_W,...
    optoEvokedEMG_W_corrected,prestimSigma,prestimTheta,epIdxNRStim,epIdxRStim,preStimSpectra,nrSpectras,NRBoutDuration] = deal(cell(nIter,nMice)); % innitialise containers
iteration = 0;

% extract opto-stim related features of EMG and EEG 
for EMGthresh = EMGThresholds'
    for nEpochsForSWA = nEpochsSWA'
        iteration = iteration + 1
        for mouseCnt = 1:nMice
            
            % load data
            currEMG= EMGs{mouseCnt} ;
            currStim=stims{mouseCnt} ;
            currScoring=scorings{mouseCnt} ;
            EMGVarRatio=EMGVarRatios{mouseCnt};
            EMGVarRatio(isinf(EMGVarRatio))=nan;
            EMGVarBins=EMGVarBinning{mouseCnt};
            
            % process sleepscoring output
            nREpochs = currScoring.nr;
            wEpochs = sort([currScoring.w;currScoring.w1]);
            rEpochs = currScoring.r;
            
            % if we are analysisng a sedation experiment, ignore REM and include sedation scoring
            if isfield(settings,'sedation')
               if  settings.sedation == 1;
                   rEpochs = sort([currScoring.s;currScoring.s4]);
               end
            end
                


            nRSpectra = currScoring.spectr(nREpochs,:);
            nR_SWA = mean(nRSpectra(:,3:17),2); % defining SWA as 0.5 - 4Hz
            nR_sigma = mean(nRSpectra(:,41:61),2); % 10-15 Hz
            nR_theta = mean(nRSpectra(:,21:41),2); % 5-10
            
            % make trial-wise opto-evoked EMG (for stimulations beginning in NREM or REM sleep)
            stimEpochs = floor(currStim(:,1)/(sR*4)); % get epoch idx for each stimulation onset
            [~,stimIdxs_NR,EpIdx_NRStim] = intersect(stimEpochs-1,nREpochs); %  indices of stimulations/NREM epochs containing NREM/stim
            opto_EMG_NR = cell2mat(arrayfun(@(x) currEMG(x-sR*tPre+1:x+sR*tPost),currStim(stimIdxs_NR,1),'un',0)); % trial-wise opto-evoked EMG
            [~,stimIdxs_r,EpIdx_RStim] = intersect(stimEpochs-1,rEpochs); %  indices of stimulations/REM epochs containing REM/stim
            opto_EMG_Rem = cell2mat(arrayfun(@(x) currEMG(x-sR*tPre+1:x+sR*tPost),currStim(stimIdxs_r,1),'un',0)); % trial-wise opto-evoked EMG
            [~,stimIdxs_w,EpIdx_WStim] = intersect(stimEpochs-1,wEpochs); %  indices of stimulations/REM epochs containing REM/stim
            opto_EMG_wake = cell2mat(arrayfun(@(x) currEMG(x-sR*tPre+1:x+sR*tPost),currStim(stimIdxs_w,1),'un',0)); % trial-wise opto-evoked EMG

            nEpochsPre = settings.tPre/4; % how many epochs preceeding stimulation are we looking at?
            priorWEpWstim = arrayfun(@(x) find(and(wEpochs<=wEpochs(EpIdx_WStim(x)),wEpochs>wEpochs(EpIdx_WStim(x))-nEpochsPre)),1:numel(EpIdx_WStim),'un',0);% find closest preceeding wake
            trlsCompleteWaking=cellfun(@numel,priorWEpWstim)==nEpochsPre; % indices of wake stims where animal was continuously awake
            
            % get SWA preceeding stim
            preSWA=arrayfun(@(x) mean(nR_SWA(x-nEpochsForSWA:x-1)),EpIdx_NRStim); % average across x preceedig episode (makes no difference)
%             epsTooShort=arrayfun(@(x) sum(diff(nREpochs(x-nEpochsForSWA:x)))==nEpochsForSWA,EpIdx_NRStim); % flag stims
%             with insufficient nr sleep preceeding stim 

            epsTooShort=arrayfun(@(x) sum(diff(nREpochs(x-nEpochsForSWA:x))),EpIdx_NRStim); % only flag stims with NO nr sleep before stim
            epsTooShort=epsTooShort>0;     
            preSWA(epsTooShort==0)=nan;
            
            % get power in spindle band perceding stim
            preSigma =arrayfun(@(x) mean(nR_sigma(x-nEpochsForSWA:x-1)),EpIdx_NRStim); % average across x preceedig episode (makes no difference)
            preSigma(epsTooShort==0)=nan;
            
            % get power in theta band preceeding stim
            preTheta =arrayfun(@(x) mean(nR_theta(x-nEpochsForSWA:x-1)),EpIdx_NRStim); % average across x preceedig episode (makes no difference)
            preTheta(epsTooShort==0)=nan;
            
            preSpectra = cell2mat(arrayfun(@(x) mean(nRSpectra(x-nEpochsForSWA:x-1,:),1),EpIdx_NRStim,'un',0));
            preSpectra(epsTooShort==0,:)=nan;
            % get time to last waking episode preceeding stim
            priorWEp = arrayfun(@(x) find(wEpochs<nREpochs(EpIdx_NRStim(x)),1,'last'),1:numel(EpIdx_NRStim),'un',0);% find closest preceeding wake
            priorWEp(cellfun(@isempty,priorWEp))={nan}; priorWEp = cell2mat(priorWEp);
            nNans = sum(isnan(priorWEp));
            if isempty(priorWEp)
                tToLastWake = nan;
            elseif nNans>0 % in case the first waking period is not in the recording
                tToLastWake = 4*(stimEpochs(stimIdxs_NR(nNans+1:end))-wEpochs(priorWEp(nNans+1:end))); % ignore first stim and calculate difference in epochs * epochlength
                tToLastWake = [nan(nNans,1);tToLastWake]; % add nan value for first stim
            else
                tToLastWake = 4*(stimEpochs(stimIdxs_NR)-wEpochs(priorWEp)); % difference in epochs * epochlength
            end
            
            % get time since beginning of the NR bout 
            NREMEpochs =  sort([-1;currScoring.nr;currScoring.nr2;currScoring.mt]); % combine all the relevant epochs and add -1 in the beginning to make sure the first bout gets measured properly in the next line
            [epochs,~,EpIdx_NREMStim] = intersect(stimEpochs-1,NREMEpochs); %  indices of stimulations/NREM epochs containing NREM/stim
            %possible that some stims have now been detected that happen during brief awakenings, get rid of that below
            if length(EpIdx_NREMStim)~=length(EpIdx_NRStim)
               nonMTStims = arrayfun(@(x) find((nREpochs(EpIdx_NRStim(x))-epochs)==0),1:numel(EpIdx_NRStim)); % detect stims during mt
               EpIdx_NREMStim = EpIdx_NREMStim(nonMTStims); % drop the wrongly detected stims
            end
            NRstartIdx=arrayfun(@(x) find(diff(NREMEpochs(1:x))>1,1,'last'),EpIdx_NREMStim); % only flag stims with NO nr sleep before stim
            NRBoutDur = EpIdx_NREMStim - NRstartIdx; % calcualte duration of NR bout
                        
            % get trial-wise EMG variance in 250 (+/- <1) ms bins
            opto_EMGVar_NR = squeeze(var(reshape(opto_EMG_NR,size(opto_EMG_NR,1),binSize,size(opto_EMG_NR,2)/binSize),[],2));
            opto_EMGVar_R = squeeze(var(reshape(opto_EMG_Rem,size(opto_EMG_Rem,1),binSize,size(opto_EMG_Rem,2)/binSize),[],2));
            opto_EMGVar_W = squeeze(var(reshape(opto_EMG_wake,size(opto_EMG_wake,1),binSize,size(opto_EMG_wake,2)/binSize),[],2));
            opto_EMGVar_W_corrected = opto_EMGVar_W(trlsCompleteWaking,:); % only takle trials with complete waking in tPre
            
                        % set "arousal threshold based on thresholding the variance"
%                         preStim = opto_EMGVar_NR(:,1:tPre/binSize_EMGVar);
%                         postStim = opto_EMGVar_NR(:,tPre/binSize_EMGVar+1:end);
%                         thresholds = mean(preStim,2)+EMGthresh*std(preStim,[],2); % define a threshold @ mean + 5* sd
            
            % define threshold based on statistics of EMG variance during waking compared to sleep
            % normally, the threshold is calculated for each condition sepparately, but it can make sense to use the same
            % threhsold for several conditions, in that case settings.forceThreshold can be set 
            if isfield(settings,'forceThreshold')
                if iscell(settings.forceThreshold)  % monkeypatch, this field can come as a cell or matrix, gotta fix this someday
                ratioBasedThresh = settings.forceThreshold{find(EMGThresholds==EMGthresh),mouseCnt};
                else
                     ratioBasedThresh = settings.forceThreshold(find(EMGThresholds==EMGthresh),mouseCnt);
                end    
            else
                ratioBasedThresh = EMGVarBins(find(EMGVarRatio>EMGthresh,1));
                ratioBasedThresh = 10.^ratioBasedThresh; % convert back from log
            end
            thresholds = ratioBasedThresh; % use the emg threshold based on
            
            threshCrossings_NR = opto_EMGVar_NR>thresholds; % apply threshold to NREM
            threshCrossings_NR(:,1:tPre*(1/binSize_EMGVar)) = 0; % avoid pre-stim threshold crossings
            thresholds = mean(thresholds);
            threshCrossings_R = opto_EMGVar_R>thresholds; % apply threshold to REM
            threshCrossings_R(:,1:tPre*(1/binSize_EMGVar)) = 0; % avoid pre-stim threshold crossings  
            
            % reject arousal durs <x seconds (by smoothing)
            smoothKernel = [zeros(1,minDur), ones(1,minDur)/minDur 0]; % make
            longThreshCrossings_NR=conv2(threshCrossings_NR,fliplr(smoothKernel),'same')==1; % reject arousal durs <2 seconds (by smoothing)
            longThreshCrossings_R=conv2(threshCrossings_R,fliplr(smoothKernel),'same')==1; % reject arousal durs <2 seconds (by smoothing)

            
            % calculate time of first arousal relative to stim onset
            tArous_NR= arrayfun(@(x) find(longThreshCrossings_NR(x,:),1),1:size(opto_EMGVar_NR,1),'un',0); % get first threshold crossing
            tArous_NR(cellfun(@ isempty,tArous_NR))={nan};  % replace empty entries by nan
            tArous_NR = cell2mat(tArous_NR)*0.25 - tPre; % convert to matrix and from "bins" back to seconds
            tArous_R= arrayfun(@(x) find(longThreshCrossings_R(x,:),1),1:size(opto_EMGVar_R,1),'un',0); % get first threshold crossing
            tArous_R(cellfun(@ isempty,tArous_R))={nan};  % replace empty entries by nan
            tArous_R = cell2mat(tArous_R)*0.25 - tPre; % convert to matrix and from "bins" back to seconds
            tmpTArous = tArous_NR; tmpTArous(isnan(tmpTArous))=0; %for the below line only
            arousInt{iteration,mouseCnt} = sum(threshCrossings_NR,2)./((tPost-tmpTArous-0.25)*4)'; % calculate arousal intensity:proabability (vaguely) of being above threshold after threshold was crossed
            
            % store arousal delays and SWA from this mouse
            tArousal_NR{iteration,mouseCnt} = tArous_NR;
            tArousal_R{iteration,mouseCnt} = tArous_R;
            epIdxNRStim{iteration,mouseCnt}=nREpochs(EpIdx_NRStim); % store the timing of the stimulation
            epIdxRStim{iteration,mouseCnt}=EpIdx_RStim; % store the timing of the stimulation
            thresholdsUsed{iteration,mouseCnt}= thresholds; 
            
            preStimSWA{iteration,mouseCnt} = preSWA;
            prestimSigma{iteration,mouseCnt} = preSigma;
            prestimTheta{iteration,mouseCnt} = preTheta;
            preStimSpectra{iteration,mouseCnt} = preSpectra;
            tSinceWake{iteration,mouseCnt} = tToLastWake;
            optoEvokedEMG_NR{iteration,mouseCnt}=opto_EMGVar_NR;
            optoEvokedEMG_R{iteration,mouseCnt}=opto_EMGVar_R;
            optoEvokedEMG_W{iteration,mouseCnt}=opto_EMGVar_W;
            optoEvokedEMG_W_corrected{iteration,mouseCnt}=opto_EMGVar_W_corrected;
            nrSpectras{iteration,mouseCnt}=mean(nRSpectra,1);
            NRBoutDuration{iteration,mouseCnt} = NRBoutDur; 
        end
    end
end

output.tArousal_NR = tArousal_NR;
output.tArousal_R = tArousal_R;
output.epIdxNRStim = epIdxNRStim;
output.epIdxRStim = epIdxRStim;
output.thresholdsUsed = thresholdsUsed;
output.preStimSWA = preStimSWA;
output.prestimSigma = prestimSigma;
output.prestimTheta = prestimTheta;
output.prestimSpectra = preStimSpectra;
output.tSinceWake = tSinceWake;
output.optoEvokedEMG_NR = optoEvokedEMG_NR;
output.optoEvokedEMG_R = optoEvokedEMG_R;
output.optoEvokedEMG_W=optoEvokedEMG_W;
output.optoEvokedEMG_W_corrected=optoEvokedEMG_W_corrected;
output.nrSpectra = nrSpectras;
output.NRBoutDur= NRBoutDuration; 

% if we are analysing a sedation experiment
if isfield(settings,'sedation')
    if  settings.sedation == 1
        output.tArousal_R = [];
        output.optoEvokedEMG_R = [];
        output.tArousal_Sed = tArousal_R;
        output.optoEvokedEMG_Sed = optoEvokedEMG_R;  
    end
end


