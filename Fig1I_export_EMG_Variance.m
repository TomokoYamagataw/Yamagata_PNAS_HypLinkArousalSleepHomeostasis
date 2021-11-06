% path to files etc
% paths.rawSignals = 'F:\Tomoko_OptoStim\OutputSignals';
% folderToSave = 'F:\Tomoko_OptoStim\OutputSignals';
paths.rawSignals = 'I:\OptoMod\OutputSignals';
folderToSave = 'I:\OptoMod\outputEMGs';

% set parameters
settings.sR = 256; %sampling rate
settings.anmPrefix = 'GDCh'; % animal name prefix
settings.binSize_EMGVar = 0.25; % over how many seconds will emg variance be calculated

% animal IDs 
mouseNr.TenHz24h = [1,5,6,7,8,9,12,15:21];
mouseNr.TenHz24hWithExclusion = [5,15:19,21];
mouseNr.TenHz24hBaseline = [5,15:19,21];
mouseNr.fiveHz24h = 15:21;
mouseNr.fiveHz24hWithExclusion = [15:19,21];
mouseNr.twoHz24h = 15:21;
mouseNr.oneHz24h = 15:21;
mouseNr.HSP = [5,15:19,21]; % high sleep pressure
mouseNr.LSP = [5,15:19,21]; % low sleep pressure
mouseNr.Dex = [1,5:9,12,16:21]; % Dex

% dates 
dates.TenHz24h = [140218 030418 120418 130518 050618 150618 090618 191218 191218 191218 160819 160819 180919 160819];  %160819 
dates.twoHz24h = [ 251218 251218 251218 190819 190819 190819 190819];
dates.fiveHz24h = [261218 261218 261218	200819 200819 200819 200819];
dates.fiveHz24hWithExclusion = [261218 261218 261218 200819 200819 200819];
dates.TenHz24hWithExclusion = [030418 191218 191218 191218 160819 160819 160819];
dates.TenHz24hBaseline = [020418 181218 181218 181218 100819 100819 100819];
dates.oneHz24h = [241218 241218 241218 190819 190819 190819 190819];
dates.HSP = [310318 211218 211218 211218 060919 060919 060919 060919];
dates.LSP = [300318 201218 201218 201218 040919 040919 040919 040919];
dates.Dex = [010518 010518 050518 150518 220618 220618 220618 150119 150119 220819 220819 220819 220819];

% choose which condition to use
currCond = 'Dex';  %TenHz24h %TenHz24hBaseline


%% calculate EMG variance and save
nMice = numel(mouseNr.(currCond));
for mouseCnt = 1:nMice %13:13 %1:nMice
    currMouse = [settings.anmPrefix,num2str(mouseNr.(currCond)(mouseCnt))];
    currDate = dates.(currCond)(mouseCnt);
    if isnan(currDate); error('trying to access nonexistent recording'); end
    fileNames = structfun(@(x)dir(fullfile(x,[currMouse,'*',num2str(currDate,'%06d'),'*'])),paths,'un',0); % search all relevant folders for mousename & date

    % load emg signal (sR: 256)
    currEMGIdx = arrayfun(@ (x)contains(fileNames.rawSignals(x).name,'emg'),1:numel(fileNames.rawSignals)); % find emg (as opposed to EEG etc)
    currEMGIdx = fileNames.rawSignals(find(currEMGIdx,1)).name;
    currEMG = load(fullfile(paths.rawSignals,currEMGIdx)); % load emg
    currEMG = currEMG.resampled_sig;
    
    % filter EMG using Vlad's filter
    sR =settings.sR;
    p1=5; p2=100; s1=4; s2=120;
    Wp=[p1 p2]/(sR/2); Ws=[s1 s2]/(sR/2); Rp=3; Rs=30; [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
    [bb1,aa1]=cheby2(n,Rs,Wn);
    currEMG=filtfilt(bb1,aa1,currEMG);
    
    % make EMG variance 
    binSize = sR/(1/settings.binSize_EMGVar);
    % it can happen that the last bin is a few samples too short, fill these samples with nan
    if rem(numel(currEMG),binSize)~=0
        currEMG(numel(currEMG):ceil(numel(currEMG)./binSize)*binSize)=nan;
    end
    EMGVar = squeeze(nanvar(reshape(currEMG,size(currEMG,1),binSize,size(currEMG,2)/binSize),[],2));
    
    %%%%%%%% added by Tomoko
    currEMGsec=currEMG;
    if rem(numel(currEMGsec),sR)~=0
        currEMGsec(numel(currEMGsec):ceil(numel(currEMGsec)./sR)*sR)=nan;
    end
    EMGVarSec = (nanvar(reshape(currEMGsec,sR,[])))';    %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % save
    fileName = [currMouse,'-',num2str(currDate),'-EMGVar'];
    save(fullfile(folderToSave,fileName),'EMGVar','EMGVarSec')
end
