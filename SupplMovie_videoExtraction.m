%% note: this is just to compare PIR with EMG variance, everything else use V2
clear all

%colormap
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 240 228 66;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue,yellow
npgCMap = [npgCMap;npgCMap];
tomokoColors =[0 0 255;0 192 0; 249 64 64]./255;

% USER INPUT
pathToVideo = 'I:\optogenetics\Movie' ; % where the videos are stored
pathToData = 'F:\Tomoko_OptoStim\EMGAnalysis_matfiles'; % where the mat file called "TomokoVideoEMGs" is saved 
currCond_EMG ='TenHz24h'; % name of matfile
fileNames{1}='20190816_GDCh18_24h_10Hz.mp4';
fileNames{2}='20190816_GDCh19_24h_10Hz.mp4';
fileNames{3} = '20190816_GDCh21_24h_10Hz.mp4';

% set start time of video
% GDCh18:
outside_t_start_video(1) = 8*3600 + 51*60 + 11; % outside time is 08:51:11
% GDCh19:
outside_t_start_video(2) = 8*3600 + 51*60 + 18; % outside time is 08:51:18
% GDCh21:
outside_t_start_video(3) = 8*3600 + 51*60 + 33; % outside time is 08:51:33
mouseNumberingForEMG = [5,6,8];

%% load EMG and sleep scoring
test2=load(fullfile(pathToData,currCond_EMG));

test2= struct(); 
test2.EMGs = test.EMGs; test2.EMGs(1:4)=cell(4,1);
test2.scorings  = test.scorings;
test2.settings = test.settings;
save(fullfile(pathToData,'TomokoVideoEMGs'),'test2');

mouseNames = {'GDCh18','GDCh19','GDCh21'};
%% extract EMG variance and spectra in spontaneous wake copared to wake with stim
sR_EMG = test2.settings.sR;
sRVideo = 14.9990;
nMice = 3;
[vidVar_stim_L,vidVar_Nostim_L,w_EMGVar_Nostim_L,w_EMGVar_stim_L,...
    vid_noStim_L,vid_noStim,vid_Stim_L,vidPos_stim_L,vidPos_noStim_L]= deal(cell(nMice,1));
spectraStim = struct();
outside_t_start_Ephys = 9*3600  + 14;% tdt recording start T is 09:00:14 for all mice

tPreStim = 30; % in seconds
tPostStim = 30; % in seconds


[allVids,allEMGs] = deal(cell(3,1));
for mouseCnt = 1:2
    % extract the relevant data from the loaded mat file
    currEMG= test2.EMGs{mouseNumberingForEMG(mouseCnt)} ;
    currCond_EMG=test2.stims{mouseNumberingForEMG(mouseCnt)} ;
    currScoring=scorings{mouseNumberingForEMG(mouseCnt)} ;
    
    % use sleepscoring wake epochs and epochs with stimulation
    wEpochs = currScoring.w;
    nrEpochs = currScoring.nr;
    stimEpochs = floor(currCond_EMG/(sR_EMG*4)); % get epoch idx for each stimulation onset
    stimEpochs = cell2mat(arrayfun(@(x)stimEpochs(x,1):stimEpochs(x,2)',1:size(stimEpochs,1),'un',0));
    stimEpochs = stimEpochs(stimEpochs<(12*3600)./4); % only use epochs during light phase
    stimOnsetEpochs = stimEpochs([1,find(diff(stimEpochs)~=1)+1]); % get epoch onset
    
    % Tomoko: the below line determines which epochs will be selected, specifically, what vigilance state the preceding epoch has to be 
    % it is currently set to only extract videos where the mouse was in NREM before stimulation, to change this, change "nrEpochs"
    % to "wEpochs"
    stimOnsetEpochs_w = arrayfun(@(x) sum(nrEpochs==(x-1))>0,stimOnsetEpochs);
    stimOnsetEpochs(~stimOnsetEpochs_w)=[];
    
    % LOAD VIDEO
    vidobj = VideoReader(fullfile(pathToVideo,fileNames{mouseCnt}));
    nFrames = floor(vidobj.Duration.*vidobj.FrameRate);
    videoTime = 2/sRVideo:1/sRVideo:12*3600; %
    startT_Ephys = outside_t_start_Ephys - outside_t_start_video(mouseCnt);
    [allVids{mouseCnt},allEMGs{mouseCnt}]= deal(cell(numel(stimOnsetEpochs),1));
    for stimCnt = 1:numel(stimOnsetEpochs)
        currStim =  stimOnsetEpochs(stimCnt)*4 - 4;
        currStart = currStim - tPreStim + startT_Ephys;
        currStop = currStim + tPostStim + startT_Ephys;
        vidobj.CurrentTime = currStart;
        
        while hasFrame(vidobj)
            % checkpoint to stop collecting frames
            currEndSample = vidobj.CurrentTime;
            if currEndSample > currStop
                break
            end
            
            image = readFrame(vidobj);
            image = squeeze(image(:,:,1)); % remove non-existent RGB components
            
            % add marker for stim onset
            if and(currEndSample>(currStart+tPreStim),currEndSample<(currStart+tPreStim+120))
                image(440:480,600:640)=ones(41,41);
            else
                image(440:480,600:640)=ones(41,41).*255;
            end
            
            if currEndSample ==currStart
                currVid = image;
            else
                currVid = cat(3,currVid,image);
            end
        end
        allVids{mouseCnt}{stimCnt}= currVid;
        extractedEMG = currEMG(sR_EMG.*(stimOnsetEpochs(stimCnt)*4 - 5 - tPreStim):sR_EMG.*(stimOnsetEpochs(stimCnt)*4 - 5 + tPostStim));
        allEMGs{mouseCnt}{stimCnt}= extractedEMG;
    end
    
end
%% save video if you want 
save(fullfile(pathToVideo,'wakeStimVideos'),'allVids','mouseNames','-v7.3');
%% make a video of EMG + camera 

for mouseCnt = 2%1:nMice
    for vidCnt = 1:1;%numel(allVids{mouseCnt})
        currVid = allVids{mouseCnt}{vidCnt};
        currEMG = allEMGs{mouseCnt}{vidCnt};
        nFrames = size(currVid,3);
        
        v = VideoWriter(fullfile(pathToVideo,[mouseNames{mouseCnt},'_vid_',num2str(vidCnt),'.avi']));
        open(v);
        
        figure;
        ax1 = subplot(2,1,1);
        displayedEMG = nan(size(currEMG));
        displayedEMG(end)= currEMG(end);
        ax2=subplot(2,1,2); hold on;
        timeInS = 0:1/sR_EMG:floor(nFrames./sRVideo);
        plot(timeInS,displayedEMG); xlabel('time [s]'); ylabel('EMG amplitude [uV]');
        set(gcf,'Position',[1142         254         417         685]);

        previousSample = 0;
        for frameCnt = 1:nFrames
            imshow(currVid(:,:,frameCnt),'parent',ax1)
            currEndSample = floor(sR_EMG.*frameCnt./sRVideo);
            if currEndSample>length(currEMG)
                currEndSample=length(currEMG);
            end 
            plot(timeInS(previousSample+1:currEndSample),currEMG(previousSample+1:currEndSample),'k','parent',ax2);
            drawnow()
            currFrame = getframe(gcf);
            writeVideo(v,currFrame);
        end
        close(v);
        close all
    end
end