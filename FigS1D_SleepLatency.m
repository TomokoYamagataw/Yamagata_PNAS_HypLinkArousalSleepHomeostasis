clear all
close all

path='I:\optogenetics\'; 

% mousenames1=[1 2 3 4 5 6 7 8]; % names of GFP ctrl mice indicated in the file name
% days1=['250218 260218';'290418 120418';'300618 260618';'300618 040718';'031018 061018';'081018 061018';'081018 061018';'080119 070119'];

% mousenames2=[5 6 8 15 16 17];
% days2=['020418 030418';'290418 120418';'100618 050618';'181218 191218';'181218 191218';'181218 191218']; 

% LPO
mousenames2=[5 15 16 17 18 19 20 21];
days2=['030418';'191218';'191218';'191218';'160819';'160819';'180919';'160819'];

% nonLPO
mousenames4=[1 6 7 8 9 12]; %non-LPO %GDCh8 - 140518
days4=['140218';'120418';'130518';'150618';'150618';'090618'];


geno = 2;

ders=strvcat('fro','occ');
der=1; 
deri=ders(der,:); 

maxep=21600; 
epochs=1:maxep; 
epochs=epochs./900; 

epochl=4; 

zermat=zeros(1,maxep); %creates an array of zeros with 1 row and as long as the total number of epoch

fs=256;  %sampling rate (you sample 256 times in each 1second)


%we want to analyze period from stimulus onset to the next stimulus (20min
%in between stimuli)
after=18;


ne=15; %there are 15 epochs in 1min

x=0:1:after*ne-1; %creates an array of epoch indices in the 2min before to after period
x=x*epochl; %converts the epoch indices to seconds

pathvs1=[path,'outputVSgfp\']; %defines the folder with output files of vigilance states for mice1 (GFP ctrl)
pathvs2=[path,'outputVSchr\']; %defines the folder with output files of vigilance states for mice2 (ChR2 exp)
pathstim=[path,'STIMs\']; %gives the folder with stimulation inf
%in the files there are two variables: sigstim and startend (giving the
%beggining and end timepoint of each stimulation)

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

if geno==1 %geno=1 ctrl group
    mousenames=mousenames1;days=days1; pathvs=pathvs1; gn='GFP'; %firststims=firststims1;
elseif geno==2
    mousenames=mousenames2; days=days2; pathvs=pathvs2; gn='GDCh';  %firststims=firststims2;
elseif geno==4
    mousenames=mousenames4; days=days4; pathvs=pathvs2; gn='GDCh';  %firststims=firststims2;
end
numanim=length(mousenames); %how many animals were recorded

WallMice=[]; %wake state of all mice
%empty array that will be filled with wake data from each
%recorded mouse
NallMice=[];
RallMice=[];

steps=[0:30:600];
%this loop goes though all the animals
%Question: what is the latency to go to sleep after stim end

%empty array to be filled with latency values of all mice
LatNallMice=[];
LatRallMice=[];
LatWallMice=[];

H1=[];
H2=[];

CorrN=NaN(72,11);CorrR=NaN(72,11);CorrW=NaN(72,11);
PrevSleep=NaN(72,11);

figure(1)
hold on 

SleepLatencyW=[]; SleepLatencyN=[]; SleepLatencyR=[]; 
SleepLatencyAll=NaN(72,numanim*3); SleepLatencyAllv2=NaN(72,numanim*3); 
for anim=1:numanim %go through this loop as many times as there are animals
    
    mousename=[gn,num2str(mousenames(anim))];
    %gives the name for each mouse as it would be written in the output file name
    %Genotype + no. of the mouse (e.g. GDCh1)
    day=days(anim,:); day(isspace(day))=[];
    
    fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
    eval(['load ',pathstim,fn0,'.mat startend ']);
    % --> loads the startend variable
    
    %converts startend from sample index to seconds
    stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
    %fs*epoch1 = how many samples in total in one epoch (256samples/sec * 4s = 1024)
    %converts sample no. to epoch index
    stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2)
    
    fn1=[mousename,'-',day,'-',deri]; %makes the full file name for the vigilance state and derivation output desired
    
    eval(['load ',pathvs,fn1,'.mat w nr r w1 nr2 r3 mt ma bastend ']);
    %loads the vigilance states file and its stated variable
    % w/nr/r and others give an indices of epochs with that state
    
    VS=zeros(7,maxep); %creates a zeros array with as many fields as max. no. of epochs
    %the rows will present each of the vigilance states
    VS(1,w)=1; %row1=WAKE and at each of the epochs indicated by w put 1 to signify +ve state
    VS(2,w1)=1; %row2=WAKE ARTIFACT
    VS(3,nr)=1; %row3=NREM
    VS(4,nr2)=1; %row4=NREM ARTIFACT
    VS(5,r)=1; %row5=REM
    VS(6,r3)=1; %row6=REM ARTIFACT
    VS(7,mt)=1; %row7=MOMENT ARTIFACT (BRIEF MOVEMENT)
    clear w nr r w1 nr2 r3 mt ma bastend; %clears the variables for the next loop
    
    wake=sum(VS([1:2 7],:)); %extract rows 1 to 2 and all columns;
    %sum to get single row of epochs spent in wake including artifacts
    %here we only care about the state not the power
    %spectra so can include the artifacts as well
    nrem=sum(VS(3:4,:));
    rems=sum(VS(5:6,:));
    
    sleep=nrem+rems; %gives total sleep state epochs
    
    stimep=stimep1; %start of stimulation episodes
    stimend=stimep2; %end of stimulation episodes
    
    %exclude values outside of the possible range 
    out0=find(stimep+after*ne>maxep); stimep(out0)=[];

    nums=length(stimep); %how many stimulation episodes there are
   
    Ws=[]; Ns=[]; Rs=[]; Ts=[]; %creates empty arrays to be filled with vigilance states and for the desired observation period
    Ltotal=[]; Ln1=[]; Ln2=[]; Lw1=[]; Lw2=[]; Lr2=[];     

    for s=1:nums %iterating over all the stimulations
        
        step=stimep(s);
        stepend=stimend(s);
        eps=stepend+1:stepend+after*ne; %returns the desired period
        stateEps20sec=step-5:step-1; %for determining the initial state
        if step-15>0
        stateEps1min=step-15:step-1; 
        end
        
        ws=wake(eps);  Ws=[Ws;ws];
        ns=nrem(eps);  Ns=[Ns;ns];
        rs=rems(eps);  Rs=[Rs;rs];
        ts=sleep(eps); Ts=[Ts;ts];
        %prevSleep=[prevSleep; NaN NaN];
        
        if sum(ts)==0; mints=length(ts)+1; else; mints=min(find(ts>0));  end
        mints=mints*4/60; % min
        Ltotal=[Ltotal; mints];

        if nanmean(nrem(stateEps1min))==1; Ln1=[Ln1; s mints]; end
        if nanmean(wake(stateEps1min))==1; Lw1=[Lw1; s mints]; end
        if nanmean(nrem(stateEps20sec))==1; Ln2=[Ln2; s mints]; end
        if nanmean(wake(stateEps20sec))==1; Lw2=[Lw2; s mints]; end
        if nanmean(rems(stateEps20sec))==1; Lr2=[Lr2; s mints]; end

    end

    SleepLatencyW=[SleepLatencyW;Lw1];
    SleepLatencyN=[SleepLatencyN;Ln1];
    SleepLatencyR=[SleepLatencyR;Lr2];
    
    SleepLatencyAll(1:size(Lw1,1),anim)=Lw1(:,2); 
    SleepLatencyAll(1:size(Ln1,1),numanim*1+anim)=Ln1(:,2); 
    SleepLatencyAll(1:size(Lr2,1),numanim*2+anim)=Lr2(:,2);
    SleepLatencyAllv2(1:size(Lw1,1),anim)=Lw1(:,2); 
    SleepLatencyAllv2(1:size(Ln1,1),numanim*1+anim)=Ln1(:,2); 
    SleepLatencyAllv2(1:size(Lr2,1),numanim*2+anim)=Lr2(:,2);
end

Prism_SleepLatencyW_withnoSleep=SleepLatencyW(:,2);
Prism_SleepLatencyN_withnoSleep=SleepLatencyN(:,2);
Prism_SleepLatencyR_withnoSleep=SleepLatencyR(:,2);

%Light
LW=SleepLatencyW;
LN=SleepLatencyN;
LR=SleepLatencyR;

out1=find(LW(:,1)>36); LW(out1,:)=[];
out2=find(LN(:,1)>36); LN(out2,:)=[];
out3=find(LR(:,1)>36); LR(out3,:)=[];

aPrism_all_SleepLatency_L=NaN(max([length(LW),length(LN),length(LR)]),3);
aPrism_all_SleepLatency_L(1:length(LW),1)=LW(:,2);
aPrism_all_SleepLatency_L(1:length(LN),2)=LN(:,2);
aPrism_all_SleepLatency_L(1:length(LR),3)=LR(:,2);
samplesize_L=size(aPrism_all_SleepLatency_L,1);
avrLatency_L=nanmean(aPrism_all_SleepLatency_L)
seLatency_L=nanstd(aPrism_all_SleepLatency_L)/sqrt(samplesize_L)


%Dark
LW=SleepLatencyW;
LN=SleepLatencyN;
LR=SleepLatencyR;

out1=find(LW(:,1)<=36); LW(out1,:)=[];
out2=find(LN(:,1)<=36); LN(out2,:)=[];
out3=find(LR(:,1)<=36); LR(out3,:)=[];

aPrism_all_SleepLatency_D=NaN(max([length(LW),length(LN),length(LR)]),3);
aPrism_all_SleepLatency_D(1:length(LW),1)=LW(:,2);
aPrism_all_SleepLatency_D(1:length(LN),2)=LN(:,2);
aPrism_all_SleepLatency_D(1:length(LR),3)=LR(:,2);
samplesize_D=size(aPrism_all_SleepLatency_D,1);
avrLatency_D=nanmean(aPrism_all_SleepLatency_D)
seLatency_D=nanstd(aPrism_all_SleepLatency_D)/sqrt(samplesize_D)





% NREM
outN=find(SleepLatencyN==18); SleepLatencyN(outN,:)=[];
numdidnotreturntosleep_NREM=size(outN,1)
samplesize_N=size(SleepLatencyN,1)
avrLatencyNREM=mean(SleepLatencyN(:,2))
seLatencyNREM=std(SleepLatencyN(:,2))/sqrt(samplesize_N)

% REM
outR=find(SleepLatencyR==18); SleepLatencyR(outR,:)=[];
numdidnotreturntosleep_REM=size(outR,1)
samplesize_R=size(SleepLatencyR,1)
avrLatencyREM=mean(SleepLatencyR(:,2))
seLatencyREM=std(SleepLatencyR(:,2))/sqrt(samplesize_R)

% Wake
outW=find(SleepLatencyW==18); SleepLatencyW(outW,:)=[];
numdidnotreturntosleep_Wake=size(outW,1)
samplesize_W=size(SleepLatencyW,1)
avrLatencyWake=mean(SleepLatencyW(:,2))
seLatencyWake=std(SleepLatencyW(:,2))/sqrt(samplesize_W);

aPrism_all_SleepLatency_L;
aPrism_all_SleepLatency_D;
aPrism_all_SleepLatency_D