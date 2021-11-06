clear all
close all

path='G:\OptoHDbackup\optogenetics\'; 

mousenames1=[1 2 3 4 5 6 7 8]; % names of GFP ctrl mice indicated in the file name 
days1=['120518';'050518';'180718';'180718';'111018';'221018';'221018';'150119'];

% LPO
mousenames2=[5 16 17 18 19 20 21]; 
days2=['010518';'150119';'150119';'220819';'220819';'220819';'220819']; 

% nonLPO
mousenames4=[1 6 7 8 9 12]; 
days4=['010518';'050518';'150518';'220618';'220618';'220618']; 

geno = 2;

tfstim=9; %1 %9
tnstim=12; %4 %12  %%%% stim per 225, jitter 22.5, stimlength=30, - 4times: 675+jitter=730, 5times: 900+jitter=955; 6times: 1125+jitter=1180, 8times: 1800+jitter=1630, 16times: 1800+jitter=3430

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
before=2;
after=0;


ne=15; %there are 15 epochs in 1min

x=-before*ne:1:after*ne; %creates an array of epoch indices in the 2min before to after period
x=x*epochl; %converts the epoch indices to seconds

pathvs1=[path,'outputVSgfp\']; %defines the folder with output files of vigilance states for mice1 (GFP ctrl)
pathvs2=[path,'outputVSchr\']; %defines the folder with output files of vigilance states for mice2 (ChR2 exp)
pathstim=[path,'STIMs\']; %gives the folder with stimulation inf
pathout=[path,'Dex\']; %mkdir(pathout)
pathfig=[path,'Figures_Dex\']; %mkdir(pathout)
%in the files there are two variables: sigstim and startend (giving the
%beggining and end timepoint of each stimulation)

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

if geno==1 %geno=1 ctrl group
    mousenames=mousenames1;days=days1; pathvs=pathvs1; gn='GFP'; %firststims=firststims1;
elseif geno==2
    mousenames=mousenames2; days=days2; pathvs=pathvs2; gn='GDCh';  %firststims=firststims2;
elseif geno==4
    mousenames=mousenames4; days=days2; pathvs=pathvs2; gn='GDCh';  %firststims=firststims2;
end
numanim=length(mousenames); %how many animals were recorded

WallMice=[];NallMice=[];RallMice=[];
TallMice=[];MallMice=[];DallMice=[];

Ls=NaN(numanim*tnstim, 4);Lmeans=[];
 
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

WakeLatencyW=[]; WakeLatencyN=[]; WakeLatencyR=[]; 
WakeLatencyAll=NaN(72,numanim*3); WakeLatencyAllv2=NaN(72,numanim*3); 
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
    
    eval(['load ',pathvs,fn1,'.mat w nr r s w1 nr2 r3 s4 mt ma bastend s ']);
    %loads the vigilance states file and its stated variable
    % w/nr/r and others give an indices of epochs with that state
    
    VS=zeros(9,maxep); %creates a zeros array with as many fields as max. no. of epochs
    %the rows will present each of the vigilance states
    VS(1,w)=1; %row1=WAKE and at each of the epochs indicated by w put 1 to signify +ve state
    VS(2,w1)=1; %row2=WAKE ARTIFACT
    VS(3,nr)=1; %row3=NREM
    VS(4,nr2)=1; %row4=NREM ARTIFACT
    VS(5,r)=1; %row5=REM
    VS(6,r3)=1; %row6=REM ARTIFACT
    VS(7,mt)=1; %row7=MOMENT ARTIFACT (BRIEF MOVEMENT)
    VS(8,s)=1;
    VS(9,s4)=1;
    clear w nr r s w1 nr2 r3 s4 mt ma bastend; %clears the variables for the next loop
    
    wake=sum(VS([1:2 7],:));
    dexs=sum(VS(8:9,:));
    
    stimep=stimep1(tfstim:tfstim+tnstim-1);; %start of stimulation episodes
    stimend=stimep2(tfstim:tfstim+tnstim-1);; %end of stimulation episodes
    
    %exclude values outside of the possible range
    %w/o this --> error: Array indices must be positive integers or logical values.
    out1=find(stimep-before*ne<1);  stimep(out1)=[];  
    out2=find(stimep+after*ne>maxep); stimep(out2)=[];

    nums=length(stimep); %how many stimulation episodes there are
   
    Ws=[]; Ns=[]; Rs=[]; Ds=[]; %creates empty arrays to be filled with vigilance states and for the desired observation period
    Ltotal=[]; Ln1=[]; Ln2=[]; Lw1=[]; Lw2=[]; Lr2=[];     

    for s=1:nums %iterating over all the stimulations
        
        step=stimep(s);
%         stepend=stimend(s);
        eps=step-1:step+2*ne-1; %returns the desired period
        stateEps20sec=step-5:step-1; %for determining the initial state
        stateEps1min=step-15:step-1; 
        
        ws=wake(eps);  Ws=[Ws;ws];
%         ns=nrem(eps);  Ns=[Ns;ns];
%         rs=rems(eps);  Rs=[Rs;rs];
        ds=dexs(eps); Ds=[Ds;ds];
        %prevSleep=[prevSleep; NaN NaN];
        
        if sum(ws)==0; mints=length(ws)+1; else; mints=min(find(ws>0));  end
        mints=(mints-1)*4; % sec
        Ltotal=[Ltotal; mints];
        
        if nanmean(nrem(stateEps1min))==1; Ln1=[Ln1; s mints]; end
        if nanmean(wake(stateEps1min))==1; Lw1=[Lw1; s mints]; end
        if nanmean(nrem(stateEps20sec))==1; Ln2=[Ln2; s mints]; end
        if nanmean(wake(stateEps20sec))==1; Lw2=[Lw2; s mints]; end
        if nanmean(rems(stateEps20sec))==1; Lr2=[Lr2; s mints]; end

    end

    WakeLatencyW=[WakeLatencyW;Lw1];
    WakeLatencyN=[WakeLatencyN;Ln1];
    WakeLatencyR=[WakeLatencyR;Lr2];
    
    WakeLatencyAll(1:size(Lw1,1),3*anim-2)=Lw1(:,2); 
    WakeLatencyAll(1:size(Ln1,1),3*anim-1)=Ln1(:,2); 
    WakeLatencyAll(1:size(Lr2,1),3*anim)=Lr2(:,2);
    WakeLatencyAllv2(1:size(Lw1,1),anim)=Lw1(:,2); 
    WakeLatencyAllv2(1:size(Ln1,1),numanim*1+anim)=Ln1(:,2); 
    WakeLatencyAllv2(1:size(Lr2,1),numanim*2+anim)=Lr2(:,2);
end

Prism_WakeLatencyW_withnoWake=WakeLatencyW(:,2);
Prism_WakeLatencyN_withnoWake=WakeLatencyN(:,2);
Prism_WakeLatencyR_withnoWake=WakeLatencyR(:,2);

% NREM
samplesize_N=size(WakeLatencyN,1)
avrLatencyNREM=mean(WakeLatencyN(:,2))
seLatencyNREM=std(WakeLatencyN(:,2))/sqrt(samplesize_N)
% REM
samplesize_R=size(WakeLatencyR,1)
avrLatencyREM=mean(WakeLatencyR(:,2))
seLatencyREM=std(WakeLatencyR(:,2))/sqrt(samplesize_R)

