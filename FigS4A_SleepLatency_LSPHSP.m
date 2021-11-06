clear all
close all

path='I:\optogenetics\'; 

mousenames2=[5 15 16 17 18 19 20 21];
days2=['300318 310318';'201218 211218';'201218 211218';'201218 211218';...
    '040919 060919';'040919 060919';'040919 060919';'040919 060919'];


geno = 2;

ders=strvcat('fro','occ');
der=1; 
deri=ders(der,:); 

maxep=10800; 
epochs=1:maxep; 
epochs=epochs./900; 

epochl=4; 

zermat=zeros(1,maxep); %creates an array of zeros with 1 row and as long as the total number of epoch

fs=256;  %sampling rate (you sample 256 times in each 1second)

%%%%%%%%% number of stimulation
nstim=12;

before=0;
after=18;


ne=15; %there are 15 epochs in 1min

x=-before*ne:1:after*ne; %creates an array of epoch indices in the 2min before to after period
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

CorrN=NaN(12,11);CorrR=NaN(12,11);CorrW=NaN(12,11);
PrevSleep=NaN(12,11);

hold on 

SWApp=[];
for dd=1:2
    Ls=[];
    SleepLatencyW=[]; SleepLatencyN=[]; SleepLatencyR=[]; 
    SleepLatencyAll=NaN(12,numanim*3); SleepLatencyAllv2=NaN(12,numanim*3); 
    for anim=1:numanim %go through this loop as many times as there are animals

        mousename=[gn,num2str(mousenames(anim))];
        
        daysi=days(anim,:); is=find(isspace(daysi)); %which dates for that mouse?
        if dd==1 day=daysi(1:is(1)-1); else day=daysi(is(1):end); end         
        day(isspace(day))=[]; %stimulation day

        fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
        eval(['load ',pathstim,fn0,'.mat startend ']);

        stimep1=round(startend(:,1)./(fs*epochl)); 
        stimep2=round(startend(:,2)./(fs*epochl));

        stimep1=stimep1(1:nstim,:); 
        stimep2=stimep2(1:nstim,:); 
        
        fn1=[mousename,'-',day,'-',deri]; 
        
        eval(['load ',pathvs,fn1,'.mat w nr r w1 nr2 r3 mt ma bastend ']);

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

        wake=sum(VS([1:2 7],:));
        nrem=sum(VS(3:4,:));
        rems=sum(VS(5:6,:));

        sleep=nrem+rems; %gives total sleep state epochs

        stimep=stimep1; %start of stimulation episodes
        stimend=stimep2; %end of stimulation episodes

        %exclude values outside of the possible range
        %w/o this --> error: Array indices must be positive integers or logical values.
        out1=find(stimep-before*ne<1);  stimep(out1)=[];  
        out2=find(stimep+after*ne>maxep); stimep(out2)=[];

        nums=length(stimep); %how many stimulation episodes there are

        Ws=[]; Ns=[]; Rs=[]; Ts=[]; %creates empty arrays to be filled with vigilance states and for the desired observation period
        Ltotal=[]; 
        Ln1=NaN(nstim, 2); Ln2=NaN(nstim, 2); Lw1=NaN(nstim, 2); Lw2=NaN(nstim, 2); Lr1=NaN(nstim, 2); Lr2=NaN(nstim, 2);     

        for s=1:nums %iterating over all the stimulations

            step=stimep(s);
            stepend=stimend(s);
            eps=stepend+1:stepend+after*ne; %returns the desired period
%             stateEps1min=step-15:step-1; 
            stateEps1epoch=step-1; 
            stateEps20sec=step-5:step-1; %for determining the initial state


            ws=wake(eps);  Ws=[Ws;ws];
            ns=nrem(eps);  Ns=[Ns;ns];
            rs=rems(eps);  Rs=[Rs;rs];
            ts=sleep(eps); Ts=[Ts;ts];
            %prevSleep=[prevSleep; NaN NaN];

            if sum(ts)==0; mints=length(ts); else; mints=min(find(ts>0));  end
            mints=mints*4/60; % min
            Ltotal=[Ltotal; mints];

            if nanmean(nrem(stateEps1epoch))==1; Ln1(s,:)=[s mints]; end
            if nanmean(wake(stateEps1epoch))==1; Lw1(s,:)=[s mints]; end
            if nanmean(rems(stateEps1epoch))==1; Lr1(s,:)=[s mints]; end
            if nanmean(nrem(stateEps20sec))==1; Ln2(s,:)=[s mints]; end
            if nanmean(wake(stateEps20sec))==1; Lw2(s,:)=[s mints]; end
            if nanmean(rems(stateEps20sec))==1; Lr2(s,:)=[s mints]; end

        end
    
        Ls=[Ls Ltotal];
        
        SleepLatencyW=[SleepLatencyW;Lw1]; %Lw2
        SleepLatencyN=[SleepLatencyN;Ln1]; %Ln2
        SleepLatencyR=[SleepLatencyR;Lr1]; %Lr2

        SleepLatencyAll(1:size(Lw1,1),anim)=Lw1(:,2); 
        SleepLatencyAll(1:size(Ln1,1),numanim*1+anim)=Ln1(:,2); 
        SleepLatencyAll(1:size(Lr2,1),numanim*2+anim)=Lr1(:,2);
        SleepLatencyAllv2(1:size(Lw2,1),anim)=Lw2(:,2); 
        SleepLatencyAllv2(1:size(Ln2,1),numanim*1+anim)=Ln2(:,2); 
        SleepLatencyAllv2(1:size(Lr2,1),numanim*2+anim)=Lr2(:,2);
    end


    LWs=SleepLatencyW(:,2);
    LNs=SleepLatencyN(:,2);
    LRs=SleepLatencyR(:,2);
    aPrism_SleepLatencyW_withnoSleep=reshape(LWs,[],anim);
    aPrism_SleepLatencyN_withnoSleep=reshape(LNs,[],anim);
    aPrism_SleepLatencyR_withnoSleep=reshape(LRs,[],anim);
    aPrism_Ls=Ls;
    
    LN_bin=reshape(LNs,3,[]); LN_bin=nanmean(LN_bin); LN_bin=reshape(LN_bin,[],anim);
    
    Ls_bin=reshape(Ls,3,4,numanim); Ls_bin=nanmean(Ls_bin,1); 
    Ls_bin=reshape(Ls_bin,[],size(Ls_bin,3),1);
    
    aPrism_LN_bin=LN_bin; % the animal was in NREM sleep when stimulation happens
    aPrism_Ls_bin=Ls_bin; % the animal was in NREM/REM sleep when stimulation happens
    
    if dd==1; aP_1st=[];aP_2nd=[];aP_3rd=[];aP_4th=[]; end
    aP_1st(:,dd)=Ls_bin(1,:);
    aP_2nd(:,dd)=Ls_bin(2,:);
    aP_3rd(:,dd)=Ls_bin(3,:);
    aP_4th(:,dd)=Ls_bin(4,:);    
    
    % NREM
    outN=find(SleepLatencyN==18); %SleepLatencyN(outN)=[];
    numdidnotreturntosleep_NREM=size(outN,1)
    samplesize_N=size(SleepLatencyN,2)
    avrLatencyNREM=nanmean(SleepLatencyN)
    seLatencyNREM=nanstd(SleepLatencyN)/sqrt(samplesize_N)
    
    % REM
    outR=find(SleepLatencyR==18); %SleepLatencyR(outR)=[];
    numdidnotreturntosleep_REM=size(outR,1)
    samplesize_R=size(SleepLatencyR,2)
    avrLatencyREM=nanmean(SleepLatencyR)
    seLatencyREM=nanstd(SleepLatencyR)/sqrt(samplesize_R)

    % Wake
    outW=find(SleepLatencyW==18); %SleepLatencyW(outW)=[];
    numdidnotreturntosleep_Wake=size(outW,1)
    samplesize_W=size(SleepLatencyW,2)
    avrLatencyWake=nanmean(SleepLatencyW)
    seLatencyWake=nanstd(SleepLatencyW)/sqrt(samplesize_W)
end