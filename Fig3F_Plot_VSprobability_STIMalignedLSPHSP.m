
clear all
close all

path='I:\optogenetics\'; % the working directory with the data folders

% mousenames1=[1 2 3 4 5 6 7 8]; % names of GFP ctrl mice indicated in the file name 
% days1=['260218';'120418';'260618';'040718';'061018';'061018';'061018';'070119'];

% 
% mousenames1=[5 15 16 17 18 19 21];
% days1=['300318 310318';'201218 211218';'201218 211218';'201218 211218';...
%     '040919 060919';'040919 060919';'040919 060919'];

mousenames2=[5 15 16 17 18 19 20 21];
days2=['300318 310318';'201218 211218';'201218 211218';'201218 211218';...
    '040919 060919';'040919 060919';'040919 060919';'040919 060919'];

% % 10Hz
% mousenames5=[15 16 17 18 19 20 21];
% days5=['191218';'191218';'191218';'160819';...
%     '160819';'160819';'160819'];

% % Sham
% mousenames6=[15 16 17 18 19 20 21]; 
% days6=['181218';'181218';'181218';...
%     '100819';'100819';'070819';'100819']; 



geno=2; %states which animal group (genotype) I want to look at 

tfstim=1; %2 %1
tnstim=12; %6 %5

ne=15; %there are 15 epochs in 1min 
before=2;  %we want to analyze the 2min before stimulus onset and 2min after
after=8;   %gives me 2min before, 2min of stim, and the 4min after 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defining all the derivations
ders=strvcat('fro','occ');
der=1; 
deri=ders(der,:);

maxep=10800; %the total no of 4s epochs in 24h period 
epochs=1:maxep; %creates an array with 1 row and columns from 1 to 21600(maxep)
epochs=epochs./900; %edits x such that each value is divided by 900 to convert the values to hours 
            %15 epochs per 1min and 60min per 1h so 900 epochs per 1h
            %%% ./ = right-array division by dividing each element of A by the corresponding element of B
fs=256; 
epochl=4; %states how long the epoch is
zermat=zeros(1,maxep); %creates an array of zeros with 1 row and as long as the total number of epoch

x=-before*ne:1:after*ne; %creates an array of epoch indices in the 2min before to after period
x=x*epochl; %converts the epoch indices to seconds 


pathvs1=[path,'outputVSgfp\']; %defines the folder with output files of vigilance states for mice1 (GFP ctrl)
pathvs2=[path,'outputVSchr\']; %defines the folder with output files of vigilance states for mice2 (ChR2 exp)
pathstim=[path,'STIMs\']; %gives the folder with stimulation inf

if geno==1 %geno=1 ctrl group
    mousenames=mousenames1;days=days1; pathvs=pathvs1; gn='GFP'; %firststims=firststims1;
elseif geno==2
    mousenames=mousenames2; days=days2; pathvs=pathvs2; gn='GDCh';  %firststims=firststims2;
end

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

numanim=length(mousenames); %how many animals were recorded

%wake/NREM/REM/Sleep/Dex state of all mice 
      
% % this loop goes though all the animals
% % for each 2min before and after we look at amount of sleep
% % Question: how efficient is stimulation in waking up the mice
%  Ls=NaN(numanim, 5);

for dd=1:2

WallMice=[];NallMice=[];RallMice=[];MallMice=[];
Lmeans=[];
    for anim=1:numanim %go through this loop as many times as there are animals
        mousename=[gn,num2str(mousenames(anim))]; mousename(isspace(mousename))=[];
        
        daysi=days(anim,:); is=find(isspace(daysi)); %which dates for that mouse?
        if dd==1 day=daysi(1:is(1)-1); else day=daysi(is(1):end); end         
        day(isspace(day))=[]; %stimulation day
        
        fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
        eval(['load ',pathstim,fn0,'.mat startend ']);

        stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
        stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2) 

        outs=find((stimep2-stimep1)<25);stimep1(outs)=[];stimep2(outs)=[];

        fn=[mousename,'-',day,'-',deri]; %makes the full file name for the vigilance state and derivation output desired

        eval(['load ',pathvs,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);%

        VS=zeros(7,maxep); 
        VS(1,w)=1; VS(2,w1)=1; VS(3,nr)=1; VS(4,nr2)=1; VS(5,r)=1; VS(6,r3)=1; VS(7,mt)=1;
        clear w nr r w1 nr2 r3 mt ma bastend ; %clears the variables for the next loop 

        wake=sum(VS(1:2,:));     
        nrem=sum(VS(3:4,:));     
        rem=sum(VS(5:6,:)); 
        move=VS(7,:);


        stimep=stimep1(tfstim:tfstim+tnstim-1); %start of stimulation episodes
%         stimep=stimep1; %start of stimulation episodes
        stimep4=stimep2(tfstim:tfstim+tnstim-1);   %after stimulation episodes

        numstim=size(stimep,1) %how many stimulation episodes there are 

        Ws=[]; Ns=[]; Rs=[]; Ms=[]; LmatN=[]; LmatR=[]; LmatW=[];%creates empty arrays to be filled 
        for s=1:numstim %iterating over all the stimulations
            step=stimep(s); 
            eps=step-before*ne:step+after*ne-1; %returns the desired period
            eps1=step-before*ne:step-1; %returns the desired period
            eps2=step+1:step+2*ne; %returns after stim

            if min(eps)<1; 
            else

            ms=move(eps);  
    %         Ms=[Ms;ms]; 
            ws=wake(eps);  %for this stimulation period only
            Ws=[Ws;ws+ms];  %for each of the 72 stims  
            LWmat=ws(1,31:60);
            LW1=min(find(LWmat>0)); 


            ns=nrem(eps); 
            Ns=[Ns;ns]; 
            LNmat=ns(1,31:60);
            LN1=min(find(LNmat>0));        

            rs=rem(eps); 
            Rs=[Rs;rs];
            LRmat=rs(1,31:60);
            LR1=min(find(LRmat>0));

            if LW1<60; 
            else LW1=60; 
                if LN1<60; else; LN1=60; end
            end

    %         Ls(numstim*(anim-1)+s,:)=[mousenames(anim) s LN1 LR1 LW1];


            LmatN=[LmatN LN1]; LmatR=[LmatR LR1]; LmatW=[LmatW LW1]; %LmatM=[LmatM LM1];
            end

        end

       %mean(VS) will average the observation period across all stimulations
       %allMice arrays combine the mean data from all the mice 
       WallMice=[WallMice;nanmean(Ws)]; %adds the Ws of this mouse to the group array 
       NallMice=[NallMice;nanmean(Ns)];    
       RallMice=[RallMice;nanmean(Rs)];       
    %    MallMice=[MallMice;nanmean(Ms)];

       LmeanN=nanmean(LmatN); LmeanR=nanmean(LmatR);  LmeanW=nanmean(LmatW);
       Lmeans=[Lmeans; mousenames(anim) LmeanN LmeanR LmeanW];

%        outw=nansum(Ws(:,16:30),2)<15; 
%        Wsw=Ws; Wsw(outw,:)=[]; WallW=[WallW;nanmean(Wsw)];
%        Nsw=Ns; Nsw(outw,:)=[]; NallW=[NallW;nanmean(Nsw)];
%        Rsw=Rs; Rsw(outw,:)=[]; RallW=[RallW;nanmean(Rsw)];
%        outn=nansum(Ns(:,16:30),2)<15; 
%        Wsn=Ws; Wsn(outn,:)=[]; WallN=[WallN;nanmean(Wsn)];
%        Nsn=Ns; Nsn(outn,:)=[]; NallN=[NallN;nanmean(Nsn)];
%        Rsn=Rs; Rsn(outn,:)=[]; RallN=[RallN;nanmean(Rsn)];   
%        outr=nansum(Rs(:,26:30),2)<5; 
%        Wsr=Ws; Wsr(outr,:)=[]; WallR=[WallR; nanmean(Wsr)]; 
%        Nsr=Ns; Nsr(outr,:)=[]; NallR=[NallR;nanmean(Nsr)];
%        Rsr=Rs; Rsr(outr,:)=[]; RallR=[RallR;nanmean(Rsr)]; 
%        nw=size(Wsw,1); nn=size(Wsn,1);  nr=size(Wsr,1); 
%        WNRnum=[WNRnum; nw nn nr];
    end

    aPrism_Wake=WallMice';
    aPrism_NREM=NallMice';
    aPrism_REM=RallMice';
    
    plot(nanmean(NallMice),'o-','LineWidth',2)
    hold on
    
%     aPrism_contW=[WallW' NallW' RallW'];
%     aPrism_contN=[WallN' NallN' RallN'];
%     aPrism_contR=[WallR' NallR' RallR'];


    aPrism_Latency=Lmeans;
    LatencyW=[mousenames; WallMice(:,30:61)'];
    aPrism_Latency=LatencyW;

   
end

