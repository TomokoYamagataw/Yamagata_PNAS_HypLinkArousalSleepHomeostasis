
clear all
% close all

path='I:\optogenetics\'; % the working directory with the data folders
% path='D:\VVlab\Opto\SDstim\'; % the working directory with the data folders


%Data from 24h stimulation of 2min every 20min with 10Hz blue light
mousenames1=[1 2 3 4 5 6 7 8]; % names of GFP ctrl mice indicated in the file name 
days1=['260218';'120418';'260618';'040718';'061018';'061018';'061018';'070119'];

% 
% 
% % 1Hz
% mousenames2=[15 16 17 18 19 20 21];
% days2=['241218';'241218';'241218';'190819';...
%     '190819';'190819';'190819'];
% 
% % 2Hz
% mousenames3=[15 16 17 18 19 20 21];
% days3=['251218';'251218';'251218';'200819';...
%     '200819';'200819';'200819'];
% 
% % 5Hz
% mousenames4=[15 16 17 18 19 20 21];
% days4=['261218';'261218';'261218';'210819';...
%     '210819';'210819';'210819'];
% 
% % 10Hz
% mousenames5=[15 16 17 18 19 20 21];
% days5=['191218';'191218';'191218';'160819';...
%     '160819';'160819';'160819'];
% 
% % Sham
% mousenames6=[15 16 17 18 19 20 21]; 
% days6=['181218';'181218';'181218';...
%     '100819';'100819';'070819';'100819']; 



% 10Hz
%%%%%%%% 10Hz 2min LPO
mousenames7=[5 15 16 17 18 19 20 21];
days7=['030418';'191218';'191218';'191218';...
    '160819';'160819';'180919';'160819'];

%%%%%%%% 10Hz 2min nonLPO
mousenames8=[1 6 7 8 9 12]; %non-LPO %GDCh8 '130518 140518'
days8=['140218';'120418';'130518';'150618';...
    '150618';'090618']; 


dd=0; % 24h: 0 %Light:1 or Dark:2

geno=1; 

tfstim=1;
tnstim=72; 

ne=15; %there are 15 epochs in 1min 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defining all the derivations
ders=strvcat('fro','occ');
der=1; 
deri=ders(der,:);

maxep=21600; %the total no of 4s epochs in 24h period 
epochs=1:maxep; %creates an array with 1 row and columns from 1 to 21600(maxep)
epochs=epochs./900; %edits x such that each value is divided by 900 to convert the values to hours 
            %15 epochs per 1min and 60min per 1h so 900 epochs per 1h
            %%% ./ = right-array division by dividing each element of A by the corresponding element of B
fs=256; 
epochl=4; %states how long the epoch is
zermat=zeros(1,maxep); %creates an array of zeros with 1 row and as long as the total number of epoch



pathvs1=[path,'outputVSgfp\']; 
pathvs2=[path,'outputVSchr\']; 
pathstim=[path,'STIMs\']; 

if geno==1 %geno=1 ctrl group
    mousenames=mousenames1;days=days1; pathvs=pathvs1; gn='GFP'; %firststims=firststims1;
elseif geno==2
    mousenames=mousenames2; days=days2; pathvs=pathvs2; gn='GDCh';  %firststims=firststims2;
elseif geno==3
    mousenames=mousenames3; days=days3; pathvs=pathvs2; gn='GDCh';
elseif geno==4
    mousenames=mousenames4; days=days4; pathvs=pathvs2; gn='GDCh'; 
elseif geno==5
    mousenames=mousenames5; days=days5; pathvs=pathvs2; gn='GDCh';
elseif geno==6
    mousenames=mousenames6; daybs=days5; days=days6; pathvs=pathvs2; gn='GDCh'; %Sham stim timing: 10Hz
elseif geno==7
    mousenames=mousenames7; days=days7; pathvs=pathvs2; gn='GDCh'; %Sham stim timing: 10Hz
elseif geno==8
    mousenames=mousenames8; days=days8; pathvs=pathvs2; gn='GDCh'; %Sham stim timing: 10Hz
end
      

vsname=strvcat('Wake','NREM','REM'); 

numanim=length(mousenames); 

WallMice=[];NallMice=[];RallMice=[];MallMice=[];
WallW=[];WallN=[];WallR=[];
NallW=[];NallN=[];NallR=[];
RallW=[];RallN=[];RallR=[];
WNRnum=[];

 Lmeans=[];
for anim=1:numanim
    mousename=[gn,num2str(mousenames(anim))]; mousename(isspace(mousename))=[];
    day=days(anim,:); day(isspace(day))=[]; 
    
    if geno==6
        dayb=daybs(anim,:); dayb(isspace(dayb))=[];
        fn0=[mousename,'_',dayb,'_stim']; 
        eval(['load ',pathstim,fn0,'.mat startend ']);
    else        
        fn0=[mousename,'_',day,'_stim']; 
        eval(['load ',pathstim,fn0,'.mat startend ']);
    end

    
    stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
    stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2) 

    outs=find((stimep2-stimep1)<25);
    stimep1(outs)=[];stimep2(outs)=[];
    
    if dd==1
        stimep1(stimep2>=10800)=[]; stimep2(stimep2>=10800)=[];
    elseif dd==2
        stimep1(stimep2<10800)=[]; stimep2(stimep2<10800)=[];
    else
    end

    fn=[mousename,'-',day,'-',deri];  
    eval(['load ',pathvs,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);%
    
    VS=zeros(8,maxep); 
    VS(1,w)=1; VS(2,w1)=1; VS(3,nr)=1; VS(4,nr2)=1; VS(5,r)=1; VS(6,r3)=1; VS(7,mt)=1; VS(8,ma)=1;
    clear w nr r w1 nr2 r3 mt ma bastend ; %clears the variables for the next loop 
    
    wake=sum(VS(1:2,:));     
    nrem=sum(VS(3:4,:));     
    rem=sum(VS(5:6,:)); 
    move=sum(VS(7:8,:));
    
    numstim=size(stimep1,1) %how many stimulation episodes there are 

% remove bag of continuous movement    
%     if find(sum(VS(:,stimep2))<1)>0 
%         stimep2(find(sum(VS(:,stimep2))<1))
%     end
    
    Ws=[]; Ns=[]; Rs=[]; % Ms=[];
    for s=1:numstim %iterating over all the stimulations
        eps=stimep1(s)-1; 
        eps2=stimep2(s); 
        
        ms=move(eps);  ms2=move(eps2);  
        ws=wake(eps);  ws2=wake(eps2);    
        Ws=[Ws;ws+ms ws2+ms2]; 

        ns=nrem(eps); ns2=nrem(eps2);
        Ns=[Ns;ns ns2];     
        
        rs=rem(eps); rs2=rem(eps2); 
        Rs=[Rs;rs rs2];
    end

   %mean(VS) will average the observation period across all stimulations
   %allMice arrays combine the mean data from all the mice 
   WallMice=[WallMice;nanmean(Ws)]; %adds the Ws of this mouse to the group array 
   NallMice=[NallMice;nanmean(Ns)];    
   RallMice=[RallMice;nanmean(Rs)];       
%    MallMice=[MallMice;nanmean(Ms)];

   outw=nansum(Ws(:,1),2)<1; 
   Wsw=Ws; Wsw(outw,:)=[]; WallW=[WallW;nanmean(Wsw)];
   Nsw=Ns; Nsw(outw,:)=[]; NallW=[NallW;nanmean(Nsw)];
   Rsw=Rs; Rsw(outw,:)=[]; RallW=[RallW;nanmean(Rsw)];
   
   outn=nansum(Ns(:,1),2)<1; 
   Wsn=Ws; Wsn(outn,:)=[]; WallN=[WallN;nanmean(Wsn)];
   Nsn=Ns; Nsn(outn,:)=[]; NallN=[NallN;nanmean(Nsn)];
   Rsn=Rs; Rsn(outn,:)=[]; RallN=[RallN;nanmean(Rsn)];   
   
   outr=nansum(Rs(:,1),2)<1; 
   Wsr=Ws; Wsr(outr,:)=[]; if size(Wsr,1)>1; WallR=[WallR;nanmean(Wsr)]; else;  WallR=[WallR;Wsr]; end
   Nsr=Ns; Nsr(outr,:)=[]; if size(Nsr,1)>1; NallR=[NallR;nanmean(Nsr)]; else;  NallR=[NallR;Nsr]; end
   Rsr=Rs; Rsr(outr,:)=[]; if size(Rsr,1)>1; RallR=[RallR;nanmean(Rsr)]; else;  RallR=[RallR;Rsr]; end
   

   WNRnum=[WNRnum; numstim size(Wsw,1) size(Wsn,1) size(Wsr,1)];
end

aPrism_Wake=WallMice';
aPrism_NREM=NallMice';
aPrism_REM=RallMice';
aPrism_VSprobability=[nanmean(WallMice)*100; nanmean(NallMice)*100; nanmean(RallMice)*100];

aPrism_contW=[WallW' NallW' RallW'];
aPrism_contN=[WallN' NallN' RallN'];
aPrism_contR=[WallR' NallR' RallR'];

aPrism_initial_state=WNRnum'
aPrism_initial_state_ratio=WNRnum./WNRnum(:,1); 
aPrism_initial_state_ratio=aPrism_initial_state_ratio';

aPrism_initial_state_probability=WNRnum./WNRnum(:,1)*100; 
aPrism_initial_state_probability=aPrism_initial_state_probability';

% if geno==7; WNRnum=WNRnum([1 5 6 7 8 9 10 11], :); end
WNRnumsum=sum(WNRnum,1);
% aPrism_initial_state_ratio=WNRnumLPOsum./WNRnumLPOsum(:,1); aPrism_initial_state_ratio=aPrism_initial_state_ratio'
aPrism_initial_state_probabilities=WNRnumsum./WNRnumsum(:,1)*100; 
aPrism_initial_state_probabilities=aPrism_initial_state_probabilities'


% aPrism_Wake=WallMice';
% aPrism_NREM=NallMice';
% aPrism_REM=RallMice';
% 
% aPrism_contW=[WallW' NallW' RallW'];
% aPrism_contN=[WallN' NallN' RallN'];
% aPrism_contR=[WallR' NallR' RallR'];


