
clear all
close all

path='E:\Optoinh\'; 

geno=2; %1, 2
dd=2; % 24h: 0 %Light:1 or Dark:2

% % % % Arch
mousenames1=[1 2 4 5 8 9]; % names of mice indicated in the file name 
days1=['160421';'150421';'160421';'140521';'140521';'150521'];

% % % % GFP
mousenames2=[1 2 4 5]; 
days2=['150421';'150421';'150521';'140521'];


tfstim=1;
tnstim=48; 

ne=15; %there are 15 epochs in 1min 
before=2;  %we want to analyze the 2min before stimulus onset and 2min after
after=8;   %gives me 2min before, 5min of stim, and the 2min after 

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

x=-before*ne:1:after*ne; %creates an array of epoch indices in the 2min before to after period
x=x*epochl; %converts the epoch indices to seconds 

pathvs=[path,'outputVS\']; %defines the folder with output files of vigilance states for mice2 (ChR2 exp)
pathstim=[path,'STIMs\']; %gives the folder with stimulation inf

if geno==1 %geno=1 ctrl group
    mousenames=mousenames1;days=days1; pathvs=pathvs; gn='Ar'; %firststims=firststims1;
elseif geno==2
    mousenames=mousenames2; days=days2; pathvs=pathvs; gn='Gf';  %firststims=firststims2;
end
      

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

numanim=length(mousenames); %how many animals were recorded

WallMice=[];NallMice=[];RallMice=[];MallMice=[];
SPallMice=NaN*ones(150,121,numanim);
WallW=[];WallN=[];WallR=[];
NallW=[];NallN=[];NallR=[];
RallW=[];RallN=[];RallR=[];
WNRnum=[];
%wake/NREM/REM/Sleep/Dex state of all mice 

 Lmeans=[];
for anim=1:numanim %go through this loop as many times as there are animals
    mousename=[gn,num2str(mousenames(anim))]; mousename(isspace(mousename))=[];
    day=days(anim,:); day(isspace(day))=[]; %which days of recordings correspond to that mouse (both baseline and stim)
    %     tfs=firststims(anim,:);
    
    fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
    eval(['load ',pathstim,fn0,'.mat startend ']);
    
    stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
    stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2) 

    outs=find((stimep2-stimep1)<70);stimep1(outs)=[];stimep2(outs)=[];
    if dd==1
        stimep1(stimep2>=10800)=[]; stimep2(stimep2>=10800)=[];
    elseif dd==2
        stimep1(stimep2<10800)=[]; stimep2(stimep2<10800)=[];
    else
    end
    
    fn=[mousename,'-',day,'-',deri]; %makes the full file name for the vigilance state and derivation output desired
    
    eval(['load ',pathvs,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);%
    
    VS=zeros(7,maxep); 
    VS(1,w)=1; VS(2,w1)=1; VS(3,nr)=1; VS(4,nr2)=1; VS(5,r)=1; VS(6,r3)=1; VS(7,mt)=1;
    clear w nr r w1 nr2 r3 mt ma bastend ; %clears the variables for the next loop 
    
    wake=sum(VS(1:2,:));     
    nrem=sum(VS(3:4,:));     
    rem=sum(VS(5:6,:)); 
    move=VS(7,:);
    art=sum(VS([2 4 6 7],:)); 
        
%     stimep=stimep1(tfstim:tfstim+tnstim-1); %start of stimulation episodes
    stimep=stimep1; %start of stimulation episodes
        
    numstim=size(stimep,1) %how many stimulation episodes there are 
    
    Ws=[]; Ns=[]; Rs=[]; Ms=[]; LmatN=[]; LmatR=[]; LmatW=[]; 
    SPs=NaN*ones(150,121,numstim); SPmat=[];
    
    for s=1:numstim %iterating over all the stimulations
        step=stimep(s); 
        eps=step-before*ne:step+after*ne-1; %returns the desired period
        eps1=step-before*ne:step-1; %returns the desired period
        eps2=step+1:step+2*ne; %returns the desired period
        
        if min(eps)<1; 
        elseif max(eps)>maxep; 
        else
        
        ms=move(eps);  
%         Ms=[Ms;ms]; 
        ws=wake(eps);  %for this stimulation period only
        Ws=[Ws;ws+ms];  %for each of the 48 stims  
        LWmat=ws(1,31:105);
%         LW1=min(find(LWmat>0)); 


        ns=nrem(eps); 
        Ns=[Ns;ns]; 
        LNmat=ns(1,31:105);
%         LN1=min(find(LNmat>0));        
        
        rs=rem(eps); 
        Rs=[Rs;rs];
        LRmat=rs(1,31:105);
%         LR1=min(find(LRmat>0));
        
        spectrna=spectr;
        art=logical(art);
        spectrna(art,:)=NaN;
        spectrm=spectr(eps,:); 
        SPs(:,:,s)=spectrm;

        end

    end

   WallMice=[WallMice;nanmean(Ws)]; %adds the Ws of this mouse to the group array 
   NallMice=[NallMice;nanmean(Ns)];    
   RallMice=[RallMice;nanmean(Rs)];       
   SPallMice(:,:,anim)=nanmean(SPs,3);


end

   
freq = 0.5:0.25:30;
timeInS = -116:4:120; 

figure;

toPlot = {occStimGFP,occStimArch};
minMax = cell2mat(cellfun(@(x)x(:,1:32),toPlot,'un',0)');
 minMax = [min(minMax(:)) max(minMax(:))];
for plotCnt = 1:2
    subplot(2,2,plotCnt)
    curSpec =smoothdata(squeeze(toPlot{plotCnt}(:,3:121)),1,'movmedian',5);
    pcolor(timeInS,freq,curSpec'); shading flat; colorbar();
    xlabel('time to stim [s]'); ylabel('hz');colormap('jet'); caxis(minMax)
end



aPrism_Wake=WallMice';
aPrism_NREM=NallMice';
aPrism_REM=RallMice';
aPrism_Spectr=SPallMice';

% aPrism_contW=[WallW' NallW' RallW']
% aPrism_contN=[WallN' NallN' RallN']
% aPrism_contR=[WallR' NallR' RallR']

% 
% aPrism_Latency=Lmeans
% LatencyW=[mousenames; WallMice(:,30:105)']
% aPrism_Latency=LatencyW



