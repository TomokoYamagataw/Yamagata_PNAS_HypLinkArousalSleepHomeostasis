
clear all
close all

path='I:\optogenetics\'; % the working directory with the data folders
% path='D:\VVlab\Opto\SDstim\'; % the working directory with the data folders


geno=1; %states which animal group (genotype) I want to look at 


%Data from 24h stimulation of 2min every 20min with 10Hz blue light
%%%%%%%% GFP
mousenames1=[1 2 3 4 5 6 7 8]; % names of GFP ctrl mice indicated in the file name 
days1=['120518';'050518';'180718';'180718';'111018';'221018';'221018';'150119'];
% savefilename='VSspec_sedation_fro_stimeffect_percentage_GFP_stim1to4'


% %%%%%%% LPO fro
% mousenames2=[1 5 6 7 8 9 12 16 18 19 20 21]; %fro
% days2=['010518';'010518';'050518';'150518';'220618';'220618';'220618';'150119';'220819';'220819';'220819';'220819']; 
% savefilename='VSspec_sedation_fro_stimeffect_percentage_GDCh_all_stim1to4'
% der=1; 
mousenames2=[5 16 18 19 20 21]; %fro
days2=['010518';'150119';'220819';'220819';'220819';'220819']; 
% savefilename='VSspec_sedation_fro_stimeffect_percentage_GDCh_LPO_stim1to4'

%%%%%%% nonLPO
mousenames4=[1 6 7 8 9 12]; %fro
days4=['010518';'050518';'150518';'220618';'220618';'220618']; 
der=1; 

% %%%%%%% GDCh occ
% mousenames2=[1 5 6 7 9 12 16 17 18 19 20 21];  %occ
% days2=['010518';'010518';'050518';'150518';'220618';'220618';'150119';'150119';'220819';'220819';'220819';'220819']; 
% savefilename='VSspec_sedation_occ_stimeffect_percentage_GDCh_all_stim1to4'
% der=2; 




% der=1; %%%%%%%%% GDCh grouping
% mousenames2=[5 16 19 21]; 
% days2=['010518';'150119';'220819';'220819']; 
% savefilename='VSspec_sedation_fro_stimeffect_percentage_GDCh5_16_19_21_stim1to4'
% 
% mousenames2=[6 12 18 20]; 
% days2=['050518';'220618';'220819';'220819']; 
% savefilename='VSspec_sedation_fro_stimeffect_percentage_GDCh6_12_18_20_stim1to4'

% mousenames2=[1 7 9 18]; 
% days2=['010518';'150518';'220618';'220819']; 
% savefilename='VSspec_sedation_fro_stimeffect_percentage_GDCh1_7_9_18_stim1to4'



%%%%%%%%%%%
% mousenames2=[1 6 7 9 12 18 20]; 
% days2=['010518';'050518';'150518';'220618';'220618';'220819';'220819']; 
% savefilename='VSspec_sedation_fro_stimeffect_percentage_GDCh1_6_7_9_12_18_20_stim1to4'

% mousenames2=[5 6 16 18 19 21]; 
% days2=['010518';'050518';'150119';'220819';'220819';'220819']; 
% savefilename='VSspec_sedation_fro_stimeffect_percentage_GDCh5_6_16_18_19_21_stim1to4'
%
% mousenames2=[1 7 9 12 20]; 
% days2=['010518';'150518';'220618';'220618';'220819']; 
% savefilename='VSspec_sedation_fro_stimeffect_percentage_GDCh1_7_9_12_20_stim1to4'





tfstim=1;
tnstim=4; %%%% stim per 225, jitter 22.5, stimlength=30, - 4times: 675+jitter=730, 5times: 900+jitter=955; 6times: 1125+jitter=1180, 8times: 1800+jitter=1630, 16times: 1800+jitter=3430

ne=15; %there are 15 epochs in 1min 
before=2;  %we want to analyze the 2min before stimulus onset and 2min after
after=2;   %gives me 2min before, 2min of stim, and the 2min after 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defining all the derivations
ders=strvcat('fro','occ');
deri=ders(der,:);

freqband=0:0.25:30;
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
pathout=[path,'Dex\']; mkdir(pathout)
pathfig=[path,'Figures_Dex\']; mkdir(pathout)

if geno==1 %geno=1 looks at the ctrl group
    mousenames=mousenames1;days=days1; pathvs=pathvs1; gn='GFP'; %firststims=firststims1;
elseif geno==2 % looks at the exp group
    mousenames=mousenames2; days=days2; pathvs=pathvs2; gn='GDCh';  %firststims=firststims2; 
elseif geno==4 
    mousenames=mousenames4; days=days4; pathvs=pathvs2; gn='GDCh';  %firststims=firststims2; 
end
      

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

numanim=length(mousenames); %how many animals were recorded

WallMice=[];NallMice=[];RallMice=[];TallMice=[];MallMice=[];DallMice=[];
%wake/NREM/REM/Sleep/Dex state of all mice 
      
% % this loop goes though all the animals
% % for each 2min before and after we look at amount of sleep
% % Question: how efficient is stimulation in waking up the mice
Diffs=[];Msp1s=[]; Msp2s=[];
for anim=1:numanim %go through this loop as many times as there are animals
    mousename=[gn,num2str(mousenames(anim))]; mousename(isspace(mousename))=[];
    day=days(anim,:); day(isspace(day))=[]; %which days of recordings correspond to that mouse (both baseline and stim)
%     tfs=firststims(anim,:);
    
    fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
    eval(['load ',pathstim,fn0,'.mat startend ']);
    %eval(EXPRESSION) evaluates the code in EXPRESSION that is a
    %string scalar 
    % --> loads the startend variable
    
    stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
    stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2) 

    outs=find((stimep2-stimep1)<25);stimep1(outs)=[];stimep2(outs)=[];
    
    fn=[mousename,'-',day,'-',deri]; %makes the full file name for the vigilance state and derivation output desired
    
    eval(['load ',pathvs,fn,'.mat spectr w nr r s w1 nr2 r3 s4 mt ma bastend -mat']);%

    VS=zeros(9,maxep); 
    VS(1,w)=1; VS(2,w1)=1; VS(3,nr)=1; VS(4,nr2)=1; VS(5,r)=1; VS(6,r3)=1; VS(7,mt)=1; VS(8,s)=1; VS(9,s4)=1;
    wake=sum(VS(1:2,:)); 
    nrem=sum(VS(3:4,:));
    rems=sum(VS(5:6,:));
    move=VS(7,:);
    dex=sum(VS(8:9,:));
    
    
    cleanSp=spectr;
    cleanSp([w1; nr2; r3; s4; mt; ma],:)=NaN;
    clear w nr r s w1 nr2 r3 s4 mt ma bastend ; %clears the variables for the next loop 
    
    Delta=nanmean(cleanSp(:,3:17),2);

    
    S=nrem+rems; %gives total sleep state epochs
        
    stimep=stimep1(tfstim:tfstim+tnstim-1); %start of stimulation episodes
    
    numstim=size(stimep,1) %how many stimulation episodes there are 
    
    meanSps1=[]; meanSps2=[]; Ds=[]; Ms=[]; 
    Sps1=zeros(30,121,numstim);Sps2=zeros(30,121,numstim);%creates empty arrays to be filled 
    for s=1:numstim %iterating over all the stimulations
        
        step=stimep(s); 
        
        eps1=step-before*ne:step-1; %returns the desired period
        eps2=step+1:step+2*ne; %returns the desired period
        
        sp1=cleanSp(eps1,:); sp2=cleanSp(eps2,:);
        Sps1(:,:,s)=sp1; Sps2(:,:,s)=sp2;
        meansp1=nanmean(sp1,1); meansp2=nanmean(sp2,1);
        meanSps1=[meanSps1;meansp1]; meanSps2=[meanSps2;meansp2];
        
    end
    Msp1=mean(meanSps1,1); Msp2=mean(meanSps2,1);
    Msp1s=[Msp1s; Msp1]; Msp2s=[Msp2s; Msp2];
    
    diffmeanSps=meanSps2./meanSps1*100;
    
    mousename
    aPrism_meanSps1=[meanSps1'];
    aPrism_meanSps2=[meanSps2'];
    aPrism_diffmeanSps=[diffmeanSps'];
   
    
    diffmeanSps(:,1:2)=NaN; %1-2:0-0.5, 
    diffmeanSps(:,39:43)=NaN; %39: 9.5, 40:9.75, 41:10, 42:10.25 43:10.5
    diffmeanSps(:,80:82)=NaN; %79: 19.5, 80:19.75, 81:20, 82:20.25 83:20.5
    diffmeanSps(:,121)=NaN; %120:29.75 121:30,

    Diff=nanmean(diffmeanSps);
    Diffs=[Diffs; Diff];

    

%     meanSps1(:,1:2)=NaN; meanSps2(:,1:2)=NaN;    
    
%     figure
%     subplot(1,2,1); semilogy(meanSps1'); legend('pre1','pre2','pre3','pre4');title (mousename);ylim([10^0 8*10^4])
%     hold on
%     subplot(1,2,2); semilogy(meanSps2'); legend('post1','post2','post3','post4');title (mousename);ylim([10^0 8*10^4])
%     
% 
%     figure
%     plot(diffmeanSps'); 
%     hold on;
%     legend('stim1','stim2','stim3','stim4'); 
%     title(mousename);
%     line([0 121],[100 100],'color','k');
%     ylim([0 500]);
%     
%     saveas(gcf,[pathfig,'\',mousename,'-',day,'-',deri,'_VSspec_sedation_stimeffect'],'tiff')

%     spDpre=mean(Sps1,3); 
%     spDpost=mean(Sps2,3);
    
   
end
    
aPrism_absSp1=[Msp1s'];
aPrism_absSp2=[Msp2s'];
aPrism_Diff_mean=Diffs';
    
figure
SE=std(Diffs)./sqrt(size(Diffs,1));
errorbar(freqband,mean(Diffs),SE)
line([0 30],[100 100],'color','k');
ylim([0 300]);
% saveas(gcf,[pathfig,savefilename],'tiff')

DiffsDelta=mean(Diffs(:,3:17),2);
aPrism_DeltaEffect=[mousenames' DiffsDelta]




% data=DallMice; %states what I want to plot
% 
% data_mean = mean(data,1); %1 specifies the dimensions, in this case across rows
% data_std  = std(data,0,1); %0 is weight (default=0) and 1=dimensions (across rows)
% err = (data_std./sqrt(size(data,1))); %SEM error 
% 
% options.color_area = [243 169 114]./255;
% options.color_line = [236 112  22]./255;
% 
% x2 = [x, fliplr(x)]; %creates a vector of x and then flipped x again, required for the error bars
% inBetween = [data_mean+err, fliplr(data_mean-err)]; %error bar vectors 
% figure
% patch(x2, inBetween, options.color_area,'EdgeColor', 'none', 'FaceAlpha', 0.5); %fills the area between error bars
% hold on 
% plot(x,data_mean,'Color', options.color_line, 'LineWidth', 2) %plots the single average for all mice 
% xticks(-120:60:240)
% yticks(0:0.2:1)
% xlim([-120 240])
% ylim([0 1])
% xlabel('Time [s]')
% ylabel('Sedation Probability')
% line([0 0], [0 1], 'Color', 'b', 'LineWidth', 1.5)
% title(['Sedation probaility ',gn,' Average ',num2str(tnstim),'stimulations from the ',num2str(tfstim),'nd stimulus'])
% 
% fn8=['Average_sedation_probaility_',gn,'_',num2str(tnstim),'_stim_from_stim',num2str(tfstim)];
% % fn8=['Average_sedation_probaility_',gn,'POA_',num2str(tnstim),'_stim_from_stim',num2str(tfstim)];
% saveas(gcf,[pathfig,fn8],'tiff') 
  

