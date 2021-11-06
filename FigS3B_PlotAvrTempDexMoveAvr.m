
clear all
% close all

path='I:\optogenetics\'; % the working directory with the data folders
% path='D:\VVlab\Opto\SDstim\'; % the working directory with the data folders



geno=2; %states which animal group (genotype) I want to look at 




%Data from 24h stimulation of 2min every 20min with 10Hz blue light
% mousenames1=[6 7 8]; % names of GFP ctrl mice indicated in the file name 
% days1=['221018';'221018';'210119'];%GFP8 2min:'150119'
% firststims1=[1110;1110;640]; %810;
% stimends1=[4570;4580;1600]; %%%%%GFP:[8190;8190;7741]; %3:1710
% injections1=[980;1010;500];
% 
% % % LPO with temp
% mousenames2=[16 17 18 19 20 21]; 
% days2=['180119';'150119';'120919';'120919';'220819';'220819']; 
% %%%%%%%% 5 min,   2 min,   5 min,   5 min,   2 min,   2min
% firststims2=[995;740;1010;1010;1125;1125]; % roughly
% stimends2=[1600;2180;4115;4115;4530;4530]; %%%% 16:1940
% injections2=[780;521;750;755;991;1002]; %confirmed



mousenames1=[6 7]; % names of GFP ctrl mice indicated in the file name 
days1=['221018';'221018';'210119'];%GFP8 2min:'150119'
firststims1=[1110;1110]; %810;
stimends1=[4570;4580]; %%%%%GFP:[8190;8190;7741]; %3:1710
injections1=[980;1010];

% % LPO with temp
mousenames2=[18 19 20 21]; 
days2=['120919';'120919';'220819';'220819']; 
%%%%%%%% 5 min,   2 min,   5 min,   5 min,   2 min,   2min
firststims2=[1010;1010;1125;1125]; % roughly
stimends2=[4115;4115;4530;4530]; %%%% 16:1940
injections2=[750;755;991;1002]; %confirmed

basetemp=22; 

% % non LPO
% % mousenames2=[1 6 7 9 12]; 
% % days2=['010518';'050518';'150518';'220618';'220618']; 
% % injtime=[];
% % firststims2=[600;755;537;810;810];
% % stimends2=[7769;7854;8515;8178;8178]; %%%%if >3430+firststim, >16times stim


maxep=21600; %the total no of 4s epochs in 24h period 
epochs=1:maxep; %creates an array with 1 row and columns from 1 to 21600(maxep)
epochs=epochs./900; %edits x such that each value is divided by 900 to convert the values to hours 
            %15 epochs per 1min and 60min per 1h so 900 epochs per 1h
            %%% ./ = right-array division by dividing each element of A by the corresponding element of B
fs=256; 
epochl=4; %states how long the epoch is
zermat=zeros(1,maxep); %creates an array of zeros with 1 row and as long as the total number of epoch


pathvs1=[path,'outputVSgfp\']; %defines the folder with output files of vigilance states for mice1 (GFP ctrl)
pathvs2=[path,'outputVSchr\']; %defines the folder with output files of vigilance states for mice2 (ChR2 exp)
pathstim=[path,'STIMs\']; %gives the folder with stimulation inf
pathTemp = [path,'temp\';];
pathout=[path,'Dex\']; mkdir(pathout)
pathfig=[path,'Figures_Dex\']; mkdir(pathout)

if geno==1 %geno=1 looks at the ctrl group
    mousenames=mousenames1;days=days1; pathvs=pathvs1; gn='GFP';  injections=injections1; firststims=firststims1; stimends=stimends1;
else %geno=2 looks at the exp group
    mousenames=mousenames2; days=days2; pathvs=pathvs2; gn='GDCh';  injections=injections2; firststims=firststims2; stimends=stimends2;
end
      

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

numanim=length(mousenames); %how many animals were recorded

WallMice=[];NallMice=[];RallMice=[];TallMice=[];MallMice=[];DallMice=[];
%wake/NREM/REM/Sleep/Dex state of all mice 
      
Showtime=[0.05 4.25];
axistemp=[0.05 4.25 23 35];

Ts=[];Temps=[];
for anim=1:numanim %go through this loop as many times as there are animals
    mousename=[gn,num2str(mousenames(anim))]; mousename(isspace(mousename))=[]
    day=days(anim,:); day(isspace(day))=[]; %which days of recordings correspond to that mouse (both baseline and stim)
    firststim=firststims(anim,:); stimend=stimends(anim,:); injection=injections(anim,:);
    
    fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
    eval(['load ',pathstim,fn0,'.mat startend ']);

    stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
    stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2) 

%     outs=find((stimep2-stimep1)<25);stimep1(outs)=[];stimep2(outs)=[];

%     
% epstim=zeros(1,maxep);
% for s=1:length(stimep1); %iterating over the stims 
%     s1=stimep1(s); s2=stimep2(s); epstim(s1:s2)=1; %marks the epoch with stims (1vs0) like a VS array
% end
% ns=length(stimep1);
% v1=zeros(ns,1);v2=ones(ns,1);
% fsx=[stimep1 stimep2 stimep2 stimep1]; fsx=fsx./900;
% fsy=[v1 v1 v2 v2]; %fsy=fsy*10^5;

   
%     stimep=stimep1(tfstim:tfstim+tnstim-1); %start of stimulation episodes
%     numstim=size(stimep,1)
    



    dex_temp=zeros(1,86400);
    % %%%%%%%%%%%%%%%%%%% Temperature data aquisition

    fnt=['temp_',mousename,'_',day]; %fnt=[mousename,'_temp_calibrate_',day]; 
    eval(['load ',pathTemp,fnt,'.mat rectemp corrtemp avrcoldtemp roomtemp -mat'])
    resampled_temp=rectemp';
    corrected_temp=corrtemp';
    
%     % %%%%%%%%%%%%%%%%%%% or 
%     
%     fnt=[mousename,'_temp_calibrate_',day]; 
%     if geno==1&&anim<3
%     eval(['load ',pathTemp,fnt,'.mat resampled_temp1 corrected_temp -mat']);
%     resampled_temp=resampled_temp1';
%     corrected_temp=corrected_temp';
%     else
%     eval(['load ',pathTemp,fnt,'.mat resampled_temp corrected_temp -mat']);
%     end
%     % %%%%%%%%%%%%%%%%%%% 
    
    
    
    %%%% decide whether you use temperature data with/without correction
%     if geno==2 && anim==2
%     tanalyzed=corrected_temp;
%     else
    tanalyzed=resampled_temp;
    min(tanalyzed)
    if min(tanalyzed)<basetemp
        tanalyzed=tanalyzed+basetemp-min(resampled_temp);
    end
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if length(tanalyzed)<86400
    dex_temp(1,1:length(tanalyzed))=tanalyzed;    
    else
    dex_temp=tanalyzed(1,1:86400);  
    end

    dex_temp=reshape(dex_temp,4,[]);

    %%%% take maximum value in 4 sec, rather than 4 sec (30sec) average: by Vincent
    d2=length(dex_temp);
    dex_temp4s=zeros(1,d2);
    for tt=1:d2
    dex_temp4s(1,tt)=max(dex_temp(:,tt)); 
    end
    for tt1=1:injection
        if dex_temp4s(1,tt1)<28
            dex_temp4s(1,tt1)=dex_temp4s(1,tt1-1); 
        end
    end
    
    for tt2=injection+680:stimend %%%%% remove LOR check noise
        if dex_temp4s(1,tt2)>30
           dex_temp4s(1,tt2)=dex_temp4s(1,tt2-1); 
        end
    end

    %dex_temp4s=mean(dex_temp,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    TempD=dex_temp4s(1:6400);
    
    %%%%%%%%% moving average %%%%%%%%%%
    data=TempD;
    lag = 450;
    Tempmoveavr = movmean(data,lag);
    
    figure
    plot(TempD) %
    figure
    plot(Tempmoveavr) %
    % 
%     Tbs=TempD(1,injection-450:injection-150);
%     Td1=TempD(1,injection+1:injection+450);
%     Td2=TempD(1,injection+451:injection+900);
    
    t3h=NaN(1,3150);
                
    if stimend-injection>3000 
        t3h=Tempmoveavr(injection-449:injection+2700);
%         Td4=TempD(1,injection+901:stimend);
    else
        if injection<450 
            t3h(1,450-injection:3150)=Tempmoveavr(1,1:3150+injection-449); 
            t3h(1,stimend-injection+449:end)=NaN;
        else
            t3h=Tempmoveavr(injection-449:injection+2700); 
            t3h(1,stimend-injection+449:end)=NaN;
%         Td4=TempD(1,1801:stimend);
        end
    end
%     
    Temps=[Temps; t3h];
%     figure
%     plot(Tbs)
%     hold on
%     plot(Td1)
%     plot(Td2)
%     plot(Td4)
%     
%     Ts=[Ts; nanmean(Tbs) nanmean(Td1) nanmean(Td2) min(Td4)];    
end

figure
plot(Temps')
aPrism=Temps';

x=1/15-30:1/15:length(Temps)/15-30;
aPrismx=x';

err=nanstd(Temps)./sqrt(size(Temps,1));
tmean=nanmean(Temps);

options.color_area = [243 169 114]./255;
options.color_line = [236 112  22]./255;

x2 = [x, fliplr(x)]; %creates a vector of x and then flipped x again, required for the error bars
inBetween = [tmean+err, fliplr(tmean-err)]; %error bar vectors 
figure
patch(x2, inBetween, options.color_area,'EdgeColor', 'none', 'FaceAlpha', 0.5); %fills the area between error bars
hold on 
plot(x,tmean,'Color', options.color_line, 'LineWidth', 2) %plots the single average for all mice 
xticks(-30:15:180)
xlim([-30 180])
xlabel('Time from injection[min]')
ylabel('Temperature[C]')
% line([0 0], [0 1], 'Color', 'm', 'LineWidth', 1.5)
% title(['Average sedation probaility ',gn,' ',num2str(tnstim),'stimulations from the ',num2str(tfstim),'st stimulus'])

Ts
        