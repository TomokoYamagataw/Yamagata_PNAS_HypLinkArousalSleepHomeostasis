
clear all
close all

path='I:\optogenetics\'; % the working directory with the data folders

%Data from 24h stimulation of 2min every 20min with 10Hz blue light
mousenames1=[1 2 3 4 5 6 7 8]; % names of GFP ctrl mice indicated in the file name 
days1=['120518';'050518';'180718';'180718';'111018';'221018';'221018';'150119'];
% firststims1=[888;755;753;753;809;1110;1110;810;];
% stimends1=[8230;7854;7952;7952;7954;8190;8190;7741;]; %%%%if >3430+firststim, >16times stim
% 
% % LPO
mousenames2=[5 16 17 18 19 20 21]; 
days2=['010518';'150119';'150119';'220819';'220819';'220819';'220819']; 
% firststims2=[600;758;758;1125;1125;1125;1125];
% stimends2=[7769;4064;4064;4530;4530;4530;4530]; %%%%if >3430+firststim, >16times stim

% non LPO
mousenames4=[1 6 7 8 9 12]; 
days4=['010518';'050518';'150518';'220618';'220618';'220618']; 
% injtime=[];
% firststims2=[600;755;537;810;810];
% stimends2=[7769;7854;8515;8178;8178]; %%%%if >3430+firststim, >16times stim

geno=4; %states which animal group (genotype) I want to look at 

tfstim=9;
tnstim=4; %%%% stim per 225, jitter 22.5, stimlength=30, - 4times: 675+jitter=730, 5times: 900+jitter=955; 6times: 1125+jitter=1180, 8times: 1800+jitter=1630, 16times: 1800+jitter=3430

ne=15; %there are 15 epochs in 1min 
before=2;  %we want to analyze the 2min before stimulus onset and 2min after
after=4;   %gives me 2min before, 2min of stim, and the 2min after 



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




pathvs1=[path,'outputVSgfp\']; %defines the folder with output files of vigilance states for mice1 (GFP ctrl)
pathvs2=[path,'outputVSchr\']; %defines the folder with output files of vigilance states for mice2 (ChR2 exp)
pathstim=[path,'STIMs\']; %gives the folder with stimulation inf
pathout=[path,'Dex\']; mkdir(pathout)
pathfig=[path,'Figures_Dex\']; mkdir(pathout)

if geno==1 %geno=1 looks at the ctrl group
    mousenames=mousenames1;days=days1; pathvs=pathvs1; gn='GFP'; %firststims=firststims1;
elseif geno==2
    mousenames=mousenames2; days=days2; pathvs=pathvs2; gn='GDCh';  %firststims=firststims2; 
else
    mousenames=mousenames4; days=days4; pathvs=pathvs2; gn='GDCh';  %firststims=firststims2; 
end
      

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

numanim=length(mousenames); %how many animals were recorded

WallMice=[];NallMice=[];RallMice=[];TallMice=[];MallMice=[];DallMice=[];
%wake/NREM/REM/Sleep/Dex state of all mice 
      
% % this loop goes though all the animals
% % for each 2min before and after we look at amount of sleep
% % Question: how efficient is stimulation in waking up the mice
 Ls=NaN(numanim*tnstim, 3);Lmeans=[]; Wn=[];
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

    outs=find((stimep2-stimep1)<28);stimep1(outs)=[];stimep2(outs)=[];
    
    fn=[mousename,'-',day,'-',deri]; %makes the full file name for the vigilance state and derivation output desired
    
    eval(['load ',pathvs,fn,'.mat spectr w nr r s w1 nr2 r3 s4 mt ma bastend -mat']);%
    
    VS=zeros(9,maxep); 
    VS(1,w)=1; VS(2,w1)=1; VS(3,nr)=1; VS(4,nr2)=1; VS(5,r)=1; VS(6,r3)=1; VS(7,mt)=1; VS(8,s)=1; VS(9,s4)=1;
    clear w nr r s w1 nr2 r3 s4 mt ma bastend ; %clears the variables for the next loop 
    
    wake=sum(VS(1:2,:)); 
    move=VS(7,:);
    dex=sum(VS(8:9,:));
    
        
    stimep=stimep1(tfstim:tfstim+tnstim-1); %start of stimulation episodes
    
    numstim=size(stimep,1) %how many stimulation episodes there are 
    
    Ws=[]; Ds=[]; Ms=[]; LmatD=[]; Wamounts=[];%creates empty arrays to be filled 
    for s=1:numstim %iterating over all the stimulations
        step=stimep(s); 
        eps=step-before*ne:step+after*ne; %returns the desired period
        eps1=step-before*ne:step-1; %returns the desired period
        eps2=step+1:step+2*ne; %returns the desired period
        
        ms=move(eps);  %total sleep (nrem and rem)
%         Ms=[Ms;ms]; 
        ws=wake(eps);  %for this stimulation period only
        aw=ws+ms;
        Ws=[Ws;ws+ms];  %for each of the 72 stims 
        
        ds=dex(eps);  %total dex
        Ds=[Ds;ds];
            
        if ds(1,30)==1
            LDmat=ds(1,31:60);
            LD1=min(find(LDmat<1))*4;
            Dmat=nansum(LDmat,2);
            Wmat=30-Dmat;
        
            if LD1<120; 
            else; LD1=120;
            end
        else
            LD1=NaN;
        end
        
        Ls(numstim*(anim-1)+s,:)=[mousenames(anim) s LD1];

        Wamounts=[Wamounts Wmat];
        LmatD=[LmatD LD1];

    end
% % 
% %     %%%% when program run for 16 animals
% %     fn4=[mousename,'_',day,'_VS_sedation_16times_2minpremidpost'];
% %     eval(['save ',pathout,fn4,'.mat tfs Ws Ms Ds -mat']);
   
   %mean(VS) will average the observation period across all stimulations
   %allMice arrays combine the mean data from all the mice 
   WallMice=[WallMice;mean(Ws)]; %adds the Ws of this mouse to the group array 
%    MallMice=[MallMice;mean(Ms)];
   DallMice=[DallMice;mean(Ds)];    
   
           
   LmeanD=nanmean(LmatD);
   Lmeans=[Lmeans; mousenames(anim) LmeanD];
   mWn=nanmean(Wamounts);
   Wn=[Wn; mousenames(anim) mWn];
end

aPrism_Wake=WallMice';
aPrism_Dex=DallMice';

aPrism_Latency=Lmeans';
aPrism_wake_amount_2min=Wn';
LatencyDex=[mousenames; DallMice(:,31:60)']


%plot(x, mean(WallMice)) %plots the average curve for all mice in total 
                            %(no error bars)

%errorbar(x,mean(WallMice),err) %plots y versus x and draws a vertical error bar at each data point.
                                %--> very messy
                                
data=DallMice; %states what I want to plot
data_mean = mean(data,1); %1 specifies the dimensions, in this case across rows
data_std  = std(data,0,1); %0 is weight (default=0) and 1=dimensions (across rows)
err = (data_std./sqrt(size(data,1))); %SEM error 

options.color_area = [243 169 114]./255;
options.color_line = [236 112  22]./255;

x2 = [x, fliplr(x)]; %creates a vector of x and then flipped x again, required for the error bars
inBetween = [data_mean+err, fliplr(data_mean-err)]; %error bar vectors 
figure
patch(x2, inBetween, options.color_area,'EdgeColor', 'none', 'FaceAlpha', 0.5); %fills the area between error bars
hold on 
plot(x,data_mean,'Color', options.color_line, 'LineWidth', 2) %plots the single average for all mice 
xticks(-120:20:120)
yticks(0:0.2:1)
xlim([-120 120])
ylim([0 1])
xlabel('Time (s)')
ylabel('Sedation Probability')
line([0 0], [0 1], 'Color', 'b', 'LineWidth', 1.5)
title(['Average sedation probaility ',gn,' ',num2str(tnstim),'stimulations from the ',num2str(tfstim),'st stimulus'])

% fn8=['Average_sedation_probaility_',gn,'_',num2str(tnstim),'_stim_from_stim',num2str(tfstim)];
% saveas(gcf,[pathfig,fn8],'tiff') 
  






% options.color_area = [128 193 219]./255;    % Blue theme
% options.color_line = [ 52 148 186]./255;
% options.color_area = [243 169 114]./255;    % Orange theme
% options.color_line = [236 112  22]./255;
