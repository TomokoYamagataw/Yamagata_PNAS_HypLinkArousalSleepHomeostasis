
clear all
close all

path='I:\optogenetics\'; % the working directory with the data folders

% Showtime=[1.917 2.083]; %2, 4, 12, 24 1min:0.0167, 4min:0.0667, 5min:0.0833 6min:0.1
Showtime=[6 7];
figt1=int2str(round(min(Showtime))); figt2=int2str(round(max(Showtime)));
der=1;

%%%%%%% GDCh fro
pathvs=[path,'outputVSchr\']; gn='GDCh';

% mousename='GDCh15' %
% day='191218'; 
% 
% mousename='GDCh18' %
% day='160819'; 
% axisEMG=[Showtime 0 7*10^5];
% axisSWA=[Showtime 0 5.5*10^3]; swaid=1;
% 
mousename='GDCh19' %
day='160819'; 
axisEMG=[Showtime 0 1*10^4];
axisSWA=[Showtime 0 4*10^3]; swaid=1;
% 
% mousename='GDCh20' %
% day='160819'; 
% axisEMG=[Showtime 0 1*10^5];
% axisSWA=[Showtime 0 2.7*10^3]; swaid=1;

% mousename='GDCh21' %
% day='200819'; 
% axisEMG=[Showtime 0 1*10^5];
% axisSWA=[Showtime 0 3.2*10^3]; swaid=1;


% %  reference
% % 1Hz
% mousenames2=[15 16 17 18 19 20 21];
% days2=['241218';'241218';'241218';'190819';...
%     '190819';'190819';'190819'];
% % 2Hz
% mousenames3=[15 16 17 18 19 20 21];
% days3=['251218';'251218';'251218';'200819';...
%     '200819';'200819';'200819'];
% % 5Hz
% mousenames4=[15 16 17 18 19 20 21];
% days4=['261218';'261218';'261218';'210819';...
%     '210819';'210819';'210819'];
% % 10Hz
% mousenames5=[15 16 17 18 19 20 21];
% days5=['191218';'191218';'191218';'160819';...
%     '160819';'160819';'160819'];
% % Sham
% mousenames6=[15 16 17 18 19 20 21]; 
% days6=['181218';'181218';'181218';...
%     '100819';'100819';'070819';'100819']; 





% fnt=[mousename,'_temp_calibrate_',day]; % axistemp=[0.755 2.255 25 35];



% 
% mousename='GFP7' %or GFP6
% day='221018'; 
% gn='GFP';
% der=1; 
% pathvs=[path,'outputVSgfp\'];
% fnt=[mousename,'_temp_',day]; 
% Showtime=[0.75 2.25];
% axistemp=[0.741 2.241 25 35];
% axisEMG=[Showtime 0 2*10^8];
% axisSWA=[Showtime 0 5*10^3];%GFP6&7
% swaid=1;



ne=15; %there are 15 epochs in 1min 
before=2;  %we want to analyze the 2min before stimulus onset and 2min after
after=2;   %gives me 2min before, 2min of stim, and the 2min after 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ders=strvcat('fro','occ');
deri=ders(der,:);

freqband=0:0.25:30;
maxep=21600; %
epochs=1:maxep; 
epochs=epochs./900; 
fs=256; 
epochl=4; %states how long the epoch is
zermat=zeros(1,maxep); %creates an array of zeros with 1 row and as long as the total number of epoch

x=-before*ne:1:after*ne; %creates an array of epoch indices in the 2min before to after period
x=x*epochl; %converts the epoch indices to seconds 


pathstim=[path,'STIMs\']; %gives the folder with stimulation inf
pathout=[path,'24h10Hz\']; %mkdir(pathout)
pathTemp = [path,'temp\';];
pathfig=[path,'Figures_Spectrogram_10Hz\']; mkdir(pathfig)
% pathV=['I:\OptoMod\OutputSIGvar\archive\'];
pathV=['I:\OptoMod\OutputSIGvar\'];

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states


%%%%%%%%%% stimulation timing
fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
eval(['load ',pathstim,fn0,'.mat startend ']);

stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2) 

outs=find((stimep2-stimep1)<25);stimep1(outs)=[];stimep2(outs)=[];

epstim=zeros(1,maxep);
for s=1:length(stimep1); %iterating over the stims 
    s1=stimep1(s); s2=stimep2(s); epstim(s1:s2)=1; %marks the epoch with stims (1vs0) like a VS array
end
ns=length(stimep1);
v1=zeros(ns,1);v2=ones(ns,1);
fsx=[stimep1 stimep2 stimep2 stimep1]; fsx=fsx./900;
fsy=[v1 v1 v2 v2]; %fsy=fsy*10^5;



fn=[mousename,'-',day,'-',deri]; %makes the full file name for the vigilance state and derivation output desired
eval(['load ',pathvs,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);%

    
if size(nr,1)==1 nr=nr'; w=w'; r=r'; w1=w1'; nr2=nr2'; r3=r3'; end
W=zermat; W([w;w1;mt])=1; N=zermat; N([nr;nr2])=1; R=zermat; R([r;r3])=1; % W([w;w1])=1; 
art=[w1;nr2;r3;ma];%art=[w1;nr2;r3;mt;ma];


% swa=zermat';
% swa(1:length(spectr),1)=mean(spectr(:,3:17),2);
% swa(art)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot noiseless Delta
cleanSp=spectr;
cleanSp(art)=NaN;
Delta=nanmean(cleanSp(:,3:17),2); Delta=Delta';

    swa=mean(spectr(:,3:17),2);
%     W=zermat; W(w)=1;
%     N=zermat; N(nr)=1;
%     R=zermat; R(r)=1;

    swaW=swa; swaW(W==0)=NaN;
    swaN=swa; swaN(N==0)=NaN;
    swaR=swa; swaR(R==0)=NaN;

    swa=[swaW swaN swaR];
    swa=swa(1:maxep,:);
    
clear w nr r s w1 nr2 r3 s4 mt ma bastend ; %clears the variables for the next loop 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot noiseless EMG
fnout1=[mousename,'-EMGv-optstim-',day];
eval(['load ',pathV,fnout1,'.mat EMGv -mat']);
cleanEMG=[0 0 EMGv(1:maxep-2)]; 
cleanEMG(art)=NaN;
%%%%%%%%%%%%%%%%%%% Temperature data aquisition
% eval(['load ',pathTemp,fnt,'.mat resampled_temp corrected_temp -mat']);
% dex_temp=corrected_temp(1,1:86400);
% dex_temp=reshape(dex_temp,4,[]);
% dex_temp4s=mean(dex_temp,1);
% TempD=dex_temp4s(1:maxep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 
x1=1:maxep; x1=x1./900;
options.color_line1 = [217 20 108]./255; % SWA
options.color_line2 = [80 80 80]./255; % EMG
options.color_line3 = [243 169 114]./255; % Temp

figure
% Spectrogram
c1=2; c2=10; %NREM: c1=2; c2=7; %REM: c1=0; c2=6; %ALL: c1=2; c2=8;
subplot ('position',[0.1 0.94-0.27 0.8 0.25])
pcolor(x1,freqband,log(cleanSp')); 
shading flat
caxis([c1 c2])
set(gca, 'YTick', []); set(gca, 'XTick', []);%set(gca,'TickDir','out')
% ylabel('Frequency [Hz]')% xlabel('Time(sec)')
% colorbar
axis ([Showtime 0.5 30]) %axis([-60*before 60*after 0.5 30])
colormap(jet)

% Hypnogram
for vs=1:3
    if vs==1 v=W; cl=[0 114/255 178/255]; elseif vs==2; v=N; cl=[0 158/255 115/255]; elseif vs==3; v=R; cl=[230/255 159/255 0]; end
    subplot('position',[0.1 0.67-0.012*(vs-1)-0.012 0.8 0.012]);
    bp=bar(x1,v,1);
    bp.FaceColor = cl;
    set(gca, 'visible', 'off');
    axis([Showtime 0 1]);%axis ([0 12 0 1]);
    set(gca, 'YTick', []); set(gca, 'XTick', []);
    set(gca, 'visible', 'off');
end;

% EMG
subplot ('position',[0.1 0.64-0.05 0.8 0.04])
hold on
for ii=1:ns
figpatch=patch(fsx(ii,:),fsy(ii,:)*1*10^4,'c');
hold on
figpatch.EdgeColor='none'; figpatch.FaceColor=[0 0.55 1]; %figpatch.LineWidth=2;
figpatch.FaceAlpha=0.0;%0.2
axis(axisEMG)
end
semilogy(x1,cleanEMG,'LineWidth',1,'Color',options.color_line2)
axis(axisEMG)
set(gca, 'YTick', []); set(gca, 'XTick', []); % set(gca,'XTick',[0:0.5:12],'TickDir','out')
set(gca, 'visible', 'off');
% ylabel('EMG')

% % SWA
% subplot ('position',[0.1 0.54-0.174 0.8 0.16])
% for ii=1:ns
% figpatch=patch(fsx(ii,:),fsy(ii,:)*10^5,'c');
% hold on
% figpatch.EdgeColor='none'; figpatch.FaceColor=[0 0.55 1]; %figpatch.LineWidth=2;
% figpatch.FaceAlpha=0.2;
% end
% swap=plot(x1,Delta,'LineWidth',0.5,'Color',options.color_line1);
% if swaid==1; axis(axisSWA);  elseif swaid==2; else axis([Showtime 0 ceil(max(Delta)/1000)*1000]);  end
% set(gca, 'YTick', []); set(gca, 'XTick', []); % set(gca,'XTick',[0:0.5:12],'TickDir','out')
% ylabel('Delta [uV^2/0.25Hz]')%ylabel('Delta')% 
% set(gca, 'XTick', [0:0.1:max(x1)]); xlabel('ZT [hour]'); 
%  
% color coded SWA
%     co = [0 .27 .8;    
%       .2 .6 0;
%       .86 .1 0;
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840];
    co = [0 114/255 178/255;    
      0 158/255 115/255; %230/255 159/255 0;
      230/255 159/255 0; %213/255 94/255 0;
    204/255 121/255 167/255
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
    set(groot,'defaultAxesColorOrder',co)
    
subplot ('position',[0.1 0.59-0.1 0.8 0.10])
for ii=1:ns
figpatch=patch(fsx(ii,:),fsy(ii,:)*10^5,'c');
hold on
figpatch.EdgeColor='none'; figpatch.FaceColor=[0 0.55 1]; %figpatch.LineWidth=2;
figpatch.FaceAlpha=0.0;%0.2
end
swap=plot(x1,swa);
if swaid==1; axis(axisSWA);  elseif swaid==2; else axis([Showtime 0 ceil(max(Delta)/1000)*1000]);  end
set(gca, 'YTick', []); 
% ylabel('SWA [uV^2/0.25Hz]')%ylabel('Delta')% 
set(gca, 'XTick', [0:0.25:max(x1)],'TickDir','out'); 
set(gca, 'XTicklabels', []);
%xlabel('ZT [hour]'); 

%     xticklabels({})
%     title ([tl,' ',fn]); 
%     grid on
%     hold on











% set(gca, 'XTick', [0:0.5:max(x1)]); xlabel('Time [min]'); 
% set(gca, 'visible', 'off');

% subplot ('position',[0.1 0.32-0.2 0.8 0.2])
% for ii=1:ns
% figpatch=patch(fsx(ii,:),fsy(ii,:),'c');
% hold on
% figpatch.EdgeColor='none'; figpatch.FaceColor=[0 0.55 1]; %figpatch.LineWidth=2;
% figpatch.FaceAlpha=0.2;
% end
% axis([Showtime 0 1])
% set(gca, 'visible', 'off');


% % % Temperature
% subplot ('position',[0.1 0.64 0.8 0.12])
% semilogy(x1,TempD,'LineWidth',1,'Color',options.color_line3)
% hold on
% axis(axistemp)
% set(gca, 'YTick', []); set(gca, 'XTick', []); % set(gca,'XTick',[0:0.5:12],'TickDir','out')
% ylabel('Temperature[°C]')
% % set(gca, 'visible', 'off');



figname=[mousename,'_Spectrogram_Hypno_EMG_Delta_Stim10Hz_ZT',figt1,'-ZT',figt2,'_v7'];
saveas(gcf,[pathfig,figname],'svg')
saveas(gcf,[pathfig,figname],'eps')









%%%%%%%%%%%%%%%%%%%%%%%%%
    

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
  

