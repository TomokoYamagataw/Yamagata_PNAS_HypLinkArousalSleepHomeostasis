
clear all
close all


%%%%%%%% run PlotEEGandEMGprofilesDaysWriteOutVar to make fnout1=[mousename,'-EEGfrontal-EMGv-',day];


path='I:\optogenetics\'; % the working directory with the data folders

%%%%%%% GDCh fro
mousename='GDCh21'
day='220819'; 
gn='GDCh';
der=1;
pathvs=[path,'outputVSchr\'];
fnt=[mousename,'_temp_calibrate_',day]; 
% fnt=[mousename,'_temp_calibrate_',day]; 

%%%%%%% 4h
Showtime=[0.05 4.25];
axistemp=[0.05 4.25 23 35];
%%%%%%%%%
%%%%%%%%2h
% Showtime=[0.75 2.25];
% axistemp=[0.755 2.255 22 36];
%%%%%%%%%
axisEMG=[Showtime 0 8*10^4];
geno=2;

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
% geno=1;


ne=15; %there are 15 epochs in 1min 
before=2;  %we want to analyze the 2min before stimulus onset and 2min after
after=2;   %gives me 2min before, 2min of stim, and the 2min after 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ders=strvcat('fro','occ');
deri=ders(der,:);

freqband=0:0.25:30;
maxep=10800; %
epochs=1:maxep; 
epochs=epochs./900; 
fs=256; 
epochl=4; %states how long the epoch is
zermat=zeros(1,maxep); %creates an array of zeros with 1 row and as long as the total number of epoch

x=-before*ne:1:after*ne; %creates an array of epoch indices in the 2min before to after period
x=x*epochl; %converts the epoch indices to seconds 


% pathvs=[path,'outputVSchr\']; %defines the folder with output files of vigilance states for mice2 (ChR2 exp)
pathstim=[path,'STIMs\']; %gives the folder with stimulation inf
pathout=[path,'Dex\']; %mkdir(pathout)
pathTemp = [path,'temp\';];
pathfig=[path,'Figures_Dex\']; %mkdir(pathout)
pathV=['I:\OptoMod\OutputSIGvar\archive\'];



vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

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
eval(['load ',pathvs,fn,'.mat spectr w nr r s w1 nr2 r3 s4 mt ma bastend -mat']);%
    
    
if size(nr,1)==1 nr=nr'; w=w'; r=r'; w1=w1'; nr2=nr2'; r3=r3'; end
W=zermat; W([w;w1])=1; N=zermat; N([nr;nr2])=1; R=zermat; R([r;r3])=1; D=zermat; D([s;s4])=1;
art=[w1;nr2;r3;s4;mt;ma];
% swa=zermat';
% swa(1:length(spectr),1)=mean(spectr(:,3:17),2);
% swa(art)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot noiseless Delta
cleanSp=spectr;
cleanSp(art)=NaN;
clear w nr r s w1 nr2 r3 s4 mt ma bastend ; %clears the variables for the next loop 
Delta=nanmean(cleanSp(:,3:17),2); Delta=Delta';

    swa=mean(spectr(:,3:17),2);
%     W=zermat; W(w)=1;
%     N=zermat; N(nr)=1;
%     R=zermat; R(r)=1;

    swaW=swa; swaW(W==0)=NaN;
    swaN=swa; swaN(N==0)=NaN;
    swaR=swa; swaR(R==0)=NaN;
    swaD=swa; swaD(D==0)=NaN;

    swa=[swaW swaN swaR swaD];
    swa=swa(1:maxep,:);
    
clear w nr r s w1 nr2 r3 s4 mt ma bastend ; %clears the variables for the next loop 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot noiseless EMG
fnout1=[mousename,'-EEGfrontal-EMGv-',day];
eval(['load ',pathV,fnout1,'.mat EMGv -mat']);
cleanEMG=EMGv(1:maxep); %EMG=EMGv(1:maxep);
cleanEMG(art)=NaN;


% %%%%%%%%%%%%%%%%%%% Temperature data aquisition
eval(['load ',pathTemp,fnt,'.mat resampled_temp corrected_temp -mat']);
%%%% decide whether you use temperature data with/without correction
dex_temp=resampled_temp(1,1:86400);
dex_temp=reshape(dex_temp,4,[]);

%%%% take maximum value in 4 sec, rather than 4 sec (30sec) average: by Vincent
d2=length(dex_temp);
dex_temp4s=zeros(1,d2);
for tt=1:d2
dex_temp4s(1,tt)=max(dex_temp(:,tt)); 
end
%dex_temp4s=mean(dex_temp,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TempD=dex_temp4s(1:maxep);

%%%%%%%%% moving average %%%%%%%%%%
data=TempD;
lag = 450;
Tempmoveavr = movmean(data,lag);
% dex_temp(4082:4216)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 
x1=1:maxep; x1=x1./900;
options.color_line1 = [217 40 180]./255; % SWA
options.color_line2 = [80 80 80]./255; % EMG
options.color_line3 = [243 169 100]./255; % Temp
%[40 200 64]./255; % green
%[243 169 114]./255; % orange
% [40 40 40]./255; % grayish black

%%%%% color coded SWA
    co = [0  40 220    
         80 200   0
        240  20   0
        217  50 160 %126  47 142
        119 172  48
        76  190 238
        162  20  47]./255;
    set(groot,'defaultAxesColorOrder',co)
    

figure
% Spectrogram
c1=2; c2=9; %NREM: c1=2; c2=7; %REM: c1=0; c2=6; %ALL: c1=2; c2=8;
subplot ('position',[0.1 0.64 0.8 0.24])
pcolor(x1,freqband,log(cleanSp')); 
shading flat
caxis([c1 c2])
set(gca, 'YTick', []); 
set(gca, 'XTick', []);%set(gca,'TickDir','out')
ylabel('Frequency [Hz]','FontSize',16)% xlabel('Time(sec)')
axis([Showtime 0.5 30]) %axis([-60*before 60*after 0.5 30])
colormap(jet)

% Hypnogram
for vs=1:4
    if vs==1 v=W; cl=co(1,:); elseif vs==2; v=N; cl=co(2,:); elseif vs==3; v=R; cl=co(3,:); elseif vs==4; v=D; cl=co(4,:); end % options.color_line1
    subplot('position',[0.1 0.6288-0.011*(vs-1) 0.8 0.011]);
    bar(x1,v,1,'FaceColor',cl);
    hold on
    set(gca, 'visible', 'off');
    axis([Showtime 0 1]);%axis ([0 12 0 1]);
    set(gca, 'YTick', []); set(gca, 'XTick', []);
    set(gca, 'visible', 'off');
end;

% Temperature
subplot ('position',[0.1 0.6-0.125 0.8 0.12])
for ii=1:ns
figpatch=patch(fsx(ii,:),fsy(ii,:)*10^5,'c');
hold on
figpatch.EdgeColor='none'; figpatch.FaceColor=[0 0.55 1]; %figpatch.LineWidth=2;
figpatch.FaceAlpha=0.2;
end
plot(x1,TempD,'LineWidth',1,'Color',options.color_line3)
hold on
axis(axistemp)
set(gca,'YTick',[22:2:36],'TickDir','out','FontSize',14)% set(gca, 'YTick', []); 
set(gca, 'XTick', []); % set(gca,'XTick',[0:0.5:12],'TickDir','out')
ylabel('Temperature[°C]','FontSize',16)

% EMG
subplot ('position',[0.1 0.48-0.086 0.8 0.08])
% for ii=1:ns
% figpatch=patch(fsx(ii,:),fsy(ii,:)*10^5,'c');
% hold on
% figpatch.EdgeColor='none'; figpatch.FaceColor=[0 0.55 1]; %figpatch.LineWidth=2;
% figpatch.FaceAlpha=0.2;
% end
semilogy(x1,cleanEMG,'LineWidth',1,'Color',options.color_line2)
hold on
axis(axisEMG);
set(gca, 'YTick', []); set(gca, 'XTick', []); % set(gca,'XTick',[0:0.5:12],'TickDir','out')
ylabel('EMG')
set(gca, 'visible', 'off');

% SWA
subplot ('position',[0.1 0.392-0.16 0.8 0.16])
for ii=1:ns
figpatch=patch(fsx(ii,:),fsy(ii,:)*10^5,'c');
hold on
figpatch.EdgeColor='none'; figpatch.FaceColor=[0 0.55 1]; %figpatch.LineWidth=2;
figpatch.FaceAlpha=0.2;
end
plot(x1,swa) % SWA %plot(x1,Delta,'LineWidth',0.5,'Color',options.color_line1) 
if geno==1; axis(axisSWA);  
else;axis([Showtime 0 ceil(max(Delta)/1000)*1000]);   end
set(gca, 'YTick', []); set(gca, 'XTick', []); % set(gca,'XTick',[0:0.5:12],'TickDir','out')
ylabel('SWA [uV^2/0.25Hz]','FontSize',16)% ylabel('Delta/SWA [uV^2/0.25Hz]')

set(gca, 'XTick', [0:1:max(x1)],'FontSize',16); xlabel('ZT [hour]','FontSize',16); 
%%%% for paper figure
% set(gca, 'visible', 'off');

%%%% for paper figure
% set(gca, 'visible', 'off');

figname=[mousename,'_Dex_EMG_Hypno_RawMaxTemp_Spectrogram_Delta_3h'];

saveas(gcf,[pathfig,figname],'tiff')















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
  

