
clear all
close all

%%%%%%%% run PlotEEGandEMGprofilesDaysWriteOutVar to make fnout1=[mousename,'-EEGfrontal-EMGv-',day];
geno=2;

path='I:\optogenetics\'; % the working directory with the data folders

%%%%%%% 12h
Showtime=[0.01 12.01];
axistemp=[0.01 12.01 20 35];
%%%%%%%%%

gn='GDCh';
der=1;
pathvs=[path,'outputVSchr\'];

day='080919';  %ZT0: 080919, ZT12: 160919
% 
mousename='GDCh18' %
axisEMG=[Showtime 0 2*10^4];
axisSWA=[Showtime 0 5*10^3];

% mousename='GDCh19' %
% axisEMG=[Showtime 0 1*10^4];
% axisSWA=[Showtime 0 5*10^3];

% mousename='GDCh20' %
% axisEMG=[Showtime 0 4*10^4];
% axisSWA=[Showtime 0 4*10^3];

% mousename='GDCh21' %
% axisEMG=[Showtime 0 5*10^4];
% axisSWA=[Showtime 0 5*10^3];


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


pathstim=[path,'STIMs\']; %gives the folder with stimulation inf
pathout=[path,'1h10Hz\']; mkdir(pathout)
pathTemp = [path,'temp\';];
pathfig=[path,'1h10Hz\']; mkdir(pathout)
pathV=['I:\OptoMod\OutputSIGvar\']; %archive\

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
eval(['load ',pathstim,fn0,'.mat startend ']);

stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2) 

outs=find((stimep2-stimep1)<800);stimep1(outs)=[];stimep2(outs)=[]; % 900: 1h

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
    
    
if size(nr,1)==1; nr=nr'; w=w'; r=r'; w1=w1'; nr2=nr2'; r3=r3'; end
W=zermat; W([w;w1;mt])=1; N=zermat; N([nr;nr2])=1; R=zermat; R([r;r3])=1;
art=[w1;nr2;r;ma];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot noiseless Delta
cleanSp=NaN(maxep,121);
cleanSp(1:length(spectr),:)=spectr;
cleanSp(art)=NaN;
clear w nr r w1 nr2 r3 mt ma bastend ; %clears the variables for the next loop 
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
    
clear w nr r w1 nr2 r3 mt ma bastend ; %clears the variables for the next loop 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot noiseless EMG
fnemg=[mousename,'-EMGv-optstim-',day];
eval(['load ',pathV,fnemg,'.mat EMGv -mat']);
cleanEMG=EMGv(1:maxep); %EMG=EMGv(1:maxep);
cleanEMG(art)=NaN;


% %%%%%%%%%%%%%%%%%%% Temperature data aquisition

% attention: sampling rate is 100ms

fnt=['temp_',mousename,'_',day]; 
eval(['load ',pathTemp,fnt,'.mat rectemp roomtemp avrcoldtemp corrtemp -mat']);
temp4sec=reshape(rectemp,[],21600); 
Temp=nanmean(temp4sec,1);

%%%% baseline correction
% if avrcoldtemp<20
%     tmp1=tmp1+(20-avrcoldtemp);
% end
%%%%%%%%%%%%%%%%

% tmp4s=reshape(temp,4,[]);
% 
% %%%% take maximum value in 4 sec, rather than 4 sec (30sec) average: by Vincent
% d2=length(tmp4s);
% tmp4s=zeros(1,d2);
% for tt=1:d2
% tmp4s(1,tt)=max(tmp1(:,tt)); 
% end
% %dex_temp4s=mean(dex_temp,1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TempD=tmp4s(1:maxep);

%%%%%%%%% moving average %%%%%%%%%%
data=Temp;
lag = 450;
Tempmoveavr = movmean(data,lag);
Tempmoveavr=Tempmoveavr(1:maxep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 
x1=1:maxep; x1=x1./900;
options.color_line1 = [217 40 180]./255; % SWA
options.color_line2 = [80 80 80]./255; % EMG
options.color_line3 = [243 169 100]./255; % Temp
options.color_line4 = [255 0 0]./255; % Temp2
%[40 200 64]./255; % green
%[243 169 114]./255; % orange
% [40 40 40]./255; % grayish black

    co = [0 114/255 178/255;    
      0 158/255 115/255; %230/255 159/255 0;
      230/255 159/255 0; %213/255 94/255 0;
    204/255 121/255 167/255];
    set(groot,'defaultAxesColorOrder',co)
% %%%%% color coded SWA
%     co = [0  40 220    
%          80 200   0
%         240  20   0
%         217  50 160 %126  47 142
%         119 172  48
%         76  190 238
%         162  20  47]./255;
%     set(groot,'defaultAxesColorOrder',co)
%     

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
for vs=1:3
    if vs==1 v=W; cl=co(1,:); elseif vs==2; v=N; cl=co(2,:); elseif vs==3; v=R; cl=co(3,:); end % options.color_line1
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
figpatch.EdgeColor='none'; figpatch.FaceColor=[0 0.5 1]; %figpatch.LineWidth=2;
figpatch.FaceAlpha=0.1;
end
plot(x1,Tempmoveavr,'LineWidth',2,'Color',options.color_line4) %TrmpD
hold on
axis(axistemp)
set(gca,'YTick',[20:5:35],'TickDir','out','FontSize',14)% set(gca, 'YTick', []); 
set(gca, 'XTick', []); % set(gca,'XTick',[0:0.5:12],'TickDir','out')
ylabel('Temperature[°C]','FontSize',16)

% EMG
subplot ('position',[0.1 0.48-0.086 0.8 0.08])
for ii=1:ns
figpatch=patch(fsx(ii,:),fsy(ii,:)*10^5,'c');
hold on
figpatch.EdgeColor='none'; figpatch.FaceColor=[0 0.5 1]; %figpatch.LineWidth=2;
figpatch.FaceAlpha=0.1;
end
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
figpatch.EdgeColor='none'; figpatch.FaceColor=[0 0.5 1]; %figpatch.LineWidth=2;
figpatch.FaceAlpha=0.1;
end
plot(x1,swa) % SWA %plot(x1,Delta,'LineWidth',0.5,'Color',options.color_line1) 
% if geno==1; 
    axis(axisSWA);  
% else;
%     axis([Showtime 0 ceil(max(Delta)/1000)*1000]);   
% end
set(gca, 'YTick', []); set(gca, 'XTick', []); % set(gca,'XTick',[0:0.5:12],'TickDir','out')
ylabel('SWA [uV^2/0.25Hz]','FontSize',16)% ylabel('Delta/SWA [uV^2/0.25Hz]')

set(gca, 'XTick', [0:1:max(x1)],'FontSize',16); xlabel('ZT [hour]','FontSize',16); 
%%%% for paper figure
% set(gca, 'visible', 'off');

%%%% for paper figure
% set(gca, 'visible', 'off');



%%%%%%%%%%%%%%%%%%%%%%%%%

figname=[mousename,'_',day,'_pectrogram_Hypno_EMG_MoveAvrTemp_1h10Hz_L_Temp'];
saveas(gcf,[pathfig,figname],'tif')
% saveas(gcf,[pathfig,figname],'svg')
