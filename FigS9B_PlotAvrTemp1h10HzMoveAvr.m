
clear all
close all

path='I:\optogenetics\'; % the working directory with the data folders
% path='D:\VVlab\Opto\SDstim\'; % the working directory with the data folders

% % LPO with temp
mousenames=[18 19 20 21]; 
days=['090919 080919';'090919 080919';'090919 080919';'090919 080919']; 
%%%%%%%% 5 min,   2 min,   5 min,   5 min,   2 min,   2min

freqband=0:0.25:30;
maxep=10800; %the total no of 4s epochs in 24h period 
epochs=1:maxep; %creates an array with 1 row and columns from 1 to 21600(maxep)
epochs=epochs./900; 
fs=256; 
epochl=4; %states how long the epoch is
zermat=zeros(1,maxep); %creates an array of zeros with 1 row and as long as the total number of epoch



ne=15; %there are 15 epochs in 1min 
before=0;  %we want to analyze the 2min before stimulus onset and 2min after
after=240;   %gives me 2min before, 2min of stim, and the 2min after 
x=-before*ne:1:after*ne; %creates an array of epoch indices in the 2min before to after period
x=x*epochl; %converts the epoch indices to seconds 

pathvs=[path,'outputVSchr\']; 
gn='GDCh';
pathstim=[path,'STIMs\']; %gives the folder with stimulation inf
pathTemp = [path,'temp\';];
pathout=[path,'1h10Hz\'];  % mkdir(pathout)
pathfig=[path,'1h10Hz\']; % mkdir(pathout)

      

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

numanim=length(mousenames); %how many animals were recorded

WallMice=[];NallMice=[];RallMice=[];TallMice=[];MallMice=[];
      
Showtime=[0.05 6.05];
axistemp=[0.05 6.05 20 35];

Temps=[];normTs=[];minTs=[];
for anim=1:numanim %go through this loop as many times as there are animals

    mousename=[gn,num2str(mousenames(anim))]; mousename(isspace(mousename))=[]
    daysi=days(anim,:); is=find(isspace(daysi));  %which days of recordings correspond to that mouse (both baseline and stim)
    
   %%%%% mean baseline temp
    day=daysi(1:is(1)-1); day(isspace(day))=[];
    fnt=['temp_',mousename,'_',day]; 
    eval(['load ',pathTemp,fnt,'.mat rectemp roomtemp avrcoldtemp corrtemp -mat']);
    temp4sec=reshape(rectemp,[],21600); 
    Temp=nanmean(temp4sec,1);
    
%     %%%% take maximum value in 4 sec, rather than 4 sec (30sec) average: by Vincent
%     d2=length(tmp4s);
%     tmp4s=zeros(1,d2);
%     for tt=1:d2
%     tmp4s(1,tt)=max(tmp1(:,tt)); 
%     end
%     
    
    data=Temp(1:maxep);
    lag = 450;
    Tempmoveavr = movmean(data,lag);
    Tempmoveavr=Tempmoveavr(1:maxep);
    meanbaseT=nanmean(Tempmoveavr);

%%%%% stim day temp
    day=daysi(is(1):end); day(isspace(day))=[];  

    fnt=['temp_',mousename,'_',day]; 
    eval(['load ',pathTemp,fnt,'.mat rectemp roomtemp avrcoldtemp corrtemp -mat']);


    fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
    eval(['load ',pathstim,fn0,'.mat startend ']);

    stimep1=round(startend(:,1)./(fs*epochl)) % stim start (stim episode 1)
    stimep2=round(startend(:,2)./(fs*epochl)) % stim end (episode 2) 

    outs=find((stimep2-stimep1)<800);stimep1(outs)=[];stimep2(outs)=[]; % 900: 1h

    % %%%%%%%%%%%%%%%%%%% Temperature data aquisition

    % fnt=[mousename,'_temp_calibrate_',day]; 
    % eval(['load ',pathTemp,fnt,'.mat resampled_temp corrected_temp -mat']);
    fnt=['temp_',mousename,'_',day]; 
    eval(['load ',pathTemp,fnt,'.mat rectemp roomtemp avrcoldtemp corrtemp -mat']);


    %%%% decide whether you use temperature data with/without correction
    temp4sec=reshape(rectemp,[],21600); 
    Temp=nanmean(temp4sec,1);
    
%     if length(tanalyzed)<86400
%     dex_temp(1,1:length(tanalyzed))=tanalyzed;    
%     else
%     dex_temp=tanalyzed(1,1:86400);  
%     end

    %%%% decide whether you use temperature data with/without correction
    %%%% baseline correction
%     if avrcoldtemp<22
%         tmp1=tmp1+(22-avrcoldtemp);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     bin=4;
% %     tmps=reshape(tmp1,4,[]);
%     tmps=reshape(tmp1,bin,[]);

    %%%% take maximum value in 4 sec, rather than 4 sec (30sec) average: by Vincent
%     d2=length(tmps);
%     tmps=zeros(1,d2);
%     for tt=1:d2
%     tmps(1,tt)=max(tmp1(:,tt)); 
%     end


    if anim==3
        figure
        plot(Temp)
        
        tmpc=Temp;

        out0=find(tmpc(1,1:400)<28);
        tmpc(out0)=NaN;            
        
        out1=find(tmpc(1,401:700)<27)+400;
        tmpc(out1)=NaN;            

        out2=find(tmpc(1,1250:end)<31)+1250;
        tmpc(out2)=NaN;  
        
        figure
        plot(tmpc) 
        Temp=tmpc;
%         tmpsmooth=smoothdata(tmpc);
%         figure
%         plot(tmpsmooth);   
%         
%         Temp=tmpsmooth(1,1:21600);
    end
           
    %TempD=tmps(1:maxep);

    %%%%%%%%% moving average %%%%%%%%%%
    lag = 450; %225
    Tempmoveavr = movmean(Temp,lag,'omitnan');
    Tempmoveavr=Tempmoveavr(1:maxep);
    
    figure
    plot(Tempmoveavr)
%     t6h=NaN(1,6400);
    tmavs=Tempmoveavr(1,1:length(Tempmoveavr)/2); %640

    Temps=[Temps; tmavs];
    
    normT=tmavs-meanbaseT;
    normTs=[normTs; normT];

    minTstim=min(normT(1,1:900));
    nminTstim=find(normT==minTstim);
    
    minTall=min(normT);
    nminTall=find(normT==minTall);    
    
    minTs=[minTs; minTstim nminTstim minTall nminTall];
end

figure
plot(Temps')
figure
plot(normTs')

aPrism_minT=minTs';

Ts=normTs; %or Temps
aPrism=Ts';


x1=1/15:1/15:length(Ts)/15;
aPrismx=x1';

err=nanstd(Ts)./sqrt(size(Ts,1));
tmean=nanmean(Ts);

aPrism_mean=tmean';
aPrism_error=err';

options.color_area = [243 169 114]./255;
options.color_line = [236 112  22]./255;

x2 = [x1, fliplr(x1)]; %creates a vector of x and then flipped x again, required for the error bars
inBetween = [tmean+err, fliplr(tmean-err)]; %error bar vectors 
figure
patch(x2, inBetween, options.color_area,'EdgeColor', 'none', 'FaceAlpha', 0.5); %fills the area between error bars
hold on 
plot(x1,tmean,'Color', options.color_line, 'LineWidth', 2) %plots the single average for all mice 
% xticks(-30:15:180)
% xlim([-30 180])
xlabel('Time from stim onset [min]')
ylabel('Temperature[C]')
% line([0 0], [0 1], 'Color', 'm', 'LineWidth', 1.5)
% title(['Average sedation probaility ',gn,' ',num2str(tnstim),'stimulations from the ',num2str(tfstim),'st stimulus'])

        