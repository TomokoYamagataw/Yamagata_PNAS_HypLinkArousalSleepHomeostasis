
clear all
close all
path='I:\Optogenetics\'


%%%% 1st day: baseline, 2nd day: stimulation day
mousenames1=[1 2 3 4 5 6 7 8]; % names of GFP ctrl mice indicated in the file name 
days1=['250218 260218';'290418 120418';'300618 260618';'300618 040718';'051018 061018';'081018 061018';'081018 061018';'080119 070119'];


% mousenames1=[1 2 3 4 6 7 8]; % GFP5:vis, 051018(7:50 clash)
% days1=['250218 260218';'290418 120418';'300618 260618';'300618 040718';'081018 061018';'081018 061018';'080119 070119'];
% mousenames2=[5 6 15 16 17 19 21]; 
% days2=['020418 030418';'290418 120418';'181218 191218';'181218 191218';'181218 191218';'100819 160819';'100819 160819'];
% mousenames2=[5 6 8 12 15 16 17 18 19 20 21]; 
% days2=['020418 030418';'290418 120418';'130518 140518';'100618 090618';...
%     '181218 191218';'181218 191218';'181218 191218';...
%     '100819 160819';'100819 160819';'070819 180919';'100819 160819']; 
% mousenames3=[5 15 16 17 18 19 21];
% days3=['020418 030418';'181218 191218';'181218 191218';'181218 191218';...
%     '100819 160819';'100819 160819';'100819 160819']; 
% mousenames4=[1 6 7 8 9 12 20]; %ChR2 mice,  %GDCh8: 050618 to 150618
% days4=['130218 140218';'290418 120418';'120518 130518';'130518 140518';...
%     '100618 150618';'100618 090618';'070819 180919']; %GDCh20:160819

mousenames4=[1 6 7 8 9 12]; %ChR2 mice, n
days4=['130218 140218';'290418 120418';'120518 130518';'130518 140518';...
    '100618 150618';'100618 090618']; 

mousenames5=[5 15 16 17 18 19 20 21]; 
days5=['020418 030418';'181218 191218';'181218 191218';'181218 191218';...
    '100819 160819';'100819 160819';'070819 180919';'100819 160819']; 

%%%%% for
% mousenames8=[16 17 18 19 20 21];
% days8=['181218 060119';'181218 060119';'100819 140819';'100819 140819';'070819 140819';'100819 140819']; 
mousenames8=[16 18 19 20 21];
days8=['181218 060119';'100819 140819';'100819 140819';'070819 140819';'100819 140819']; 




geno=5;





ders=strvcat('fro','occ'); 
der=1; 
deri=ders(der,:);

maxep=21600;
x=1:maxep; x=x./900;
epochl=4;
zermat=zeros(1,maxep);
fs=256;

f=0:0.25:30; %frequency is recorded in quarter Hz resolution from 0-30Hz, --> 121 bins

pathvs1=[path,'outputVSgfp\'];
pathvs2=[path,'outputVSchr\'];
pathstim=[path,'STIMs\'];

vsname=strvcat('Wake','NREM','REM');

WS=[]; ES=[]; SWE=[]; SWAamp=[]; SWAfreq=[]; SWANREM=[];

if geno==1 mousenames=mousenames1; days=days1; pathvs=pathvs1; gn='GFP'; 
elseif geno==2; mousenames=mousenames2; days=days2; pathvs=pathvs2; gn='GDCh'; 
elseif geno==3; mousenames=mousenames3; days=days3; pathvs=pathvs2; gn='GDCh'; 
elseif geno==4; mousenames=mousenames4; days=days4; pathvs=pathvs2; gn='GDCh';
elseif geno==5; mousenames=mousenames5; days=days5; pathvs=pathvs2; gn='GDCh';
elseif geno==8; mousenames=mousenames8; days=days8; pathvs=pathvs2; gn='GDCh';
end

numanim=length(mousenames);

for dd=1:2 %iterating over the days 
    
    spec=[]; specL=[]; specD=[]; %spec is day specific, animal total
    
    for anim=1:numanim %iterating over the animals
        mousename=[gn,num2str(mousenames(anim))]; %which mouse?
        daysi=days(anim,:); is=find(isspace(daysi)); %which dates for that mouse?
       
        
        day=daysi(is(1):end); day(isspace(day))=[]; %stimulation day
        fn2=[mousename,'_',day,'_stim'];eval(['load ',pathstim,fn2,'.mat startend -mat']); %importing stimulation marks
        stimep1=round(startend(:,1)./(fs*epochl)); % stim start (in epoch indices)
        stimep2=round(startend(:,2)./(fs*epochl)); % stim end (in epoch indices)
        epstim=zeros(1,maxep);
        
        for s=1:length(stimep1); %iterating over the stims 
            s1=stimep1(s); s2=stimep2(s); epstim(s1:s2)=1; %marks the epoch with stims (1vs0) like a VS array
        end
        
        if dd==1 day=daysi(1:is(1)-1); else day=daysi(is(1):end); end %which recording to look at 
        day(isspace(day))=[];
        fn1=[mousename,'-',day,'-',deri];
        eval(['load ',pathvs,fn1,'.mat w nr r w1 nr2 r3 mt ma bastend spectr -mat']);
        % w/nr/r and others give an indices of epochs with that state
        
        VS=zeros(1,maxep);
        
        VS(1,nr)=1;  if geno==5&&dd==2&&anim==7; VS(1,nr2)=1; end
        spectr(VS(1,:)==0,:)=NaN;
        
        swa=nanmean(spectr(:,3:17),2);
        swa=swa';
       
        swe=nansum(reshape(swa,1800,12));
        %%%%%%%%%%%%% Correction for length of recordings %critical for
            maxVS=max([max(nr) max(w) max(w1) max(r)]); 
            if maxVS>1800*11
            swe(1,12)=swe(1,12)+swe(1,12)*(1800*12-maxVS)/(maxVS-1800*11);
            end
        %%%%%%%%%%%%%
        cumswe=cumsum(swe);

        for light=1:2
            if light==1; spL=spectr(1:10800,:); meanspL=nanmean(spL);
            elseif light==2; spD=spectr(10801:21600,:);  meanspD=nanmean(spD);
            end
        end
        
        %meansp=nanmean(spectr);
        if geno==1 && mousenames(anim)==5; meansp=NaN(1,121); meanspL=NaN(1,121); meanspD=NaN(1,121); end     

        specL=[specL;meanspL]; 
        specD=[specD;meanspD];  %spec=[spec;meansp];%adds average freq power
        
        SWE=[SWE; cumswe];
        
        swaA=spectr(:,3:17);
        swaBand=f(3:17);
        SWAa=[];
        for t=1800:1800:21600
            SWAa=[SWAa; nanmean(swaA(t-1799:t))];
        end
        SWANREM=[SWANREM SWAa];
        [peaks freq]=max(SWAa, [], 2);
        freq=swaBand(freq);
        freq=freq';
        SWAamp=[SWAamp, peaks];
        SWAfreq=[SWAfreq, freq];


    end

    m1=nanmean(specL); m2=nanmean(specD); %averages across animals 
    s1=nanstd(specL); s2=nanstd(specD); %SD error of m1
    WS=[WS;[m1;m2]]; %two days of averages across all animals 
    ES=[ES;[s1;s2]];
    
    if dd==1; spec1=specL; spec2=specD; else; spec3=specL; spec4=specD; end %spec1 is baseline, spec2 is stim 
    
end

SWE1=SWE(1:numanim,:)./SWE(1:numanim,end)*100;
SWE2=SWE(numanim+1:2*numanim,:)./SWE(1:numanim,end)*100; %devided by baseline ZT24
SWErel=[SWE1; SWE2];

steep1=SWE1(:,6); % per 2 hours, baseline day, Light
steep2=SWE2(:,6);% per 2 hours, STIM day, Light
steep3=(SWE1(:,12)-SWE1(:,6));% per 2 hours, baseline day, Dark
steep4=(SWE2(:,12)-SWE2(:,6));% per 2 hours, STIM day, Dark

aP_steepL=[steep1' steep2'];
aP_steepD=[steep3' steep4'];
%%%%%%%%%% Fig3E' sleep
aP_steepLD=[aP_steepL; aP_steepD];
%%%%%%%%%%%%%%%%%%%%
SWED2=nanmean(SWE(:,10:12),1);


aPrism_SWA=SWANREM;
aPrism_SWE=[(SWE(1:numanim,10:12)) (SWE(numanim+1:2*numanim,10:12))];
% % aPrism_SWED2=[(SWED2(1:numanim,1))' (SWED2(numanim+1:2*numanim,1))'];
aPrism_SWE_log=log(SWE);

%%%%%%%%%% Fig3E
aPrism_SWE_percentage=SWErel';
%%%%%%%%%%%%%%%%%%%%


aPrism_spectra_baseline_L=spec1';
aPrism_spectra_baseline_D=spec2';
aPrism_spectra_STIM_L=spec3';
aPrism_spectra_STIM_D=spec4';
aPrism_x=f';
y1=numanim*ones(121,1);

%%%%%%%%%% Fig3F
aPrism_spectra=[WS(1,:)' ES(1,:)' y1 WS(2,:)' ES(2,:)' y1 WS(3,:)' ES(3,:)' y1 WS(4,:)' ES(4,:)' y1];
%%%%%%%%%%%%%%%%%%%%
[r0,p0]=ttest(nanmean(spec1(:,3:17),2),nanmean(spec3(:,3:17),2))
tm=nanmean(nanmean(spec3(:,3:17),2)./nanmean(spec1(:,3:17),2))

[r1,p1]=ttest(log(spec1),log(spec3)); %compares baseline day and stim day for each freq bin

[r2,p2]=ttest(log(spec2),log(spec4));
[rL,pL]=ttest(nanmean(log(spec1(:,3:17)),2),nanmean(log(spec3(:,3:17)),2))


%%%%%%%%%% Fig3F p-value
aP_r1=r1';
%%%%%%%%%%%%%%%%%%%%


%plotting the average curvfigure(4*geno-3)
figure(2)
semilogy(f,WS,'LineWidth',2)
data_mean=WS;
data_std=ES;
err=(data_std./sqrt(numanim));
f2=[f, fliplr(f)];
inBetween1 = [data_mean(1,:)+err(1,:), fliplr(data_mean(1,:)-err(1,:))];
inBetween2 = [data_mean(3,:)+err(3,:), fliplr(data_mean(3,:)-err(3,:))];
patch(f2, inBetween1, 'm' ,'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on
patch(f2, inBetween2, 'r' ,'EdgeColor', 'none', 'FaceAlpha', 0.3);
semilogy(f,data_mean(1,:),'m','LineWidth',2)
semilogy(f,data_mean(2,:),'c','LineWidth',2)
semilogy(f,data_mean(3,:),'r','LineWidth',2)
semilogy(f,data_mean(4,:),'b','LineWidth',2)
xlabel('Frequency (Hz)')
%ylabel('Wake Probability')
title('Mean NREM power spectra')
hold off 
axis([0.5 30 0 2.5*10^3])



%same as PLOT(...), except a log scale is used for the Y-axis.

%Plotting all the mice
figure(4*geno-2)
for i=1:numanim
    subplot(ceil(numanim/2),2,i)
    semilogy(f,spec1(i,:),'k',f,spec3(i,:),'b','LineWidth',2)
    ylabel(['Mouse ', num2str(mousenames(i))])
    xlabel('Frequency (Hz)')
    axis([0.5 30 0 4*10^3])
end

figure(4*geno-1)
plot(nanmean(SWErel(1:numanim,:),1),'k')
hold on
plot(nanmean(SWErel(numanim+1:numanim*2,:),1),'b')
plot([0 12],[100 100],'-r')
axis([0 12 0 120])



figure(4*geno)
for i=1:numanim
    subplot(ceil(numanim/2),2,i)
    plot(SWErel(i,:),'k', 'LineWidth',2)
    hold on
    plot(SWErel(i+numanim,:),'b','LineWidth',2)
    ylabel(['Mouse ', num2str(mousenames(i))])
    xlabel('Time (h)')
    hold on
    plot([0 12],[100 100],'-r')
    axis([0 12 0 150])
    grid on

end

pause