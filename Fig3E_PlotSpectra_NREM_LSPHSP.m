
clear all
close all
path='I:\Optogenetics\'




% mousenames2=[5 15 16 17 18 19 21]; % GDCh8:vis, 18 %;'070819 160819'
% days2=['020418 030418';'181218 191218';'181218 191218';'181218 191218';...
%     '100819 160819';'100819 160819';'100819 160819']; 

mousenames2=[5 15 16 17 18 19 20 21];
days2=['300318 310318';'201218 211218';'201218 211218';'201218 211218';...
    '040919 060919';'040919 060919';'040919 060919';'040919 060919'];

geno=2;

dur=2; % 2: 2-hour after SD (ZT4-6)

ders=strvcat('fro','occ'); 
der=1; 
deri=ders(der,:);

maxep=10800;
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
end

numanim=length(mousenames);

for dd=1:2 %iterating over the days 
    
    spec=[]; %spec is day specific, animal total
    
    for anim=1:numanim %iterating over the animals
        mousename=[gn,num2str(mousenames(anim))]; %which mouse?
        daysi=days(anim,:); is=find(isspace(daysi)); %which dates for that mouse?
        
        if dd==1 day=daysi(1:is(1)-1); else day=daysi(is(1):end); end 
        day(isspace(day))=[]; %stimulation day
        
        fn2=[mousename,'_',day,'_stim'];eval(['load ',pathstim,fn2,'.mat startend -mat']); %importing stimulation marks
        stimep1=round(startend(:,1)./(fs*epochl)); % stim start (in epoch indices)
        stimep2=round(startend(:,2)./(fs*epochl)); % stim end (in epoch indices)
%         epstim=zeros(1,maxep);
%         
%         for s=1:length(stimep1); %iterating over the stims 
%             s1=stimep1(s); s2=stimep2(s); epstim(s1:s2)=1; %marks the epoch with stims (1vs0) like a VS array
%         end
        
        fn1=[mousename,'-',day,'-',deri];
        eval(['load ',pathvs,fn1,'.mat w nr r w1 nr2 r3 mt ma bastend spectr -mat']);
        
        VS=zeros(1,maxep);
        
        VS(1,nr)=1;
        
        spectr(VS==0,:)=NaN; 
        
        spectrRS=spectr(900*4+1:900*(4+dur),:);
    
        meansp=nanmean(spectrRS);
        swa=nanmean(spectrRS(:,3:17),2); %all rows from 3-17 bin = 0.5-4 Hz averaged across columns
        swa=swa';
       
        if geno==1 && mousenames(anim)==5; meansp=NaN(1,121); end     

        spec=[spec;meansp]; %adds average freq power 
        
        swaA=spectr(:,3:17);
        swaBand=f(3:17);
        SWAa=[];
        for t=1:dur
            SWAa=[SWAa; nanmean(swaA(900*(t-1)+1:900*t))];
        end
        SWANREM=[SWANREM SWAa];

    end
   

    m1=nanmean(spec); %averages across animals 
    s1=nanstd(spec); %SD error of m1
    WS=[WS;m1]; %two days of averages across all animals 
    ES=[ES;s1];
    
    if dd==1 spec1=spec; else spec2=spec; end %spec1 is baseline, spec2 is stim 
    
end


[r,p]=ttest2(log(spec1),log(spec2)); 

aPrism_SWA=SWANREM;
aPrism_spectra_LSP=spec1';
aPrism_spectra_HSP=spec2';
aPrism_x=f';

%plotting the average curvfigure(4*geno-3)
semilogy(f,WS,'LineWidth',2)
data_mean=WS;
data_std=ES;
err=(data_std./sqrt(numanim));
f2=[f, fliplr(f)];
inBetween1 = [data_mean(1,:)+err(1,:), fliplr(data_mean(1,:)-err(1,:))];
inBetween2 = [data_mean(2,:)+err(2,:), fliplr(data_mean(2,:)-err(2,:))];
patch(f2, inBetween1, 'k' ,'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on
patch(f2, inBetween2, 'b' ,'EdgeColor', 'none', 'FaceAlpha', 0.3);
semilogy(f,data_mean(1,:),'k','LineWidth',2)
semilogy(f,data_mean(2,:),'b','LineWidth',2)
xlabel('Frequency (Hz)')
%ylabel('Wake Probability')
title('Mean NREM power spectra')
hold off 
axis([0.5 30 0 4*10^3])



%same as PLOT(...), except a log scale is used for the Y-axis.

%Plotting all the mice
figure(4*geno-2)
for i=1:numanim
    subplot(ceil(numanim/2),2,i)
    semilogy(f,spec1(i,:),'k',f,spec2(i,:),'b','LineWidth',2)
    ylabel(['Mouse ', num2str(mousenames(i))])
    xlabel('Frequency (Hz)')
    axis([0.5 30 0 4*10^3])
end
% 
% figure(4*geno-1)
% plot(nanmean(SWErel(:,1:numanim),2),'k')
% hold on
% plot(nanmean(SWErel(:,numanim+1:numanim*2),2),'b')
% plot([0 12],[100 100],'-r')
% axis([0 12 0 120])


% 
% figure(4*geno)
% for i=1:numanim
%     subplot(ceil(numanim/2),2,i)
%     plot(SWErel(:,i),'k', 'LineWidth',2)
%     hold on
%     plot(SWErel(:,i+numanim),'b','LineWidth',2)
%     ylabel(['Mouse ', num2str(mousenames(i))])
%     xlabel('Time (h)')
%     hold on
%     plot([0 12],[100 100],'-r')
%     axis([0 12 0 150])
%     grid on

% end

% pause