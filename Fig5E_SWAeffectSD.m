
clear all
close all
path='I:\optogenetics\';
pathin=[path,'outputVS1\WakeEnhance\'];
pathF=[path,'Figures\WakeEnhance\']; mkdir(pathF)


% % % GFP
% mousenames1=strvcat('GFP1','GFP2','GFP3','GFP4','GFP5','GFP6','GFP7','GFP8');
% SDstim1=strvcat('160518','300418','010718','030718','041018','111018','111018','090119');  % SD+stim
% SDonly1=strvcat('140518','020518','030718','010718','011018','091018','091018','110119');  % SD only

%%%%%%%%%%%%%%%%%% without GFP5
%%%%%%%%% GFP
mousenames1=strvcat('GFP1','GFP2','GFP3','GFP4','GFP6','GFP7','GFP8');
SDstim1=strvcat('160518','300418','010718','030718','111018','111018','090119');  % SD+stim
SDonly1=strvcat('140518','020518','030718','010718','091018','091018','110119');  % SD only

%%%%%%%%%%%%%%%%%% without GFP6
%%%%%%%%% GFP
% mousenames1=strvcat('GFP1','GFP2','GFP3','GFP4','GFP5','GFP7','GFP8');
% SDstim1=strvcat('160518','300418','010718','030718','041018','111018','090119');  % SD+stim
% SDonly1=strvcat('140518','020518','030718','010718','011018','091018','110119');  % SD only

%%%%%%%%% CHR2 mice if dr=1
% mousenames2=strvcat('GDCh1','GDCh5','GDCh6','GDCh7','GDCh8','GDCh9','GDCh12','GDCh16'); % ChR2
% SDstim2=strvcat('040518','040518','300418','030718','130618','130618','130618','090119');  % SD+stim
% SDonly2=strvcat('060518','060518','020518','010718','110618','110618','110618','110119');  % SD only
% dr=1
% mousenames2=strvcat('GDCh1','GDCh5','GDCh6','GDCh7','GDCh8','GDCh9','GDCh12','GDCh16','GDCh18','GDCh19','GDCh20','GDCh21'); % ChR2
% SDstim2=strvcat('040518','040518','300418','030718','130618','130618','130618','090119','080819','110819','080819','110819');  % SD+stim
% SDonly2=strvcat('060518','060518','020518','010718','110618','110618','110618','110119','110819','080819','110819','080819');  % SD only

% dr=1
% mousenames2=strvcat('GDCh5','GDCh16','GDCh18','GDCh19','GDCh21'); % ChR2
% SDstim2=strvcat('040518','090119','080819','110819','110819');  % SD+stim
% SDonly2=strvcat('060518','110119','110819','080819','080819');  % SD only

% mousenames2=strvcat('GDCh5','GDCh16'); % ChR2
% SDs2=strvcat('040518','090119');  % SD+stim
% SDo2=strvcat('060518','110119');  % SD only

% mousenames2=strvcat('GDCh18','GDCh19','GDCh21'); % ChR2
% SDstim2=strvcat('080819','110819','110819');  % SD+stim
% SDonly2=strvcat('110819','080819','080819');  % SD only

% % CHR2 mice if dr=2,occ
% % mousenames2=strvcat('GDCh1','GDCh5','GDCh6','GDCh7','GDCh8','GDCh9','GDCh12','GDCh16','GDCh17'); % ChR2
% % SDs2=strvcat('040518','040518','300418','030718','130618','130618','130618','090119','090119');  % SD+stim
% % SDo2=strvcat('060518','060518','020518','010718','110618','110618','110618','110119','110119');  % SD only
% % dr=2

% mousenames2=strvcat('GDCh1','GDCh5','GDCh6','GDCh7','GDCh8','GDCh12'); % ChR2
% SDs2=strvcat('040518','040518','300418','030718','130618','130618');  % SD+stim
% SDo2=strvcat('060518','060518','020518','010718','110618','110618');  % SD only

f=0:0.25:30;
maxep=10800; % 21600;
zermat=zeros(1,maxep);

genoNs=strvcat('GFP','ChR2')

ders=['fro';'occ'];
% dr=1
der=ders(dr,:)


int=1; %2: per 2hour, 1: per 1 hour, 0.5:per 30min

if int==1
x=1:1:10; 
x=x./(1/int)-int/2;
elseif int==2
x=1:1:11./int; 
x=x./(1/int)-int/2;   
else
x=1:1:10./int; 
x=x./(1/int)-int/2;   
end

for geno=1:2
    if geno==1
        mousenames=mousenames1; SDo=SDonly1; SDs=SDstim1; genoN=genoNs(geno,:);
    else
        mousenames=mousenames2; SDo=SDonly2; SDs=SDstim2; genoN=genoNs(geno,:);
    end

    numanim=size(mousenames,1);

    SP1=[];
    SP2=[];

    for n=1:numanim
        
%         if geno==1 
%             if n==2
%                 continue; 
%             end
%         end

        mouse=mousenames(n,:); mouse(isspace(mouse))=[];

        day=SDo(n,:);
        fn1=[mouse,'-',day,'-epochSO'];
        eval(['load ',pathin,fn1,'.mat epochSO -mat']);

        fname=[mouse,'-',day,'-',der]
        fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
        eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
        swa=nanmean(spectr(:,3:17),2);
        swa1=swa(1:1800); swa1=nanmean(reshape(swa1,900*int,[]));
        swa2=swa(epochSO:epochSO+5400-1); 
        swa2=nanmean(reshape(swa2,900*int,[]));
        sd=[]; sd(1:2/int)=NaN;
        swaSDo=[swa1 sd swa2]; 
        %swaSDo=swaSDo./nanmean(swaSDo(1:2/int))*100;
        %swaSDo=swaSDo./swaSDo(2)*100;

        day=SDs(n,:);
        fn1=[mouse,'-',day,'-epochSO'];
        eval(['load ',pathin,fn1,'.mat epochSO -mat']);

        fname=[mouse,'-',day,'-',der]
        fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
        eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
        swa=nanmean(spectr(:,3:17),2);
        swa1=swa(1:1800); 
        swa1=nanmean(reshape(swa1,900*int,[]));
        swa2=swa(epochSO:epochSO+5400-1); 
        swa2=nanmean(reshape(swa2,900*int,[]));
        sd=[]; sd(1:2/int)=NaN;
        swaSDs=[swa1 sd swa2]; 
        %swaSDs=swaSDs./nanmean(swaSDs(1:2/int))*100;
        %swaSDs=swaSDs./swaSDs(2)*100;

        SP1=[SP1;swaSDo];
        SP2=[SP2;swaSDs];

    end

    m1=mean(SP1);s1=std(SP1)./sqrt(numanim);
    m2=mean(SP2);s2=std(SP2)./sqrt(numanim);

    figure(1)
    subplot(1,2,geno)
    errorbar(x,m1,s1,'d-k','LineWidth',2)
    hold on
    errorbar(x,m2,s2,'d-r','LineWidth',2)
    axis([0 12 400 2000])
    grid on
    xlabel('Time [hours]')
    title(genoN)
    %pause
    
    beforeSD1=mean(SP1(:,1:2/int),2);
    beforeSD2=mean(SP2(:,1:2/int),2);
    percentSP1=SP1./beforeSD1*100;
    percentSP2=SP2./beforeSD2*100;
    
    m1=mean(percentSP1);s1=std(percentSP1)./sqrt(numanim);
    m2=mean(percentSP2);s2=std(percentSP2)./sqrt(numanim);
    
    figure(2)
    subplot(1,2,geno)
    errorbar(x,m1,s1,'d-k','LineWidth',2)
    hold on
    errorbar(x,m2,s2,'d-r','LineWidth',2)
    axis([0 12 50 200])
    grid on
    xlabel('Time [hours]')
    title(genoN)

end

% saveas(gcf,[pathF,['OptoSTIMeffectsSleepTime_2h_',der]],'tif')
% saveas(gcf,[pathF,['OptoSTIMeffectsSleepTime_2h_withoutGFP6_',der]],'tif')

% figure(1)
% saveas(gcf,[pathF,['OptoSTIMeffects_absSWA_per30min_GFP-8anim_ChR-12anim_',der]],'tif')
% figure(2)
% saveas(gcf,[pathF,['OptoSTIMeffects_percentBeforeSD_per30min_GFP-8anim_ChR-12anim_',der]],'tif')

% figure(1)
% saveas(gcf,[pathF,['OptoSTIMeffects_absSWA_per1h_GFP-8anim_ChR-12anim_',der]],'tif')
% figure(2)
% saveas(gcf,[pathF,['OptoSTIMeffects_percentBeforeSD_per1h_GFP-8anim_ChR-12anim_',der]],'tif')
% 
% figure(1)
% saveas(gcf,[pathF,['OptoSTIMeffects_absSWA_per2h_GFP-8anim_ChR-12anim_',der]],'tif')
% figure(2)
% saveas(gcf,[pathF,['OptoSTIMeffects_percentBeforeSD_per2h_GFP-8anim_ChR-12anim_',der]],'tif')
% 
