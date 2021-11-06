
clear all
close all

path='E:\Optoinh\';
pathin=[path,'outputVS\'];
pathF=[path,'FiguresWakeSuppression\']; mkdir(pathF)


%%%%%% GFP
mousenames1=[1 2 4 5]; 
SDonly1=strvcat('120421','120421','060521','100521');  % SD only
SDstim1=strvcat('080421','080421','100521','060521');  % SD+stim

%%%%% Arc
mousenames2=[1 2 4 5 8 9]; % names of mice indicated in the file name 
SDonly2=strvcat('080421','120421','080421','100521','100521','060521');  % SD only
SDstim2=strvcat('120421','080421','120421','060521','060521','100521');  % SD+stim


f=0:0.25:30;
maxep=10800; % 21600;
zermat=zeros(1,maxep);

genoNs=strvcat('GFP','Arch')

ders=['fro';'occ'];
dr=1
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
        mousenames=mousenames1; SDo=SDonly1; SDs=SDstim1; genoN=genoNs(geno,:); gn='Gf';
    else
        mousenames=mousenames2; SDo=SDonly2; SDs=SDstim2; genoN=genoNs(geno,:); gn='Ar'; 
    end

    numanim=size(mousenames,2);

    SP1=[];
    SP2=[];

    for n=1:numanim
        
%         if geno==1 
%             if n==2
%                 continue; 
%             end
%         end


        mouse=[gn,num2str(mousenames(n))]; mouse(isspace(mouse))=[];

        day=SDo(n,:);
        fn1=[mouse,'-',day,'-epochSO'];
        eval(['load ',pathin,fn1,'.mat epochSO epochSDon epochSDoff -mat']);

        fname=[mouse,'-',day,'-',der]
        fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
        eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
               
        N=zermat; 
        N(nr)=1;
        swa=nanmean(spectr(:,3:17),2);
        swa(N==0,:)=NaN;
        swa=swa';

        swa1=NaN*ones(1,1800); 
        swa1(1:epochSDon)=swa(1:epochSDon);
        swa1=nanmean(reshape(swa1,900*int,[]));
        
        swa2=swa(epochSO:epochSO+5400-1); 
        swa2=nanmean(reshape(swa2,900*int,[]));
        sd=[]; sd(1:2/int)=NaN;
        swaSDo=[swa1 sd swa2]; 
        %swaSDo=swaSDo./nanmean(swaSDo(1:2/int))*100;
        %swaSDo=swaSDo./swaSDo(2)*100;

        day=SDs(n,:);
        fn2=[mouse,'-',day,'-epochSO'];
        eval(['load ',pathin,fn2,'.mat epochSO epochSDon epochSDoff -mat']);

        fname=[mouse,'-',day,'-',der]
        fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
        eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
        
        N=zermat; 
        N(nr)=1;
        swa=nanmean(spectr(:,3:17),2);
        swa(N==0,:)=NaN;
        swa=swa';        
%         swa=nanmean(spectr(:,3:17),2);
%         swa1=swa(1:1800); 
        swa1=NaN*ones(1,1800); 
        swa1(1:epochSDon)=swa(1:epochSDon);
        
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

    m1=nanmean(SP1);s1=nanstd(SP1)./sqrt(numanim);
    m2=nanmean(SP2);s2=nanstd(SP2)./sqrt(numanim);

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
    
    m1=nanmean(percentSP1);s1=nanstd(percentSP1)./sqrt(numanim);
    m2=nanmean(percentSP2);s2=nanstd(percentSP2)./sqrt(numanim);
    
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
% figure(1)
% saveas(gcf,[pathF,['OptoSTIMeffects_absSWA_per30min_GFP-8anim_ChR-12anim_',der]],'tif')
% figure(2)
% saveas(gcf,[pathF,['OptoSTIMeffects_percentBeforeSD_per30min_GFP-8anim_ChR-12anim_',der]],'tif')

figure(1)
saveas(gcf,[pathF,['OptoSTIMeffects_absSWA_per1h_GFP-4anim_Arch-6anim_',der]],'tif')
figure(2)
saveas(gcf,[pathF,['OptoSTIMeffects_percentBeforeSD_per1h_GFP-4anim_Arch-6anim_',der]],'tif')
% 
% figure(1)
% saveas(gcf,[pathF,['OptoSTIMeffects_absSWA_per2h_GFP-8anim_ChR-12anim_',der]],'tif')
% figure(2)
% saveas(gcf,[pathF,['OptoSTIMeffects_percentBeforeSD_per2h_GFP-8anim_ChR-12anim_',der]],'tif')
% 
