
clear all
close all

path='E:\Optoinh\';
pathin=[path,'outputVS\'];
pathF=[path,'FiguresWakeSuppression\']; mkdir(pathF)


%%%%%% GFP
mousenames1=[1 2 4 5]; 
SDonly1=strvcat('120421','120421','060521','100521');  % SD only
SDstim1=strvcat('080421','080421','100521','060521');  % SD+stim
%%%%%%%%%% all baseline data GFP
Base1=strvcat('070421','070421','050521','050521');

%%%%% Arch
mousenames2=[1 2 4 5 8 9]; % names of mice indicated in the file name 
SDonly2=strvcat('080421','120421','080421','100521','100521','060521');  % SD only
SDstim2=strvcat('120421','080421','120421','060521','060521','100521');  % SD+stim
%%%%%%%%%% all baseline data ChR2
Base2=strvcat('070421','070421','070421','050521','050521','050521'); 

% 
% mousenames3=strvcat('GDCh1','GDCh6','GDCh7','GDCh8','GDCh9','GDCh12'); % ChR2
% SDstim3=strvcat('040518','300418','030718','130618','130618','130618');  % SD+stim
% SDonly3=strvcat('060518','020518','010718','110618','110618','110618');  % SD only
% %%%%%%%%%% all baseline data ChR2
% Base3=strvcat('130218','290418','300618','100618','100618','100618');  % think to replace GDCh1 baseline
% 
% 


%%%%%%%% to compare to Caffeine inj

dr=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

genoNs=strvcat('GFP','Arch')

f=0:0.25:30;
maxep=10800;  %%%%21600;
zermat=zeros(1,maxep);
x=1:1:maxep; x=x./900;

ders=['fro';'occ'];

pathin0=[path,'outputVS\'];

der=ders(dr,:)

int=2;
SPbefore=[];
SPafter=[];
for geno=1:2 %1:GFP 2:ChR2
    if geno==1
        gn='Gf'; mousenames=mousenames1; Base=Base1; SDo=SDonly1; SDs=SDstim1; genoN=genoNs(geno,:);
    elseif geno==2
        gn='Ar'; mousenames=mousenames2; Base=Base2; SDo=SDonly2; SDs=SDstim2; genoN=genoNs(geno,:); 
    end
    
    numanim=size(mousenames,2);
    
    SPo_BfAf=[];
    SPn_BfAf=[];
    SPbSD=[];
    SPaSD=[];
    SPo_Af_base=[];
    SPn_Af_base=[];
    SPo_score=[];
    SPn_score=[];
    SPo_score2=[];
    SPn_score2=[];
    
    for n=1:numanim
        
        mouse=[gn,num2str(mousenames(n))]; mouse(isspace(mouse))=[];
 
        %%%%%%%%%%% baseline
        day=Base(n,:);
        epochSO=3600;
        fname=[mouse,'-',day,'-',der];
        
        fn=[mouse,'-',day,'-',der];
        eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
%         spectr(:,1:2)=NaN;
        N=zermat; N(nr)=1;N=N(1:900*12);
        spectr=spectr(1:900*12,:);spectr(N==0,:)=NaN;
        sp5=nanmean(spectr(1:1800,:));
        sp6=nanmean(spectr(epochSO:epochSO+int*900,:));
        sp7=nanmean(spectr);
        sp8=nanmean(spectr(epochSO-1800+1:epochSO+int*900,:)); %sp8=nanmean(spectr(1801:3600,:)); 
        
        %%%%%%%%%%% Only        
        day=SDo(n,:);
        fn1=[mouse,'-',day,'-epochSO'];
        eval(['load ',pathin,fn1,'.mat epochSO -mat']);
        
        fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
        eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
%         spectr(:,1:2)=NaN;
        N=zermat; N(nr)=1;N=N(1:900*12);
        spectr=spectr(1:900*12,:);spectr(N==0,:)=NaN;
        sp1=nanmean(spectr(1:1780,:));sp2=nanmean(spectr(epochSO:epochSO+int*900,:));
        
        %%%%%%%%%%% stim
        day=SDs(n,:);
        fn1=[mouse,'-',day,'-epochSO'];
        eval(['load ',pathin,fn1,'.mat epochSO -mat']);
        
        fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
        eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
%         spectr(:,1:2)=NaN;
        N=zermat; N(nr)=1;N=N(1:900*12);
        spectr=spectr(1:900*12,:);spectr(N==0,:)=NaN;
        sp3=nanmean(spectr(1:1780,:));sp4=nanmean(spectr(epochSO:epochSO+int*900,:));

        
        %%%%%%%%% sp1: 1st2h SDonly, sp2: after SDonly
        %%%%%%%%% sp3: 1st2h SD+stim, sp4: after SD+stim
        
        figure(2+geno)
        subplot(2,4,n)
        semilogy(f,sp7,'-k','LineWidth',2);
        axis([0 30 1*10^0 5*10^3])
        hold on
        semilogy(f,sp1,'-g','LineWidth',2);        
        semilogy(f,sp2,'-b','LineWidth',2);    
        semilogy(f,sp3,'-m','LineWidth',2);
        semilogy(f,sp4,'-r','LineWidth',2);
        axis([0 30 1*10^0 5*10^3])
        
%         sp1=log(sp1);sp2=log(sp2);sp3=log(sp3);sp4=log(sp4);
%         sp5=log(sp5);sp6=log(sp6);sp7=log(sp7);sp8=log(sp8);
        
%         SPo_BfAf=[SPo_BfAf;sp2./sp1*100];
%         SPn_BfAf=[SPn_BfAf;sp4./sp3*100];
        SPbSD=[SPbSD;sp3./sp1];
        SPaSD=[SPaSD;sp4./sp2];
        
        SPo_Af_base=[SPo_Af_base;sp2./sp8*100];
        SPn_Af_base=[SPn_Af_base;sp4./sp8*100];     
        
        SPo_score=[SPo_score;(sp2.*sp5)./(sp1.*sp8)*100]; % it is equal (sp2./sp8)./(sp1./sp5)
        SPn_score=[SPn_score;(sp4.*sp5)./(sp3.*sp8)*100]; % it is equal (sp4./sp8)./(sp3./sp5)    
        
        
        %         subplot(2,2,n)
        %         plot(f,sp2./sp1*100,'-b','LineWidth',2);
        %         hold on
        %         plot(f,sp4./sp3*100,'-r','LineWidth',2);
        %         axis([0 20 0 5000])
        %         plot(swa)
        %         title ([SD1(n,:),' vs ',SD2(n,:)]);
        
    end
    
    m0=100*nanmean(SPaSD./SPbSD); s0=nanstd(100*SPaSD./SPbSD)./sqrt(size(100*SPaSD./SPbSD,1));
    
    m5=nanmean(SPo_Af_base);s5=nanstd(SPo_Af_base)./sqrt(size(SPo_Af_base,1));
    m6=nanmean(SPn_Af_base);s6=nanstd(SPn_Af_base)./sqrt(size(SPn_Af_base,1));
    
    m7=nanmean(SPo_score);s7=nanstd(SPo_score)./sqrt(size(SPo_score,1));
    m8=nanmean(SPn_score);s8=nanstd(SPn_score)./sqrt(size(SPn_score,1));
    
    m1=nanmean(SPo_BfAf);s1=nanstd(SPo_BfAf)./sqrt(size(SPo_BfAf,1));
    m2=nanmean(SPn_BfAf);s2=nanstd(SPn_BfAf)./sqrt(size(SPn_BfAf,1));
    m3=nanmean(SPbSD);s3=nanstd(SPbSD)./sqrt(size(SPbSD,1));
    m4=nanmean(SPaSD);s4=nanstd(SPaSD)./sqrt(size(SPaSD,1));
    SPbefore=[SPbefore m3' s3'];  %%%% column1: GFP mean, column2: GFP std,  column3: ChR2 mean, column4: ChR2 std
    SPafter=[SPafter m4' s4'];  %%%% column1: GFP mean, column2: GFP std,  column3: ChR2 mean, column4: ChR2 std
%     
    figure(1)
    subplot(1,2,geno)
    errorbar(f,m5,s5,'d-b','LineWidth',1) %blue: SDonly
    hold on
    errorbar(f,m6,s6,'d-r','LineWidth',1) %bred: SD+stim
    title(genoN)
    axis([0 30 80 220])
    grid on
    xlabel('Frequency (Hz)')
    
    figure(2)
    subplot(1,2,geno)
    errorbar(f,m7,s7,'d-b','LineWidth',1) %blue: SDonly
    hold on
    errorbar(f,m8,s8,'d-r','LineWidth',1) %bred: SD+stim
    title(genoN)
    axis([0 30 80 220])
    grid on
    xlabel('Frequency (Hz)')   
    
    figure(7)
    errorbar(f,m4*100,s4*100,'d-','LineWidth',1) %blue: SDonly
    hold on
    axis([0 30 70 160])
    grid on
    xlabel('Frequency (Hz)')    
    
    figure(8)
    errorbar(f,m0,s0,'d-','LineWidth',1) %blue: SDonly
    hold on
    axis([0 30 70 160])
    grid on
    xlabel('Frequency (Hz)')    
    
    if geno==1; aPrismgfp1=SPo_score'; aPrismgfp2=SPn_score'; aPrismgfpdiff=(SPn_score-SPo_score)'; aPrismgfpbeforeafter=(100*SPaSD./SPbSD)';
    elseif geno==2; aPrismArch1=SPo_score'; aPrismArch2=SPn_score';  aPrismArchdiffA=(SPn_score-SPo_score)'; aPrismArch1beforeafter=(100*SPaSD./SPbSD)';
    end
    
%     spec1=nanmean(SPo_score(:,3:17),2);
%     spec2=nanmean(SPn_score(:,3:17),2);
%     [h05_t,p05_t]=ttest(spec1,spec2)
%     mean(spec2'./spec1')
%     [h_tlog,p_tlog]=ttest(nanmean(log(SPo_score(:,3:17)),2),nanmean(log(SPn_score(:,3:17)),2))    
%     nanmean(nanmean(log(SPo_score(:,3:17)),2)./nanmean(log(SPn_score(:,3:17)),2))

end

[h05,p05]=ttest2(aPrismgfpbeforeafter', aPrismArch1beforeafter');
[h01,p01]=ttest2(aPrismgfpbeforeafter', aPrismArch1beforeafter','alpha',0.01);
[h005,p005]=ttest2(aPrismgfpbeforeafter', aPrismArch1beforeafter','alpha',0.005);
[h001,p001]=ttest2(aPrismgfpbeforeafter', aPrismArch1beforeafter','alpha',0.001);

aP = [h05' 2*h01' 3*h005' 4*h001'];

saveas(gcf,[pathF,['OptoSTIMeffectsSWA_baseline_vs_afterSD2h_GFP-4anim_Arch-6anim_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_baseline_vs_afterSD2h_GFP-7anim_ChR-4anim_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_baseline_vs_afterSD2h_GFP-7nim_ChR-8anim_',der]],'tif')

% 
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_baseline2h_vs_afterSD2h_withoutGFP6']],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_baseline2h_vs_afterSD2h_withoutGFP6']],'tif')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% comment out
% %%%%%%%%% SWA increase Stim or SDonly
% figure
% subplot(1,2,1)
% errorbar(f,SPbefore(:,1),SPbefore(:,2),'d-g','LineWidth',2) % GFP before SD
% hold on
% errorbar(f,SPafter(:,1),SPafter(:,2),'d-b','LineWidth',2) % GFP after SD
% axis([0 30 80 140])
% grid on
% 
% subplot(1,2,2)
% 
% errorbar(f,SPbefore(:,3),SPbefore(:,4),'d-m','LineWidth',2) % ChR2 before SD
% hold on
% errorbar(f,SPafter(:,3),SPafter(:,4),'d-r','LineWidth',2) % ChR2 after SD
% axis([0 30 80 140])
% grid on
% xlabel('Frequency (Hz)')
% 
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_leftBefore_rightAfter_withoutGFP6']],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_leftBefore_rightAfter_withoutGFP6']],'tif')
% 
% 
% figure
% subplot(1,2,1)
% errorbar(f,SPbefore(:,1),SPbefore(:,2),'d-g','LineWidth',2) % GFP before SD
% hold on
% errorbar(f,SPbefore(:,3),SPbefore(:,4),'d-m','LineWidth',2) % ChR2 before SD
% axis([0 30 80 140])
% grid on
% 
% subplot(1,2,2)
% errorbar(f,SPafter(:,1),SPafter(:,2),'d-b','LineWidth',2) % GFP after SD
% hold on
% errorbar(f,SPafter(:,3),SPafter(:,4),'d-r','LineWidth',2) % ChR2 after SD
% axis([0 30 80 140])
% grid on
% xlabel('Frequency (Hz)')
