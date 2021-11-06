
clear all
close all
path='I:\optogenetics\';
% pathin=[path,'outputVS1\WakeEnhance\'];
pathin=[path,'outputVSchr\']; 
pathF=[path,'Figures_WakeEnhance\']; mkdir(pathF)


% % % % % GFP
% % mousenames1=strvcat('GFP1','GFP2','GFP3','GFP4','GFP5','GFP6','GFP7','GFP8');
% % SDstim1=strvcat('160518','300418','010718','030718','041018','111018','111018','090119');  % SD+stim
% % SDonly1=strvcat('140518','020518','030718','010718','011018','091018','091018','110119');  % SD only
% 
mousenames1=strvcat('GFP1','GFP2','GFP3','GFP4','GFP6','GFP7','GFP8'); 
SDstim1=strvcat('160518','300418','010718','030718','111018','111018','090119');  % SD+stim
SDonly1=strvcat('140518','020518','030718','010718','091018','091018','110119');  % SD only

% % %%%%%%%%%% for fro
% % %%%% CHR2 mice
% % mousenames2=strvcat('GDCh1','GDCh7','GDCh8','GDCh9'); % ChR2
% % SDstim2=strvcat('040518','030718','130618','130618');  % SD+stim
% % SDonly2=strvcat('060518','010718','110618','110618');  % SD only
% % 
% % mousenames3=strvcat('GDCh5','GDCh6','GDCh16','GDCh12'); % ChR2
% % SDstim3=strvcat('040518','300418','090119','130618');  % SD+stim
% % SDonly3=strvcat('060518','020518','110119','110618');  % SD only
% % 
% % mousenames4=strvcat('GDCh18','GDCh19','GDCh20','GDCh21'); % ChR2
% % SDstim4=strvcat('080819','110819','080819','110819');  % SD+stim
% % SDonly4=strvcat('110819','080819','110819','080819');  % SD only
% % dr=1;
% % genoNs=strvcat('GFP','ChR2outVLPO','ChR2aroundVLPO','newBatch')
% % % mousenames3=strvcat('GDCh5','GDCh6','GDCh8','GDCh12','GDCh16','GDCh18','GDCh19','GDCh20','GDCh21'); % ChR2
% % % SDstim3=strvcat('040518','300418','130618','130618','090119','080819','110819','080819','110819');  % SD+stim
% % % SDonly3=strvcat('060518','020518','110618','110618','110119','110819','080819','110819','080819');  % SD only
% % % dr=1;
% % mousenames2=strvcat('GDCh1','GDCh5','GDCh6','GDCh7','GDCh8','GDCh9','GDCh12','GDCh16','GDCh18','GDCh19','GDCh20','GDCh21'); % ChR2
% % SDstim2=strvcat('040518','040518','300418','030718','130618','130618','130618','090119','080819','110819','080819','110819');  % SD+stim
% % SDonly2=strvcat('060518','060518','020518','010718','110618','110618','110618','110119','110819','080819','110819','080819');  % SD only
% % dr=1;



% %%%% CHR2 mice% 
mousenames2=strvcat('GDCh5','GDCh16','GDCh18','GDCh19','GDCh20','GDCh21'); % ChR2
SDstim2=strvcat('040518','090119','080819','110819','080819','110819');  % SD+stim
SDonly2=strvcat('060518','110119','110819','080819','110819','080819');  % SD only

mousenames3=strvcat('GDCh1','GDCh6','GDCh8','GDCh12','GDCh7','GDCh9'); % ChR2
SDstim3=strvcat('040518','300418','130618','130618','030718','130618');  % SD+stim
SDonly3=strvcat('060518','020518','110618','110618','010718','110618');  % SD only
dr=1;


genoNs=strvcat('GFP','ChR2','ChR2')


%%%%%%%%%% for occ
%%% CHR2 mice
% mousenames2=strvcat('GDCh1','GDCh5','GDCh6','GDCh7','GDCh8','GDCh9','GDCh12','GDCh16','GDCh17','GDCh18','GDCh19','GDCh20','GDCh21'); % ChR2
% SDstim2=strvcat('040518','040518','300418','030718','130618','130618','130618','090119','090119','080819','110819','080819','110819');  % SD+stim
% SDonly2=strvcat('060518','060518','020518','010718','110618','110618','110618','110119','110119','110819','080819','110819','080819');  % SD only
% dr=2;
% %%%%%%%%% for occ
% %% CHR2 mice
% mousenames2=strvcat('GDCh5','GDCh16','GDCh18','GDCh19','GDCh21'); % ChR2 %'GDCh17'
% SDstim2=strvcat('040518','090119','080819','110819','110819');  % SD+stim %'090119'
% SDonly2=strvcat('060518','110119','110819','080819','080819');  % SD only %'110119'
% % dr=2;
% 
% mousenames3=strvcat('GDCh1','GDCh6','GDCh7','GDCh8','GDCh9','GDCh12','GDCh20'); % ChR2
% SDstim3=strvcat('040518','300418','030718','130618','130618','130618','080819');  % SD+stim
% SDonly3=strvcat('060518','020518','010718','110618','110618','110618','110819');  % SD only
% dr=2;

int=1; % analyzed period 2-hour after SD


ders=['fro';'occ'];
der=ders(dr,:)


f=0:0.25:30;
maxep=10800;  %%%%21600;
zermat=zeros(1,maxep);
x=1:1:maxep; x=x./900;


SPbefore=[];SPafter=[];
SPeffect=[];

for geno=1:size(genoNs,1)%2 %1:GFP 2:ChR2 3:ChR2
    if geno==1
        mousenames=mousenames1; SDonly=SDonly1; SDstim=SDstim1; genoN=genoNs(geno,:);
    elseif geno==2
        mousenames=mousenames2; SDonly=SDonly2; SDstim=SDstim2; genoN=genoNs(geno,:);
    elseif geno==3
        mousenames=mousenames3; SDonly=SDonly3; SDstim=SDstim3; genoN=genoNs(geno,:);
    elseif geno==4
        mousenames=mousenames4; SDonly=SDonly4; SDstim=SDstim4; genoN=genoNs(geno,:);
    end


    
    numanim=size(mousenames,1);
    
    SPonly_BfAf=[];
    SPstim_BfAf=[];
    SPbSD=[];
    SPaSD=[];
%     Zscore_SDstimeffect=[];
    
    for n=1:numanim
        
        mouse=mousenames(n,:); mouse(isspace(mouse))=[];
        
        %%%%%%%%% SD+stimulation        
        day=SDonly(n,:);
        fn1=[mouse,'-',day,'-epochSO'];
        eval(['load ',pathin,fn1,'.mat epochSO -mat']);
        
        %epochSO=3600;
        
        fname=[mouse,'-',day,'-',der];
        fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
        eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
        spectr(:,1:2)=NaN;
        N=zermat; N(nr)=1;N=N(1:900*8);
        spectr=spectr(1:900*8,:);spectr(N==0,:)=NaN;
        sp1=nanmean(spectr(1:1800,:)); % 2-hour before SD
        sp2=nanmean(spectr(epochSO:epochSO+int*900,:)); % 2-hour after SD
        

        %%%%%%%%% SD+stimulation 
        day=SDstim(n,:);
        fn1=[mouse,'-',day,'-epochSO'];
        eval(['load ',pathin,fn1,'.mat epochSO -mat']);
        
        %epochSO=3600;
        
        fname=[mouse,'-',day,'-',der];
        fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
        eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
        spectr(:,1:2)=NaN;
        N=zermat; N(nr)=1;N=N(1:900*8);
        spectr=spectr(1:900*8,:);spectr(N==0,:)=NaN;
        sp3=nanmean(spectr(1:1800,:));  % 2-hour before SD
        sp4=nanmean(spectr(epochSO:epochSO+int*900,:)); % 2-hour after SD
        logsp2=log(sp2);logsp4=log(sp4);
        
        %%%%%%%%% sp1: 1st2h SDonly, sp2: after SDonly
        %%%%%%%%% sp3: 1st2h SD+stim, sp4: after SD+stim
        
        sp_difference=100+(sp4-sp2)/mean(sp2(1,3:17))*100;
        figure(10*geno) % figure(10*geno+n)
        plot(sp_difference);
        hold on
        line([0 120],[100 100],'Color','red')
        title(mouse)
        axis([0 121 40 180])
        
        
        SPonly_BfAf=[SPonly_BfAf;sp2./sp1*100];
        SPstim_BfAf=[SPstim_BfAf;sp4./sp3*100];
        
%         Zscore_SDstimeffect=[Zscore_SDstimeffect; 100*(sp4./sp3)./(sp2./sp1)];

        %         subplot(2,2,n)
        %         plot(f,sp2./sp1*100,'-b','LineWidth',2);
        %         hold on
        %         plot(f,sp4./sp3*100,'-r','LineWidth',2);
        %         axis([0 20 0 5000])
        %         plot(swa)
        %         title ([SD1(n,:),' vs ',SD2(n,:)]);
        
%         SPbSD=[SPbSD;sp3./sp1*100];
%         SPaSD=[SPaSD;sp4./sp2*100];

        
    end
    
    m1=nanmean(SPonly_BfAf);s1=nanstd(SPonly_BfAf)./sqrt(size(SPonly_BfAf,1));
    m2=nanmean(SPstim_BfAf);s2=nanstd(SPstim_BfAf)./sqrt(size(SPstim_BfAf,1));
    
    figure(100)
    subplot(1,size(genoNs,1),geno)
    errorbar(f,m1,s1,'d-k','LineWidth',1) %blue: SDonly
    hold on
    errorbar(f,m2,s2,'d-r','LineWidth',1) %bred: SD+stim
    title(genoN)
    axis([0 30 80 150])
    
    grid on
    xlabel('Frequency (Hz)')
    

%     m3=nanmean(SPbSD);s3=nanstd(SPbSD)./sqrt(size(SPbSD,1));
%     m4=nanmean(SPaSD);s4=nanstd(SPaSD)./sqrt(size(SPaSD,1));

%     Zscore_SDstimeffect; 
%     Zscore=nanmean(Zscore_SDstimeffect);
    
    SP_stimVSonly=SPstim_BfAf./SPonly_BfAf*100;
    m5=nanmean(SP_stimVSonly); s5=nanstd(SP_stimVSonly)./sqrt(size(SP_stimVSonly,1)); %s5=nanstd(SP_stimVSonly);
    SPeffect=[SPeffect m5' s5'];
    
    %%%%%%%%% for Prism
    if geno==1; 
        Prism_Zscore_stimeffect_GFP=SP_stimVSonly'; 
        Prism_SPonly_GFP=SPonly_BfAf'; 
        Prism_SPstim_GFP=SPstim_BfAf';
    elseif geno==2; 
        Prism_Zscore_stimeffect_ChR2A=SP_stimVSonly';
        Prism_SPonly_ChRA=SPonly_BfAf'; 
        Prism_SPstim_ChRA=SPstim_BfAf';
    elseif geno==3; 
        Prism_Zscore_stimeffect_ChR2B=SP_stimVSonly';
        Prism_SPonly_ChRB=SPonly_BfAf'; 
        Prism_SPstim_ChRB=SPstim_BfAf';
    end    
      
%     SPbefore=[SPbefore m3' s3'];  %%%% column1: GFP mean, column2: GFP std,  column3: ChR2 mean, column4: ChR2 std
%     SPafter=[SPafter m4' s4'];  %%%% column1: GFP mean, column2: GFP std,  column3: ChR2 mean, column4: ChR2 std

    %%%%%% Plotting all the mice
    figure(110+geno)
    for i=1:numanim
        subplot(ceil(numanim/2),2,i)
        semilogy(f,SPonly_BfAf(i,:),'k',f,SPstim_BfAf(i,:),'r','LineWidth',2)
        ylabel(['Mouse ', num2str(mousenames(i,:))])
        xlabel('Frequency (Hz)')
        axis([0 30 80 200])
    end
    
end

% figure(1)
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_1st2h_vs_afterSD2h_GFP-8nim_ChR-12anim_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_1st2h_vs_afterSD2h_GFP-8nim_ChR-12anim_',der]],'tif')
% 
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_1st2h_vs_afterSD2h_GFP-8nim_ChR-13anim_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_1st2h_vs_afterSD2h_GFP-8nim_ChR-13anim_',der]],'tif')

figure(110)
errorbar(f,SPeffect(:,1),SPeffect(:,2),'d-g','LineWidth',2) % GFP after SD
hold on
errorbar(f,SPeffect(:,3),SPeffect(:,4),'d-m','LineWidth',2) % ChR2 after SD
errorbar(f,SPeffect(:,5),SPeffect(:,6),'d-r','LineWidth',2) % ChR2 after SD
axis([0 30 80 150])
grid on
xlabel('Frequency (Hz)')


% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_afterSD2h_per_1st2h_percentile_GFP-8nim_VS_ChR-12anim_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_afterSD2h_per_1st2h_percentile_GFP-8nim_VS_ChR-12anim_',der]],'tif')
% % 
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_afterSD2h_per_1st2h_percentile_GFP-8nim_VS_ChR-13anim_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_afterSD2h_per_1st2h_percentile_GFP-8nim_VS_ChR-13anim_',der]],'tif')

tt=100*ones(size(Prism_Zscore_stimeffect_ChR2B));  tt(1:2,:)=NaN;
tt2=100*ones(size(Prism_Zscore_stimeffect_GFP));  tt2(1:2,:)=NaN;
[ha1,pa1] = ttest(Prism_Zscore_stimeffect_ChR2B',tt');
[ha2,pa2] = ttest(Prism_Zscore_stimeffect_ChR2B',tt','Alpha',0.01);
[ha3,pa3] = ttest(Prism_Zscore_stimeffect_GFP',tt2');
[ha4,pa4] = ttest(Prism_Zscore_stimeffect_GFP',tt2','Alpha',0.01);
aPrism_ha=[ha1; ha2; ha3; ha4];
% [ha2,pa2] = ttest(Prism_Zscore_stimeffect_ChR2B','Alpha',0.01);
% [ha3,pa3] = ttest(Prism_Zscore_stimeffect_ChR2B','Alpha',0.005);
% [ha4,pa4] = ttest(Prism_Zscore_stimeffect_ChR2B','Alpha',0.001);

[h1,p1] = ttest2(Prism_Zscore_stimeffect_GFP',Prism_Zscore_stimeffect_ChR2B');
[h2,p2] = ttest2(Prism_Zscore_stimeffect_GFP',Prism_Zscore_stimeffect_ChR2B','Alpha',0.01);
[h3,p3] = ttest2(Prism_Zscore_stimeffect_GFP',Prism_Zscore_stimeffect_ChR2B','Alpha',0.005);
[h4,p4] = ttest2(Prism_Zscore_stimeffect_GFP',Prism_Zscore_stimeffect_ChR2B','Alpha',0.001);
aPrism_t1=h1;
aPrism_t2=h2;
aPrism_t3=h3;
aPrism_t4=h4;




















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



% %%%%%%%%% SWA increase Stim or SDonly
% subplot(1,2,1)
% errorbar(f,SPbefore(:,1),SPbefore(:,2),'d-c','LineWidth',2) % GFP before SD
% hold on
% errorbar(f,SPafter(:,1),SPafter(:,2),'d-k','LineWidth',2) % GFP after SD
% subplot(1,2,2)
% errorbar(f,SPbefore(:,3),SPbefore(:,4),'d-m','LineWidth',2) % ChR2 before SD
% hold on
% errorbar(f,SPafter(:,3),SPafter(:,4),'d-r','LineWidth',2) % ChR2 after SD
% axis([0 30 80 140])
% grid on
% xlabel('Frequency (Hz)')

% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_1st2h_vs_afterSD2h_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_1st2h_vs_afterSD2h_',der]],'tif')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_1st2h_vs_afterSD2h_withoutGFP6_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_1st2h_vs_afterSD2h_withoutGFP6_',der]],'tif')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_1st2h_vs_afterSD2h_withoutGFP6&GDCh9_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_1st2h_vs_afterSD2h_withoutGFP6&GDCh9_',der]],'tif')

% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_leftBefore_rightAfter_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_leftBefore_rightAfter_',der]],'tif')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_leftBefore_rightAfter_withoutGFP6&GDCh9_',der]],'fig')
% saveas(gcf,[pathF,['OptoSTIMeffectsSWA_leftBefore_rightAfter_withoutGFP6&GDCh9_',der]],'tif')
