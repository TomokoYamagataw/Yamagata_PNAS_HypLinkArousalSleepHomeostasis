close all;
clear all;

path='E:\Optoinh\'; 
pathin = [path,'Calorimetry\'];

mousenames=['GAD3';'GAD4';'GAD5';'GAD6'];
days=['120421';'030421';'050421';'100421'];
 
numanim=length(mousenames);

figure
rHPs=[];rTemps=[];avrT=[];avrC=[];
for anim=1:numanim
    mouse=[num2str(mousenames(anim,:))]; mouse(isspace(mouse))=[];
    day=days(anim,:); day(isspace(day))=[]; %which days of recordings correspond to that mouse (both baseline and stim)    

    fname4=['calo_',mouse]; %=Dex %_pre
    eval(['load ',pathin,fname4,'.mat aPrism_all aPrism_HPmWsnip -mat']);
   
%     fname5=[mouse,'_DexTempAligned']; %=Dex %_pre
%     eval(['load ',pathin,fname5,'.mat aPrism_Temp -mat']);
    fname5=[mouse,'_DexTempAligned2']; %=Dex %_pre
    eval(['load ',pathin,fname5,'.mat aPrism_Temp2 -mat']);

    HP=aPrism_all(:,1); % =aPrism_HPmWsnip
    HPs = reshape(HP,8,[]);
    HPs = nanmean(HPs); %HPs = max(HPs); %HPs = mean(HPs);
    reHP1 = fillgaps(HPs,3,1); reHP1=reHP1';
    
    rHPs=[rHPs reHP1];
    
% %     Temps = reshape(aPrism_Temp,18,[]);
% %     Temps = nanmean(Temps);
% %     
% %     rTemps=[rTemps Temps'];
    rTemps=[rTemps aPrism_Temp2];
    
    
    
    %%%%%%%%%%%%%%%% for statistics: t1: -1h to 0h before inj, t2: 2h to 3h after Dex inj
    t1=nanmean(aPrism_Temp2(201:300,1)); t2=nanmean(aPrism_Temp2(501:600,1));
    avrT=[avrT; t1 t2];
    c1=nanmean(reHP1(201:300,1)); c2=nanmean(reHP1(501:600,1));
    avrC=[avrC; c1 c2];

    %%%%%%%%%%%%%%%%
end
% aPrism_avrT=avrT';
% aPrism_avrC=avrC';
aPrism_Calo = rHPs;
plot(aPrism_Calo)
aPrism_Temp = rTemps;
figure
plot(aPrism_Temp)
% aPrism_both = [rHPs rTemps]:


