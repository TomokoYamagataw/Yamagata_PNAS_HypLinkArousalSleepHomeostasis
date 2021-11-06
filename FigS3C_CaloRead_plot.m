close all;
clear all;

%%%% LD cycle 12/12 -7:00/19:00
%Calo_2020_01_09.CSV
% startRow = 499; 07:00
% endRow = 9504; 18:14
% injection = 2515 - 2555, 10:04 see calo when you inejcted
% antisedan injection 5400
% GAD1

% Calo_2020_01_10.CSV
% startRow = 3; 00:00
% endRow = Inf; 17:06
% injection = 8217 - 8262, 10:03 see calo when you injected
% antisedan injection 11372 - 11420
% GAD2


% 15 min: 196 - 211

%%%% Initialize variables.
path='E:\Optoinh\'; 
pathin = [path,'Calorimetry\'];
pathout = [path,'Calorimetry\'];

% mouse='GAD1'
% recorddate='090120'
% Dst = 2515;
% Dst2 = 2555;
% Dend = 5400; % Antisedan injection

% 
% mouse='GAD2'
% recorddate='100120'
% Dst = 8217;
% Dst2 = 8252;
% % heating start: 10370
% Dend = 11374; % Antisedan injection
% 
% mouse='GAD3'
% recorddate='110421'
% Dst = 2888;
% Dst2 = 2918; 
% Dend = 6384; % Antisedan injection

% % % % mouse='GAD4'
% % % % recorddate='100421'
% % % % Dst = 10314;
% % % % Dst2 = 10345; 
% % % % Dend = 13550; % Antisedan injection

% mouse='GAD4'
% recorddate='030421'
% Dst = 9658;
% Dst2 = 9705; 
% Dend = 13127; % Antisedan injection

% mouse='GAD5'
% recorddate='050421'
% Dst = 9398;
% Dst2 = 9423; 
% Dend = 12544; % Antisedan injection

mouse='GAD6'
recorddate='100421'
Dst = 5357;
Dst2 =5416; 
Dend = 8579 % Antisedan injection


% 


fname2=['calo_',mouse,'_',recorddate,'_Dex'];
eval(['load ',pathin,fname2,'.mat Calo TempC FlowL O2mL O2mLc CO2mL CO2mLc RER RERc HPmW HPmWc -mat']);
% 

% Calo.TempC = cell2mat(rawNumericColumns(:, 1));
% Calo.FlowLmin = cell2mat(rawNumericColumns(:, 2));
% Calo.O2mL = cell2mat(rawNumericColumns(:, 3));
% Calo.CO2mLmin = cell2mat(rawNumericColumns(:, 4));
% Calo.RER = cell2mat(rawNumericColumns(:, 5));
% Calo.HPmW = cell2mat(rawNumericColumns(:, 6));


% O2mLc=O2mL;
% for ii=2:length(O2mLc)-20
%     if O2mLc(ii,1)<=0
%         O2mLc(ii-1:ii+20,1)=NaN; 
%     end
%     ii=ii+15;
% end
% out1=isnan(O2mLc);
% 
% RERc=RER;
% RERc(out1)=NaN;
% stdRER=nanstd(RERc);
% medRER=nanmedian(RERc);
% out2=find(RERc>medRER+4*stdRER);
% RERc(out2)=NaN;


HPmWc=medfilt1(HPmWc,5);
CO2mLc=medfilt1(CO2mLc,5);
O2mLc=medfilt1(O2mLc,5);
RERc=medfilt1(RERc,5);

figure
plot(HPmWc)
hold on
figure
plot(CO2mLc)
hold on
plot(O2mLc)
plot(RERc)


gap = NaN(50,1);

nanmat2=NaN*ones(2350,1);

HPmW1=HPmWc(Dst-2400:Dst-1);
HPmW2=HPmWc(Dst2:Dst2+2349);
if length(HPmWc)<Dend+2349
    HPmW3=nanmat2; 
    HPmW3(1:length(HPmWc)-Dend+1,:)=HPmWc(Dend:length(HPmWc),1);
else
    HPmW3=HPmWc(Dend:Dend+2349);
end

aPrism_HPmWsnip=[HPmW1;gap;HPmW2;gap;HPmW3];

RER1=RERc(Dst-2400:Dst-1);
RER2=RERc(Dst2:Dst2+2349);
if length(HPmWc)<Dend+2349
    RER3=nanmat2; 
    RER3(1:length(HPmWc)-Dend+1,:)=RERc(Dend:length(HPmWc),1);
else
    RER3=RERc(Dend:Dend+2349);
end
aPrism_RERsnip=[RER1;gap;RER2;gap;RER3];

COt1=CO2mLc(Dst-2400:Dst-1);
COt2=CO2mLc(Dst2:Dst2+2349);
if length(HPmWc)<Dend+2349
    COt3=nanmat2; 
    COt3(1:length(HPmWc)-Dend+1,:)=CO2mLc(Dend:length(HPmWc),1);
else
    COt3=CO2mLc(Dend:Dend+2349);
end
aPrism_COtsnip=[COt1;gap;COt2;gap;COt3];

Ot1=O2mLc(Dst-2400:Dst-1);
Ot2=O2mLc(Dst2:Dst2+2349);
if length(HPmWc)<Dend+2349
    Ot3=nanmat2; 
    Ot3(1:length(HPmWc)-Dend+1,:)=O2mLc(Dend:length(HPmWc),1);
else
    Ot3=O2mLc(Dend:Dend+2349);
end
aPrism_Otsnip=[Ot1;gap;Ot2;gap;Ot3];

aPrism_all=[aPrism_HPmWsnip,aPrism_RERsnip,aPrism_COtsnip,aPrism_Otsnip];

fname4=['calo_',mouse]; %=Dex %_pre
eval(['save ',pathout,fname4,'.mat aPrism_all aPrism_HPmWsnip -mat']);


r1=reshape(RER1,[],12);  rs1=nanstd(r1); r1=nanmean(r1);
r2=reshape([gap;RER2],[],12);  rs2=nanstd(r2); r2=nanmean(r2);
r3=reshape([gap;RER3],[],12);  rs3=nanstd(r3); r3=nanmean(r3); 


figure
plot([r1 NaN r2 NaN r3])


