
clear all
close all

path='E:\Optoinh\'; 
pathF=[path,'FiguresWakeSuppression\']; mkdir(pathF)
pathin=[path,'outputVS\']; 


%%% Left:SD+stim Rigth:SD only
%%%%%% Wake Suppression Arch

% % % % Arch
geno=1;
mousenames1=[1 2 4 5 8 9];
days1=['080421 120421';'120421 080421';'080421 120421';'100521 060521';'100521 060521';'060521 100521'];

% % % % GFP
% geno=2;
% mousenames2=[1 2 4 5]; 
% days2=['120421 080421';'120421 080421';'060521 100521';'100521 060521'];

pathvs=[path,'outputVS\']; %defines the folder with output files of vigilance states for mice2 (ChR2 exp)
% pathstim=[path,'STIMs\']; 

if geno==1 %geno=1 ctrl group
    mousenames=mousenames1;days=days1; gn='Ar'; %firststims=firststims1;
elseif geno==2
    mousenames=mousenames2; days=days2; gn='Gf';  %firststims=firststims2;
end 

% mousenames=strvcat('GDCh1','GDCh1','GDCh4','GDCh4','GDCh5','GDCh5','GDCh6','GDCh6','GDCh7','GDCh7','GDCh8','GDCh8','GDCh9','GDCh9','GDCh12','GDCh12','GDCh16','GDCh16')
% days=strvcat('040518','060518','160518','140518','040518','060518','300418','020518','030718','010718','130618','110618','130618','110618','130618','110618','090119','110119')

ders=['fro';'occ']
numanim=length(mousenames);

f=0:0.25:30;
maxep=10800; %21600
zermat=zeros(1,maxep);
x=1:1:maxep; x=x./900;

OOO=[];
dr=1; %1%2
der=ders(dr,:)

for ii=1:numanim
    
    mousename=[gn,num2str(mousenames(ii))]; mousename(isspace(mousename))=[];
    daysi=days(ii,:); is=find(isspace(daysi));
    
    for dd=1:2 %1:4 %5:
        
        if dd==1 day=daysi(1:is(1)-1); else day=daysi(is(1):end); end
        day(isspace(day))=[];

    fname=[mousename,'-',day,'-',der]

    fn=[mousename,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
    eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);

    swa=mean(spectr(:,3:17),2);
    W=zermat; W(w)=1;
    N=zermat; N(nr)=1;
    R=zermat; R(r)=1;

    swaW=swa; swaW(W==0)=NaN;
    swaN=swa; swaN(N==0)=NaN;
    swaR=swa; swaR(R==0)=NaN;

    figure(ii)
    swa=[swaW swaN swaR];
    subplot(2,1,floor(dd))
    plot(swa);
    axis([0 5400 0 7000])


    title (fn); 
    
    end

%     orient tall
    saveas(gcf,[pathF,['Hypnogram_',der,'_',mousename]],'fig')
end

