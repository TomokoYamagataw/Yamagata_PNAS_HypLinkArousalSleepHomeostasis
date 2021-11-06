
clear all
close all

path='E:\Optoinh\';

% %%%%% Arch
% gn='Ar'; 
% mousenames=[1 1 2 2 4 4 5 5 8 8 9 9]; % names of mice indicated in the file name 
% days=strvcat('080421','120421','120421','080421','080421','120421','100521','060521','100521','060521','060521','100521');
% epochSOs=[3640 4003 4010 3984 3697 3602 3813 3678 4099 3468 4203 3496]; %3398 for Ar9 100521?
% epochSDons=[1800 1800 1800 1800 1800 1800 1629 1681 1629 1681 1681 1629]; 
% epochSDoffs=[3580 3580 3580 3580 3580 3580 3409 3461 3409 3461 3461 3409]; %3398 for Ar9 100521?


%%%%% GFP
gn='Gf';
mousenames=[1 1 2 2 4 4 5 5]; 
days=strvcat('120421','080421','120421','080421','060521','100521','100521','060521');
epochSOs=[3698 3584 3627 3704 3576 3425 3794 3657];
epochSDons=[1800 1800 1800 1800 1681 1629 1629 1681]; 
epochSDoffs=[3580 3580 3580 3580 3461 3409 3409 3461]; 

f=0:0.25:30;
maxep=10800; %maxep=21600;
zermat=zeros(1,maxep);
x=1:1:maxep; x=x./900;

ders=['fro';'occ']
pathin=[path,'outputVS\'];

numanim=size(mousenames,2);
OOO=[];
dr=1 %dr=2
der=ders(dr,:)

for n=1:numanim

    figure
    mousename=[gn,num2str(mousenames(n))]; mousename(isspace(mousename))=[];
    day=num2str(days(n,:));
    epochSO=epochSOs(n)
    epochSDon=epochSDons(n)   
    epochSDoff=epochSDoffs(n)
    
    fname=[mousename,'-',day,'-',der]

    fn=[mousename,'-',day,'-',der]; %fn=[mouse,'-',day,'-',der,'-VSspec'];
    eval(['load ',pathin,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);

    swa=mean(spectr(:,3:17),2);
    W=zermat; W(w)=1;
    N=zermat; N(nr)=1;
    R=zermat; R(r)=1;

    swaW=swa; swaW(W==0)=NaN;
    swaN=swa; swaN(N==0)=NaN;
    swaR=swa; swaR(R==0)=NaN;

    swa=[swaW swaN swaR];
    
    swa=swa(1:900*6,:)
    plot(swa);
    hold on
    plot([epochSO epochSO],[0 5000],'-m','LineWidth',2) 
    plot([epochSDon epochSDon],[0 5000],'-b','LineWidth',1) 
    plot([epochSDoff epochSDoff],[0 5000],'-c','LineWidth',1) 
    axis([0 5400 0 5000])
    title (fn); 
    
    pause
    close all
    
    fn1=[mousename,'-',day,'-epochSO'];
    eval(['save ',pathin,fn1,'.mat epochSO epochSDon epochSDoff -mat']);
    

end