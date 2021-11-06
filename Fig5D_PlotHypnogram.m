
clear all
close all
path='I:\optogenetics\';
pathF=[path,'Figures\WakeEnhance\'];
% pathin=[path,'outputVS1\WakeEnhance\']; %
pathin=[path,'outputVSchr\']; 

% mousenames=strvcat('GDCh1','GDCh1','GDCh4','GDCh4','GDCh5','GDCh5','GDCh6','GDCh6','GFP1','GFP1','GFP2','GFP2');
% days=strvcat('040518','060518','140518','160518','040518','060518','020518','300418','140518','160518','020518','300418');
% ders=['fro';'occ']

%%%% Left:SD+stim Rigth:SD only
%%%%%%% Wake Enhancement ChR2
% mousenames=strvcat('GDCh1','GDCh1','GDCh4','GDCh4','GDCh5','GDCh5','GDCh6','GDCh6','GDCh7','GDCh7','GDCh8','GDCh8','GDCh9','GDCh9','GDCh12','GDCh12','GDCh16','GDCh16')
% days=strvcat('040518','060518','160518','140518','040518','060518','300418','020518','030718','010718','130618','110618','130618','110618','130618','110618','090119','110119')
% ders=['fro';'occ']

%%%%%% occ only %%%%%
% mousenames=strvcat('GDCh17','GDCh17') %%occ
% days=strvcat('090119','110119') %%occ
% der=strvcat('occ'); 
% ders=strvcat('occ')
% dr=2; %1%2

% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
% days=strvcat('080819','110819','080819','110819','080819','110819','080819','110819');
% ders=['fro';'occ']

%%%%%%%%%% Wake Enhancement GFP
% mousenames=strvcat('GFP1','GFP1','GFP2','GFP2','GFP3','GFP3','GFP4','GFP4','GFP5','GFP5','GFP6','GFP6','GFP7','GFP7','GFP8','GFP8')
% days=strvcat('160518','140518','300418','020518','010718','030718','030718','010718','041018','011018','111018','091018','111018','091018','090119','110119')
% ders=['fro';'occ']

%%%%%%% caffeine
mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
days=strvcat('250819','280819','280819','250819','250819','280819','280819','250819');
ders=['fro';'occ']

f=0:0.25:30;
maxep=10800; %21600
zermat=zeros(1,maxep);
x=1:1:maxep; x=x./900;

numanim=size(mousenames,1);
OOO=[];
dr=1; %1%2
der=ders(dr,:)

figure
for ii=1:size(mousenames,1)
    
    for n=(ii-1)*2+1:ii*2 %1:4 %5:
    
    mouse=mousenames(ii,:); mouse(isspace(mouse))=[];
    day=days(n,:);

    fname=[mouse,'-',day,'-',der]

    fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
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
    subplot(2,1,floor(n/ii))
    plot(swa);
    axis([0 5400 0 5000])
    %plot(swa)

    title (fn); 
    
    end

%     orient tall
%     saveas(gcf,[pathF,['Hypnogram_',der,'_',mouse]],'fig')
    saveas(gcf,[pathF,['Hypnogram_',der,'_',mouse,'_caffeine']],'fig')
end

