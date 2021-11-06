
clear all
close all
path='I:\optogenetics\';

% mousenames=strvcat('GDCh1','GDCh1','GDCh4','GDCh4','GDCh5','GDCh5','GDCh6','GDCh6','GFP1','GFP1','GFP2','GFP2','GFP3','GFP3','GFP4','GFP4');
% days=strvcat('040518','060518','140518','160518','040518','060518','020518','300418','140518','160518','020518','300418','010718','030718','010718','030718');
% epochSOs=[3607 3696 3675 3606 4016 3965 3549 3807 3504 3669 3620 3948 3683 3608 3841 3810];

% mousenames=strvcat('GDCh7','GDCh7','GDCh8','GDCh8','GDCh9','GDCh9','GDCh12','GDCh12','GDCh16','GDCh16');
% days=strvcat('030718','010718','130618','110618','130618','110618','130618','110618','090119','110119');
% epochSOs=[3789 3687 3765 3556 3939 3734 4049 3758 3907 3610]; %3772
% 
% mousenames=strvcat('GDCh17','GDCh17');
% days=strvcat('090119','110119');
% epochSOs=[3953 3724];
% dr=2;
% 
% mousenames=strvcat('GDCh18','GDCh18','GDCh19','GDCh19','GDCh20','GDCh20','GDCh21','GDCh21');
% days=strvcat('080819','110819','080819','110819','080819','110819','080819','110819');
% epochSOs=[3933 3834 3668 4017 3628 3764 3705 4015];
% epochSOs=[3933 3834 3668 4017 3613 3662];

%%%%%% caffeine
mousenames=strvcat('GDCh18','GDCh18','GDCh19','GDCh19','GDCh20','GDCh20','GDCh21','GDCh21');
days=strvcat('250819','280819','280819','250819','250819','280819','280819','250819');
epochSOs=[3696 3646 4158 3787 5295 4350 3817 4017];
% epochSOs=[3825 3646 4158 3787 5295 4350 3817 4017];


% mousenames=strvcat('GFP5','GFP5','GFP6','GFP6','GFP7','GFP7','GFP8','GFP8')
% days=strvcat('041018','011018','111018','091018','111018','091018','090119','110119')
% epochSOs=[3728 3777 3606 3676 3637 3669 3618 3645];

%%%% Left:SD+stim Rigth:SD only
%%%%%%% Wake Enhancement ChR2
% mousenames=strvcat('GDCh1','GDCh1','GDCh4','GDCh4','GDCh5','GDCh5','GDCh6','GDCh6','GDCh7','GDCh7','GDCh8','GDCh8','GDCh9','GDCh9','GDCh12','GDCh12','GDCh16','GDCh16')
% days=strvcat('040518','060518','160518','140518','040518','060518','300418','020518','030718','010718','130618','110618','130618','110618','130618','110618','090119','110119')
% mousenames=strvcat('GDCh17','GDCh17') %%occ
% days=strvcat('090119','110119') %%occ
% der=strvcat('occ_30Hz'); 
% ders=strvcat('occ')


f=0:0.25:30;
maxep=10800; %maxep=21600;
zermat=zeros(1,maxep);
x=1:1:maxep; x=x./900;

ders=['fro';'occ']
% pathin=[path,'outputVS1\WakeEnhance\'];
pathin=[path,'outputVSchr\']; %

numanim=size(mousenames,1);
OOO=[];
dr=1 %dr=2
der=ders(dr,:)

for n=1:numanim

    figure
    mouse=mousenames(n,:); mouse(isspace(mouse))=[];
    day=days(n,:);
    epochSO=epochSOs(n)
    
    fname=[mouse,'-',day,'-',der]

    fn=[mouse,'-',day,'-',der]; %fn=[mouse,'-',day,'-',der,'-VSspec'];
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
    axis([0 5400 0 5000])
    title (fn); 
    
    pause
    close all
    
    fn1=[mouse,'-',day,'-epochSO'];
    eval(['save ',pathin,fn1,'.mat epochSO -mat']);
    

end