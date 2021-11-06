
clear all
close all

path='E:\Optoinh\'; 


%%%%%% GFP
mousenames1=[1 2 4 5]; 
SDo1=strvcat('120421','120421','060521','100521');  % SD only
SDs1=strvcat('080421','080421','100521','060521');  % SD+stim

%%%%% Arch
mousenames2=[1 2 4 5 8 9]; % names of mice indicated in the file name 
SDo2=strvcat('080421','120421','080421','100521','100521','060521');  % SD only
SDs2=strvcat('120421','080421','120421','060521','060521','100521');  % SD+stim


genoNs=strvcat('GFP','Arch')

f=0:0.25:30;
maxep=10800; %21600
zermat=zeros(1,maxep);
x=1:1:maxep; x=x./900;

pers=4;
x1=0.5:1:pers;

int=2;

ders=['fro';'occ']
pathin=[path,'outputVS\'];

dr=1;%dr=2 %occ %dr=1 %fro
der=ders(dr,:)

cols=strvcat('db','sb','dr','sr');
% xlab=strvcat('gfpSDo','gfpSDs','archSDo','archSDs')
xlab=strvcat('1','2','3','4')

WbG=[];
NsdG=[];
SLG=[];
SaG=[];

for geno=1:2
    if geno==1
        gn='Gf'; mousenames=mousenames1; SDo=SDo1; SDs=SDs1; genoN=genoNs(geno,:);
    else
        gn='Ar';  mousenames=mousenames2; SDo=SDo2; SDs=SDs2; genoN=genoNs(geno,:);
    end

    numanim=size(mousenames,2);

    Wb=[];
    Nsd=[];
    SL=[];
    Sa=[];

    for n=1:numanim

        mouse=[gn,num2str(mousenames(n))]; mouse(isspace(mouse))=[];

        day=SDo(n,:);
        fn1=[mouse,'-',day,'-epochSO'];
        eval(['load ',pathin,fn1,'.mat epochSO epochSDon epochSDoff -mat']);

        fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
        eval(['load ',pathin,fn,'.mat w nr r w1 nr2 r3 mt ma bastend -mat']);

        W=zermat; W([w; w1])=1;
        N=zermat; N([nr;nr2])=1;

        wb1=sum(W(1:epochSDoff))./15;
        nd1=N(1801:epochSO); sg1=sum(N(epochSO:epochSO+int*900))./15;
        re=rem(length(nd1),pers); nd1(length(nd1)-re+1:end)=[];
        nd1=sum(reshape(nd1,length(nd1)/pers,pers))./15;
        sl1=(epochSO-epochSDoff)./15;

        day=SDs(n,:);
        fn2=[mouse,'-',day,'-epochSO'];
        eval(['load ',pathin,fn2,'.mat epochSO epochSDon epochSDoff -mat']);

        fn=[mouse,'-',day,'-',der];%fn=[mouse,'-',day,'-',der,'-VSspec'];
        eval(['load ',pathin,fn,'.mat w nr r w1 nr2 r3 mt ma bastend -mat']);

        W=zermat; W([w; w1])=1;
        N=zermat; N([nr;nr2])=1;

        wb2=sum(W(1:epochSDoff))./15;
        nd2=N(1801:epochSO); sg2=sum(N(epochSO:epochSO+int*900))./15;
        re=rem(length(nd2),pers); nd2(length(nd2)-re+1:end)=[];
        nd2=sum(reshape(nd2,length(nd2)/pers,pers))./15;
        sl2=(epochSO-epochSDoff)./15;

        Wb=[Wb;wb1 wb2];
        Nsd=[Nsd;nd1 nd2];
        %Nsd=[Nsd;cumsum(nd1) cumsum(nd2)];
        SL=[SL;sl1 sl2];
        Sa=[Sa;sg1 sg2];  
        
        if geno==1
                WbG1=Wb;
                SLG1=SL;
                SaG1=Sa;
                Nsd1=Nsd; 
        else
                WbG2=Wb;
                SLG2=SL;
                SaG2=Sa;
                Nsd2=Nsd; 
        end
    end
  
end
WbG1=WbG1';
WbG2=WbG2';
SLG1=SLG1';
SLG2=SLG2';
SaG1=SaG1';
SaG2=SaG2';

aPrismSleepLatency=[SLG2 SLG1];
aPrismSleepStability=[SaG2 SaG1];

%      WbG = [WbG1 WbG2]; WbG=WbG'
%      SLG = [SLG1 SLG2]; SLG=SLG'
%      SaG = [SaG1 SaG2]; SaG=SaG'
     
for jj=1:2
    figure
    if jj==1; WbG=WbG1; SLG=SLG1; SaG=SaG1; else;  WbG=WbG2; SLG=SLG2; SaG=SaG2; end
    subplot(2,2,1)
    for i=1:2 scatter(ones(1,2)*i,WbG(:,i),[cols(i,:)],'filled'); hold on; end
    hold on
    boxplot(WbG,'notch','on');
    set(gca,'XTickLabel',xlab);
    grid on
    xlabel([])
    title('Wake before SD')

    subplot(2,2,2)
    for i=1:2 scatter(ones(1,2)*i,SLG(:,i),[cols(i,:)],'filled'); hold on; end
    hold on
    boxplot(SLG,'notch','on');
    set(gca,'XTickLabel',xlab);
    xlabel([])
    grid on
    title('Sleep latency')

    m1=nanmean(Nsd1); s1=nanstd(Nsd1)./sqrt(size(Nsd1,1));
    m1=reshape(m1,pers,2);s1=reshape(s1,pers,2);

    m2=nanmean(Nsd2); s2=nanstd(Nsd2)./sqrt(size(Nsd2,1));
    m2=reshape(m2,pers,2);s2=reshape(s2,pers,2);

    subplot(2,3,4)
    errorbar(x1,m1(:,1),s1(:,1),'d-b','LineWidth',2);
    hold on
    errorbar(x1,m1(:,2),s1(:,2),'s-r','LineWidth',2);
    axis([0 4 0 2])
    grid on
    xlabel('Hours')
    ylabel('min')
    title('Sleep during SD: GFP mice')

    subplot(2,3,5)
    errorbar(x1,m2(:,1),s2(:,1),'d-b','LineWidth',2);
    hold on
    errorbar(x1,m2(:,2),s2(:,2),'s-r','LineWidth',2);
    axis([0 4 0 2])
    grid on
    xlabel('Hours')
    ylabel('min')
    title('Sleep during SD: Arch mice')

    subplot(2,3,6)
    for i=1:2 scatter(ones(1,2)*i,SaG(:,i),[cols(i,:)],'filled'); hold on; end
    hold on
    boxplot(SaG,'notch','on');
    set(gca,'XTickLabel',xlab);
    grid on
    xlabel([])
    title('Sleep after SD')

end

% hold on
% errorbar(f,SPbefore(:,3),SPbefore(:,4),'d-r','LineWidth',2)
% axis([0 30 80 140])
% grid on
%
% subplot(1,2,2)
% errorbar(f,SPafter(:,1),SPafter(:,2),'d-b','LineWidth',2)
% hold on
% errorbar(f,SPafter(:,3),SPafter(:,4),'d-r','LineWidth',2)
% axis([0 30 70 150])
% grid on
% xlabel('Frequency (Hz)')


% 
% 
% subplot(2,4,1)
% for i=1:2 scatter(ones(1,2)*i,WbG1(:,i),[cols(i,:)],'filled'); hold on; end
% hold on
% boxplot(WbG1,'notch','on');
% set(gca,'XTickLabel',xlab);
% grid on
% xlabel([])
% title('Wake before SD GFP')
% 
% subplot(2,4,2)
% for i=1:2 scatter(ones(1,2)*i,WbG2(:,i),[cols(i,:)],'filled'); hold on; end
% hold on
% boxplot(WbG2,'notch','on');
% set(gca,'XTickLabel',xlab);
% grid on
% xlabel([])
% title('Wake before SD ChR2')