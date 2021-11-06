
clear all
close all

path='I:\optogenetics\';

%%%% for martin, Please note that I excluded GFP5 on Prism, but not here. 
% mousenames1=[1 2 3 4 5 6 7 8]; %GFP controls % GFP4: data back to 040718 from 260618 % GFP5: 031018 to 051018
% days1=['250218 260218';'290418 120418';'300618 260618';'300618 040718';'031018 061018';'081018 061018';'081018 061018';'080119 070119'];

mousenames1=[1 2 3 4 6 7 8]; %GFP controls % GFP4: data back to 040718 from 260618 % GFP5: 031018 to 051018
days1=['250218 260218';'290418 120418';'300618 260618';'300618 040718';'081018 061018';'081018 061018';'080119 070119'];

%%%%% LPO
mousenames2=[5 15 16 17 18 19 20 21]; %ChR2 mice %GDCh8: 050618 to 150618
days2=['020418 030418';...
    '181218 191218';'181218 191218';'181218 191218';...
    '100819 160819';'100819 160819';'070819 180919';'100819 160819'];	%160819

%%%% nonLPO
mousenames3=[1 6 7 8 9 12]; %ChR2 mice %GDCh8: 050618 to 150618
days3=['130218 140218';'290418 120418';'120518 130518';'100618 150618';'100618 150618';'100618 090618'];

%%%%% 8sec
% mousenames2=[16 18 19 20 21];
% days2=['181218 060119';'100819 140819';'100819 140819';'070819 140819';'100819 140819']; 

%%%% original
% mousenames2=[5 6 8 12 15 16 17 18 19 20 21]; %ChR2 mice %GDCh8: 050618 to 150618
% days2=['020418 030418';'290418 120418';'130518 140518';'100618 090618';...
%     '181218 191218';'181218 191218';'181218 191218';...
%     '100819 160819';'100819 160819';'070819 180919';'100819 160819'];	%160819
% 

ders=strvcat('fro','occ');der=1;deri=ders(der,:);

maxep=21600;x=1:maxep; x=x./900;epochl=4;
zermat=zeros(1,maxep);

pathvs1=[path,'outputVSgfp\'];
pathvs2=[path,'outputVSchr\'];

vsname=strvcat('Wake','NREM','REM','SWA')

int=2;
numint=24/int;
numh=900;
x=int/2:int:24;

SWAs=[]; SWAm=[];
for geno=1:3
    if geno==1; mousenames=mousenames1; days=days1; pathvs=pathvs1; gn='GFP'; 
    elseif geno==2; mousenames=mousenames2; days=days2; pathvs=pathvs2; gn='GDCh'; 
    elseif geno==3; mousenames=mousenames3; days=days3; pathvs=pathvs2; gn='GDCh';         
    end
    numanim=length(mousenames);

    figure;

    WNR=[];

    wgeno=[];
    ngeno=[];
    rgeno=[];
    swageno=[];
        
    for anim=1:numanim
        mousename=[gn,num2str(mousenames(anim))];
        daysi=days(anim,:); is=find(isspace(daysi));
        
        SWA1=[];
        W1=[];
        N1=[];
        R1=[];

        wnr=[];

        for dd=1:2

            if dd==1 day=daysi(1:is(1)-1); else day=daysi(is(1):end); end
            day(isspace(day))=[];

            fn1=[mousename,'-',day,'-',deri];
            eval(['load ',pathvs,fn1,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);

            VS=zeros(7,maxep);
            VS(1,w)=1;VS(2,w1)=1;VS(3,nr)=1;VS(4,nr2)=1;VS(5,r)=1;VS(6,r3)=1;VS(7,mt)=1;
            clear w nr r w1 nr2 r3 mt ma bastend;

            % 24 hour amount of vigilance states
            te=sum(sum(VS,2));
            wake=sum(sum(VS(1:2,:)));nrem=sum(sum(VS(3:4,:)));rems=sum(sum(VS(5:6,:)));
            wnr=[wnr wake/te*100 nrem/te*100 rems/te*100];

            % time course
            te=sum(VS,1); te=sum(reshape(te,numh*int,numint));

            wake=sum(VS(1:2,:));nrem=sum(VS(3:4,:));rems=sum(VS(5:6,:));

            % SWA
            swa=mean(spectr(:,3:17),2);swa=swa'; swaN=swa; swaN(VS(3,:)==0)=NaN;
            swa=nanmean(reshape(swaN,numh*int,numint));
               
            w=sum(reshape(wake,numh*int,numint));
            n=sum(reshape(nrem,numh*int,numint));
            r=sum(reshape(rems,numh*int,numint));

            W1=[W1 w./te*100];N1=[N1 n./te*100];R1=[R1 r./te*100];
            SWA1=[SWA1 swa];
            
        end

        %W1=W1./W1(numint)*100;N1=N1./N1(numint)*100;R1=R1./R1(numint)*100;
        
        WNR=[WNR;wnr];

        wgeno=[wgeno;W1];
        ngeno=[ngeno;N1];
        rgeno=[rgeno;R1];
        swageno=[swageno;SWA1./nanmean(SWA1(1:numint))*100];

    end

    m1=nanmean(wgeno); s1=nanstd(wgeno)./sqrt(numanim); m1=reshape(m1,numint,2); s1=reshape(s1,numint,2);
    m2=nanmean(ngeno); s2=nanstd(ngeno)./sqrt(numanim); m2=reshape(m2,numint,2); s2=reshape(s2,numint,2);
    m3=nanmean(rgeno); s3=nanstd(rgeno)./sqrt(numanim); m3=reshape(m3,numint,2); s3=reshape(s3,numint,2);
    m4=nanmean(swageno); s4=nanstd(swageno)./sqrt(numanim);m4=reshape(m4,numint,2); s4=reshape(s4,numint,2);

    for vs=1:4
        if vs==1 m=m1; s=s1; elseif vs==2 m=m2; s=s2; elseif vs==3 m=m3; s=s3; else; m=m4; s=s4; end
        subplot(1,4,vs)
        errorbar(x,m(:,1),s(:,1),'o-b','LineWidth',2);
        hold on
        errorbar(x,m(:,2),s(:,2),'o-r','LineWidth',2);
        %axis([0 24 0 120])
        set(gca,'XTick',[0:4:24])
        grid on
        title(vsname(vs,:));
    end
    
    %pause
    SWAs=[SWAs swageno(:,1:12)' swageno(:,13:24)'];
    SWAm=[SWAm; nanmean(swageno(:,1:12)) nanmean(swageno(:,13:24))];
end
aPrism_SWA=SWAs;
aPrism_SWA_mean=SWAm;