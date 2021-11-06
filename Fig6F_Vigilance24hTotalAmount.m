
clear all
close all

path='E:\Optoinh\'; 

% Arch
mousenames1=[1 2 4 5 8 9]; % names of mice indicated in the file name 
days1=['070421 160421';'070421 150421';'070421 160421';'050521 140521';'050521 140521';'050521 150521'];

% GFP
mousenames2=[1 2 4 5]; 
days2=['070421 150421';'070421 150421';'050521 150521';'050521 140521'];


ders=strvcat('fro','occ');der=1;deri=ders(der,:);

maxep=21600;x=1:maxep; x=x./900;epochl=4;
zermat=zeros(1,maxep);

pathvs=[path,'outputVS\'];

vsname=strvcat('Wake','NREM','REM','SWA')

int=12;
numint=24/int;
numh=900;
x=int/2:int:24;

SWAs=[]; SWAm=[]; aW=[]; aN=[]; aR=[]; aS=[]; mean_sd_wnr=[]; wpt=[]; npt=[]; rpt=[];rspt=[]; spt=[];  msdper=[];
for geno=1:2
if geno==1
    mousenames=mousenames1;days=days1; gn='Ar'; 
elseif geno==2  %geno=2 ctrl group
    mousenames=mousenames2; days=days2; gn='Gf';  
end
    numanim=length(mousenames);

    figure;

    WNR=[];
    wgeno=[]; ngeno=[]; rgeno=[]; sgeno=[]; swageno=[];

    for anim=1:numanim
        mousename=[gn,num2str(mousenames(anim))];
        daysi=days(anim,:); is=find(isspace(daysi));

        SWA1=[]; W1=[];N1=[];R1=[]; S1=[];
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
            sep=sum(VS,1); sep=sum(reshape(sep,numh*int,numint));

            wake=sum(VS(1:2,:));nrem=sum(VS(3:4,:));rems=sum(VS(5:6,:));

            % SWA
            swa=mean(spectr(:,3:17),2);swa=swa'; swaN=swa; swaN(VS(3,:)==0)=NaN;
            swa=nanmean(reshape(swaN,numh*int,numint));
               
            w=sum(reshape(wake,numh*int,numint)); %L D
            n=sum(reshape(nrem,numh*int,numint)); %L D
            r=sum(reshape(rems,numh*int,numint)); %L D
            ts=n+r; %rps=r./s*100;
            
            
            W1=[W1 w./sep*100];N1=[N1 n./sep*100];R1=[R1 r./sep*100]; S1=[S1 ts./sep*100]; %%%sham-L, sham-D, stim-L, stim-D

            SWA1=[SWA1 swa];
            
        end

        %W1=W1./W1(numint)*100;N1=N1./N1(numint)*100;R1=R1./R1(numint)*100;
        
        WNR=[WNR;wnr];

        wgeno=[wgeno;W1];
        ngeno=[ngeno;N1];
        rgeno=[rgeno;R1];
        sgeno=[sgeno;S1];
        
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
    
    aW=[aW; wgeno]; aN=[aN; ngeno]; aR=[aR; rgeno]; aS=[aS; sgeno]; 
    mean_sd_wnr=[mean_sd_wnr; m1 s1 m2 s2 m3 s3];
    
    wp1=wgeno(:,3)./wgeno(:,1);  wp2=wgeno(:,4)./wgeno(:,2);
    np1=ngeno(:,3)./ngeno(:,1);  np2=ngeno(:,4)./ngeno(:,2);
    rp1=rgeno(:,3)./rgeno(:,1);  rp2=rgeno(:,4)./rgeno(:,2);
    
    rps1=(rgeno(:,3)./rgeno(:,1))./(sgeno(:,3)./sgeno(:,1)); 
    rps2=(rgeno(:,4)./rgeno(:,2))./(sgeno(:,4)./sgeno(:,2));
    
    sp1=sgeno(:,3)./sgeno(:,1); sp2=sgeno(:,4)./sgeno(:,2);
    
      wpt=[wpt; wp1 wp2];
      npt=[npt; np1 np2];
      rpt=[rpt; rp1 rp2];
      rspt=[rspt; rps1 rps2];
      spt=[spt; sp1 sp2];
      
    m5=nanmean(wp1)*100; s5=nanstd(wp1)./sqrt(numanim);    m9=nanmean(wp2)*100; s9=nanstd(wp2)./sqrt(numanim); 
    m6=nanmean(np1)*100; s6=nanstd(np1)./sqrt(numanim);    m10=nanmean(np2)*100; s10=nanstd(np2)./sqrt(numanim); 
    m7=nanmean(rp1)*100; s7=nanstd(rp1)./sqrt(numanim);    m11=nanmean(rp2)*100; s11=nanstd(rp2)./sqrt(numanim); 
    m8=nanmean(rps1)*100; s8=nanstd(rps1)./sqrt(numanim);  m12=nanmean(rps2)*100; s12=nanstd(rps2)./sqrt(numanim);
    
    m21=nanmean(sp1)*100; s21=nanstd(sp1)./sqrt(numanim);  m22=nanmean(sp2)*100; s22=nanstd(sp2)./sqrt(numanim);
    
    
    msdper=[msdper; m5 s5 m9 s9 m6 s6 m10 s10 m7 s7 m11 s11 m8 s8 m12 s12 m21 s21 m22 s22];
%     %pause
%     SWAs=[SWAs swageno(:,1:12)' swageno(:,13:24)'];
%     SWAm=[SWAm; nanmean(swageno(:,1:12)) nanmean(swageno(:,13:24))];
end
aPrism_AmountWakeLD=aW';  %Arch(n=6) L-base D-base L-stim D-stim;
                       %GFP(n=4) L-base D-base L-stim D-stim;
aPrism_AmountNREMLD=aN';                        
aPrism_AmountREMLD=aR'; 
aPrism_AmountSleeppLD=aS'; 

aPrism_wnr_persham=[wpt*100 npt*100 rpt*100 rspt*100 spt*100]; 
% aPrism_wnr_persham=aPrism_wnr_persham';
s_wnrpersham=nanmean(wgeno); s_wnrpersham=nanstd(wgeno)./sqrt(numanim);

a_mean_sd_wnr_percentage=msdper %% geno=1:Arch, geno=2:GFP
mean_sd_wnr
pause
% aPrism_SWA=SWAs;
% aPrism_SWA_mean=SWAm;