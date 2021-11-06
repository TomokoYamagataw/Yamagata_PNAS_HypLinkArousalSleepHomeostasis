clear all
close all

path='G:\OptoHDbackup\optogenetics\';%path='I:\optogenetics\'; 

% mousenames1=[1 2 3 4 5 6 7 8]; % names of GFP ctrl mice indicated in the file name 
% days1=['260218';'120418';'260618';'040718';'061018';'061018';'061018';'070119'];

%%%%%%%% 10Hz 2min LPO
mousenames2=[5 15 16 17 18 19 20 21];
days2=['020418 030418';'181218 191218';'181218 191218';'181218 191218';...
    '100819 160819';'100819 160819';'070819 180919';'100819 160819'];

mousenames4=[1 6 7 8 9 12]; %non-LPO %GDCh8 '130518 140518'
days4=['130218 140218';'290418 120418';'120518 130518';'100618 150618';...
    '100618 150618';'100618 090618']; 

%%%%%%%% 10Hz 8sec
mousenames8=[16 18 19 20 21];
days8=['181218 060119';'100819 140819';'100819 140819';'070819 140819';'100819 140819']; 



geno=4; %states which animal group (genotype) I want to look at 

before=4;  %we want to analyze the 2min before stimulus onset and 2min after
after=18;   %gives me 2min before, 2min of stim, and the 2min after 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defining all the derivations
ders=strvcat('fro','occ');
der=1; 
deri=ders(der,:);

maxep=21600; %the total no of 4s epochs in 24h period 
epochs=1:maxep; %creates an array with 1 row and columns from 1 to 21600(maxep)
epochs=epochs./900; %edits x such that each value is divided by 900 to convert the values to hours 
            %15 epochs per 1min and 60min per 1h so 900 epochs per 1h
            %%% ./ = right-array division by dividing each element of A by the corresponding element of B
fs=256; 
epochl=4; %states how long the epoch is
zermat=zeros(1,maxep); %creates an array of zeros with 1 row and as long as the total number of epoch

ne=15; %there are 15 epochs in 1min 

x=-before*ne:1:after*ne; %creates an array of epoch indices in the 2min before to after period
x=x*epochl; %converts the epoch indices to seconds 



pathvs1=[path,'outputVSgfp\']; %defines the folder with output files of vigilance states for mice1 (GFP ctrl)
pathvs2=[path,'outputVSchr\']; %defines the folder with output files of vigilance states for mice2 (ChR2 exp)
pathstim=[path,'STIMs\']; %gives the folder with stimulation inf

if geno==1 %geno=1 ctrl group
    mousenames=mousenames1;days=days1; pathvs=pathvs1; gn='GFP'; 
elseif geno==2
    mousenames=mousenames2; days=days2; pathvs=pathvs2; gn='GDCh'; 
elseif geno==4
    mousenames=mousenames4; days=days4; pathvs=pathvs2; gn='GDCh'; 
elseif geno==8
    mousenames=mousenames8; days=days8; pathvs=pathvs2; gn='GDCh'; 
end

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states

numanim=length(mousenames); %how many animals were recorded

aPrisms=[]; aPrisms4=[]; awPrisms=[]; aPrisms1=[]; aPrisms2=[];
mWp=[]; mNp=[]; mRp=[]; mSp=[];
seW=[]; seN=[]; seR=[]; seS=[];
% for pp=1:4
%     if pp==1; tfstim=1; tnstim=18; 
%     elseif pp==2;  tfstim=19; tnstim=36; 
%     elseif pp==3;  tfstim=37; tnstim=54;
%     elseif pp==4;  tfstim=55; tnstim=71;
%     end

tfstim=1; 
tnstim=36;

%%%%%%% calcuate baseline  Light and Dark separately %%% please specify Light or Dark by changing pp later
bin=2; %  

bSp=[]; bswa=[]; wswa=[];
for anim=1:numanim 

mousename=[gn,num2str(mousenames(anim))]; mousename(isspace(mousename))=[];
daysi=days(anim,:); is=find(isspace(daysi)); 

day=daysi(1:is(1)-1); day(isspace(day))=[]; 
fn=[mousename,'-',day,'-',deri];
eval(['load ',pathvs,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);%

VS=zeros(7,maxep); 
VS(1,w)=1; VS(2,w1)=1; VS(3,nr)=1; VS(4,nr2)=1; VS(5,r)=1; VS(6,r3)=1; VS(7,mt)=1;
sleep=sum(VS(3:6,:)); 

sleepp=reshape(sleep,[],bin); sleepp=nanmean(sleepp);


Ns=zeros(1,maxep);        
Ns(1,nr)=1; 
spN=spectr; spN(Ns==0,:)=NaN; 
swa=nanmean(spN(:,3:17),2); swa=swa';  

wswa=[wswa; nanmean(swa)];

swap=reshape(swa,[],bin); 
avrswap=nanmean(swap);

bSp=[bSp; sleepp]; bswa=[bswa; avrswap];
end


%%%%%%% day STIM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for pp=2   %pp=1:Light, pp=2:Dark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SWAs=[];SWAs2=[];
        WallMice=[];NallMice=[];RallMice=[];MallMice=[];SallMice=[];
        for anim=1:numanim 

            mousename=[gn,num2str(mousenames(anim))]; mousename(isspace(mousename))=[];
            daysi=days(anim,:); is=find(isspace(daysi)); 


                day=daysi(is(1):end); day(isspace(day))=[]; 

                fn=[mousename,'-',day,'-',deri];
                eval(['load ',pathvs,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);%

                VS=zeros(7,maxep); 
                VS(1,w)=1; VS(2,w1)=1; VS(3,nr)=1; VS(4,nr2)=1; VS(5,r)=1; VS(6,r3)=1; VS(7,mt)=1;

                wake=sum(VS(1:2,:));     
                nrem=sum(VS(3:4,:));     
                rem=sum(VS(5:6,:)); 
                move=VS(7,:);

                Ns=zeros(1,maxep);        
                Ns(1,nr)=1; 
                spN=spectr; spN(Ns==0,:)=NaN; 
                swa=nanmean(spN(:,3:17),2); swa=swa';            


                fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
                eval(['load ',pathstim,fn0,'.mat startend']);

                stimep1=round(startend(:,1)./(fs*epochl)); % stim start
                stimep2=round(startend(:,2)./(fs*epochl)); % stim end 

                if geno<8;  outs=find((stimep2-stimep1)<25);stimep1(outs)=[];stimep2(outs)=[]; end
                out2s=find((stimep2+after*ne)>21600); stimep1(out2s)=[];stimep2(out2s)=[];
            
                if pp==1
                stimep=stimep1(tnstim*(pp-1)+1:tnstim*pp,1); 
                stimend=stimep2(tnstim*(pp-1)+1:tnstim*pp,1); 
                else
                stimep=stimep1(tnstim*(pp-1)+1:length(stimep1)); 
                stimend=stimep2(tnstim*(pp-1)+1:length(stimep1)); 
                end


                swamat=[];
                numstim2=size(stimep,1);
                for jj=1:numstim2
                    step=stimep(jj); 
                    eps0=step+1:step+after*ne;
                    eps=step-before*ne:step+after*ne-1;
                    %%%%% eps include before, eps0 during+after
                    if min(eps)<1;  
                        swasnip=NaN(1,length(eps));
                    else
                        swasnip=swa(eps);
                    end
                    swamat=[swamat; swasnip];
                end
                SWAs=[SWAs; swamat];               
                SWAs2=[SWAs2; nanmean(swamat)];

                numstim=size(stimep,1) 
                Ws=[]; Ns=[]; Rs=[]; Ms=[]; LmatN=[]; LmatR=[]; LmatW=[];
                for s=1:numstim 
                    step=stimep(s); 
                    eps=step-before*ne:step+after*ne-1; 
                    if min(eps)<1;  
                         ms=NaN(1,length(eps));
                         ws=NaN(1,length(eps));
                         ns=NaN(1,length(eps));
                         rs=NaN(1,length(eps));
                    else
                    ms=move(eps); 
                    ws=wake(eps);   
                    ns=nrem(eps); 
                    rs=rem(eps); 
                    end
                    Ws=[Ws;ws+ms]; 
                    Ns=[Ns;ns];  
                    Rs=[Rs;rs];
                end

                NallMice=[NallMice;nanmean(Ns)];    
                RallMice=[RallMice;nanmean(Rs)];
                SallMice=[SallMice;nanmean(Ns+Rs)];

   
        end
        figure(1)
        plot(nanmean(NallMice),'o-','LineWidth',2)
        hold on
        figure(2)
        plot(nanmean(RallMice),'o-','LineWidth',2)
        hold on
        figure(3)
        plot(nanmean(SallMice),'o-','LineWidth',2)
        hold on                     

        %%%%%%%%% State per 4 sec
        Snorm=100*SallMice./bSp(:,pp);
        
        mNp=nanmean(NallMice);        
        mRp=nanmean(RallMice);
        mSp=nanmean(Snorm);

        seN=nanstd(NallMice)./sqrt(anim);      
        seR=nanstd(RallMice)./sqrt(anim);         
        seS=nanstd(Snorm)./sqrt(anim);
        
        %%%%%%%%% state per 1 min
        Snorm1=Snorm;
        Snorm1=reshape(Snorm1,size(Snorm,1),15,[]);        
        Snorm1=nanmean(Snorm1,2);  %%%% avr 1 min       
        Snorm1=squeeze(Snorm1);

        %%%%%%%%%% SWA per 4sec bin
        normSWA1=SWAs2;  %%%% normalize
        normSWA=100*normSWA1./bswa(:,pp); %%%% normalize

        normswamat=nanmean(normSWA,1)';
        seswamat=nanstd(normSWA,1)/sqrt(size(normSWA,1)); seswamat=seswamat';        

%         wnormSWA=100*normSWA1./wswa;
% 
%         wnormswamat=nanmean(wnormSWA,1)';
%         wseswamat=nanstd(wnormSWA,1)/sqrt(size(wnormSWA,1)); wseswamat=wseswamat';
%  
        %%%%%%%%%% SWA per 1min bin       
        SWAnorm=SWAs2;
        SWAnorm=reshape(SWAnorm,size(SWAnorm,1),15,[]);        
        SWAnorm=nanmean(SWAnorm,2);  %%%% avr 1 min       
        SWAnorm=squeeze(SWAnorm);
        SWAnorm1min=100*SWAnorm./bswa(:,pp);
        

        %%%%%%%%% for 4sec plot
        aPrism_n=anim*ones(length(mSp),1);
        %%%%%%%%% for 1min plot
        aPrism_n2=anim*ones(length(SWAnorm1min),1);

        aPrisms=[aPrisms normswamat seswamat aPrism_n];
%         awPrisms=[awPrisms wnormswamat wseswamat aPrism_n];        
        aPrisms1=[aPrisms1 mNp' mRp'];
        aPrisms2=[aPrisms2 mSp' seS' aPrism_n];
 
        aPrisms_SleepAmount_1min=Snorm1';        
        aPrisms_SWA_1min=SWAnorm1min';           
end
%%%%%%%%%%%%%%%%%%

aPrisms; %time course of swa after stim, normalized with baseline SWA at corresponding time period (Figure3I)
% awPrisms; %time course of swa after stim, normalized with whole baseline SWA(Figure3I)
aPrisms1; %NREM and REM nor normalized
aPrisms2; %time course of the amount of sleep after stim, normalized (Figure3D)
aPrisms_SleepAmount_1min;
aPrisms_SWA_1min;
pause
% aPrism_Wake=mWp';
% aPrism_NREM=mNp';
% aPrism_REM=mRp';
% aPrism_Sleep=mSp';
% aPrism_Wake_SE=seW';
% aPrism_NREM_SE=seN';
% aPrism_REM_SE=seR';
% aPrism_Sleep_SE=seS';
