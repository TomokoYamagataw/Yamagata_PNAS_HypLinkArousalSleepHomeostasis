close all
clear all


path='G:\OptoHDbackup\optogenetics\'; %
pathstim = [path,'STIMs\']; %pathstim ='E:\Optoinh\STIMs\'; % pathStim = 'I:\optogenetics\STIMs';
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 240 228 66;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue,yellow
npgCMap = [npgCMap;npgCMap];
pathFigs = [path,'Figures\']; %pathFigs = 'I:\optogenetics\finalFigures\';
period1=30; period2=30; %0 105:7min 
LD=1;

%%%% for martin, Please note that I excluded GFP5 on Prism, but not here.
mousenames1=[1 2 3 4 5 6 7 8]; %GFP controls % GFP4: 040718 to 260618 % GFP5: 031018 to 051018
days1=['250218 260218';'290418 120418';'300618 260618';'300618 040718';'031018 061018';'081018 061018';'081018 061018';'080119 070119'];
% mousenames1=[1 2 3 4 6 7 8]; %GFP controls % GFP4: 040718 to 260618 % GFP5: 031018 to 051018
% days1=['250218 260218';'290418 120418';'300618 260618';'300618 040718';'081018 061018';'081018 061018';'080119 070119'];

%%%% LPO
mousenames2=[5 15 16 17 18 19 20 21];
days2=['020418 030418';'181218 191218';'181218 191218';'181218 191218';...
    '100819 160819';'100819 160819';'070819 180919';'100819 160819'];	
%%%%% nonLPO
mousenames3=[1 6 7 8 9 12]; %ChR2 mice %GDCh8: 050618 to 150618
% % days3=['130218 140218';'290418 120418';'120518 130518';'100618 150618';'100618 150618';'100618 090618'];
days3=['130218 140218';'290418 120418';'120518 130518';'100618 150618';'100618 150618';'100618 090618'];



% % % % %%%%% GFP controls 
% % % % mousenames1=[1 2 4 5]; 
% % % % days1=['070421 150421';'070421 150421';'050521 150521';'050521 140521'];
% % % % %%%%% Arch
% % % % mousenames2=[1 2 4 5 8 9];
% % % % days2=['070421 160421';'070421 150421';'070421 160421';'050521 140521';'050521 140521';'050521 150521'];


ders=strvcat('fro','occ');der=1;deri=ders(der,:);
groups = {'GFP','LPO','nLPO'};%groups = {'GFP','Arch'}; %
maxep=21600;x=1:maxep; x=x./900;epochl=4; fs=256; 
zermat=zeros(1,maxep);


vsname=strvcat('Wake','NREM','REM','SWA');
pathvs=[path,'outputVS\'];% pathvs=[path,'outputVSall\'];% relic from old cold can be replaced
int=2;
numint=24/int;
numh=900; % num epochs in 1h


%%%%%%%%  
MOVEinclude = 1
%%%%%%%%  

for geno=1:3 %1:numel(groups)
    if geno==1
        mousenames=mousenames1;days=days1;gn='GFP';pathvs=[path,'outputVSgfp\'];% ;
    elseif geno ==2
        mousenames=mousenames2;days=days2;gn='GDCh';pathvs=[path,'outputVSchr\'];% 
    else
        mousenames=mousenames3;days=days3;gn='GDCh';pathvs=[path,'outputVSchr\'];% 
    end
    numanim=length(mousenames);
     
    
    aPrismWW=[]; aPrismWN=[]; aPrismWR=[]; aPrismNW=[]; aPrismNN=[]; aPrismNR=[]; aPrismRW=[]; aPrismRN=[]; aPrismRR=[];   
    for dd=[2 1]


        [wakeSpectraFro,wakeSpectraOcc,wDurs,wOnsets,wakeSpecGramOcc,wakeSpecGramFro] = deal(cell(numanim,2));
        VVmean=cell(1,numanim);

        MV=[]; TS=[]; nTS=[];   mW=[]; mN=[]; mR=[];
        for anim=1:numanim
            mousename=[gn,num2str(mousenames(anim))];
            daysi=days(anim,:); is=find(isspace(daysi));

            [SWA1,W1,N1,R1,SWA1_cs,W1_cs,mt_cs,N1_cs,R1_cs,SWA_cs]=deal([]);


            wnr=[];

            if dd==1; day=daysi(1:is(1)-1); 
            else; day=daysi(is(1):end); 
            end
            day(isspace(day))=[];
            
            
            if dd==2
            fn0=[mousename,'_',day,'_stim']; 
            eval(['load ',pathstim,fn0,'.mat startend ']);

            stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
            stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2) 

            outs=find((stimep2-stimep1)<25);
            stimep1(outs)=[];stimep2(outs)=[];
            end

%                 if LD==1
%                     stimep1(stimep1>=10800)=[]; stimep2(stimep2>=10800)=[];
%                 elseif LD==2
%                     stimep1(stimep1<10800)=[]; stimep2(stimep2<10800)=[];
%                 else
%                 end

                % load recording
                fn1=[mousename,'-',day,'-',deri];
                eval(['load ',pathvs,fn1,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
                currFro = spectr;


                % get vigilance states
                VS=zeros(7,maxep);
                VS(1,w)=1;VS(2,w1)=1;VS(3,nr)=1;VS(4,nr2)=1;VS(5,r)=1;VS(6,r3)=1;VS(7,mt)=1;
                clear w nr r w1 nr2 r3 mt ma bastend;
                VS(:,21600)=0;

                wake=sum(VS(1:2,:));nrem=sum(VS(3:4,:));rems=sum(VS(5:6,:)); move = VS(7,:);
                
                %%%%%%%%%% choose whether analysis includes move or not
                if MOVEinclude==1; wakm=sum(VS([1:2 7],:)); elseif MOVEinclude==0; wakm=wake; end
%               

                Vmat=zeros(3,maxep);
                Vmat(1,:)=wakm;
                Vmat(2,:)=nrem;
                Vmat(3,:)=rems;
                
                TSmat=sum(Vmat,2); TS=[TS TSmat];
                nTSmat=sum(TS); nTS=[nTS nTSmat];
                MVmat=TSmat./nTSmat; %MVmat=mean(Vmat,2)*100;
                MV=[MV MVmat];

                numstim=size(stimep1,1);

                Ws=[]; Ns=[]; Rs=[]; Ws2=[]; Ns2=[]; Rs2=[]; % Ms=[];

                ne=15; %there are 15 epochs in 1min 
                before=2;  
                after=4;   %gives me 10min before, 2min of stim, and the 8min after 

                VV=cell(1,numstim);
                for s=1:numstim %iterating over all the stimulations
                step=stimep1(s); 
                eps=step-before*ne:step+after*ne-1; %returns the desired period
                if eps(1)<1 
                    continue; 
                end
                if eps(end)>21600 
                    continue; 
                end
    %             eps1=step-before*ne:step-1; %returns the desired period
    %             eps2=step+1:step+2*ne; %returns the desired period
    %             ws=wake(eps);  ms=move(eps); Ws=[Ws;ws+ms];  
                ws=wakm(eps); Ws=[Ws;ws]; ws2=reshape(ws,15,[]); ws2=mean(ws2); Ws2=[Ws2;ws2];
                ns=nrem(eps); Ns=[Ns;ns]; ns2=reshape(ns,15,[]); ns2=mean(ns2); Ns2=[Ns2;ns2]; 
                rs=rems(eps); Rs=[Rs;rs]; rs2=reshape(rs,15,[]); rs2=mean(rs2); Rs2=[Rs2;rs2];
                vv=Vmat(:,eps); 
                vv2=reshape(vv,3,15,[]); vv2=nanmean(vv2,2); vv2=squeeze(vv2);
                VV(1,s)={vv2};
                end

                mWmat=mean(Ws2(:,3)); 
                mNmat=mean(Ns2(:,3)); 
                mRmat=mean(Rs2(:,3)); 
                VVmeanmat=nanmean(cat(3,VV{:}),3);
                VVmean(1,anim)={VVmeanmat};
    %         
    %             w_boutDurs = 4.*(w_boutOffsets-w_boutOnsets); %./60
    %             n_boutDurs = 4.*(n_boutOffsets-n_boutOnsets); %./60;
    %             r_boutDurs = 4.*(r_boutOffsets-r_boutOnsets); %./60;
        clear wake nrem rems move wakm;
        
        mW=[mW mWmat];  mN=[mN mNmat]; mR=[mR mRmat];
        end
        
        mmW=mean(mW); sdW=std(mW)/sqrt(anim);
        mmN=mean(mN); sdN=std(mN)/sqrt(anim);
        mmR=mean(mR); sdR=std(mR)/sqrt(anim);
        
        VVplot=mean(cat(3,VVmean{:}),3); 
        VVplot0=cat(3,VVmean{:}); 

         if dd==1; 
             geno;
             mW1=mW;  mN1=mN; mR1=mR;
             baseSD1=[mmW sdW; mmN sdN; mmR sdR];
         elseif dd==2        
             mW2=mW;  mN2=mN; mR2=mR;
             baseSD2=[mmW sdW; mmN sdN; mmR sdR];
         end
    end

    geno
    baseSD1=baseSD1*100
    baseSD2=baseSD2*100
    [hW,pW]=ttest(mW1,mW2) % [hW,pW,ciW,statsW]
    [hN,pN]=ttest(mN1,mN2)
    [hR,pR]=ttest(mR1,mR2)
end
    

