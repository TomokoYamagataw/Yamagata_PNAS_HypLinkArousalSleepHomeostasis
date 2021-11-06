
clear all
close all

path='E:\Optoinh\'; 
pathFigs = [path,'FiguresWakeSuppression\']; 



%%%%% GFP controls 
mousenames1=[1 2 4 5]; 
days1=['120421 080421';'120421 080421';'060521 100521';'100521 060521'];

%%%%% Arch
mousenames2=[1 2 4 5 8 9];
days2=['080421 120421';'120421 080421';'080421 120421';'100521 060521';'100521 060521';'060521 100521'];


tfstim=1;
tnstim=22; 

ne=15; %there are 15 epochs in 1min 
before=2;  %we want to analyze the 2min before stimulus onset and 2min after
after=3;   %gives me 2min before, 5min of stim, and the 2min after 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defining all the derivations
ders=strvcat('fro','occ');
der=1; 
deri=ders(der,:);

maxep=10800; 
epochs=1:maxep; 
epochs=epochs./900;
fs=256; 
epochl=4; %states how long the epoch is
zermat=zeros(1,maxep); %creates an array of zeros with 1 row and as long as the total number of epoch

x=-before*ne:1:after*ne; %creates an array of epoch indices in the 2min before to after period
x=x*epochl; %converts the epoch indices to seconds 

pathvs=[path,'outputVS\']; %defines the folder with output files of vigilance states for mice2 (ChR2 exp)
pathstim=[path,'STIMs\']; %gives the folder with stimulation inf

vsname=strvcat('Wake','NREM','REM'); %makes a char array with rows of the passed strings of vigilance states



for geno = 1:2
    if geno==1 %geno=1 ctrl group
        mousenames=mousenames1;days=days1; pathvs=pathvs; gn='Gf'; %firststims=firststims1;
    elseif geno==2
        mousenames=mousenames2; days=days2; pathvs=pathvs; gn='Ar';  %firststims=firststims2;
    end
    numanim=length(mousenames); %how many animals were recorded
    [SpecGramFro,SpecGramOcc] = deal(cell(numanim,2));

    for anim=1:numanim %go through this loop as many times as there are animals
        mousename=[gn,num2str(mousenames(anim))]; mousename(isspace(mousename))=[];
        daysi=days(anim,:); is=find(isspace(daysi));
        for dd=[2,1]

            if dd==1 day=daysi(1:is(1)-1); else day=daysi(is(1):end); end
            day(isspace(day))=[];

            if dd==2
            fn0=[mousename,'_',day,'_stim']; %makes the full file name mouse_date_condition
            eval(['load ',pathstim,fn0,'.mat startend ']);

            stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
            stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2) 

%             outs=find((stimep2-stimep1)<70);stimep1(outs)=[];stimep2(outs)=[];
%                 if DUR==1
%                     stimep1(stimep2>=10800)=[]; stimep2(stimep2>=10800)=[];
%                 elseif DUR==2
%                     stimep1(stimep2<10800)=[]; stimep2(stimep2<10800)=[];
%                 end
            end


            fn=[mousename,'-',day,'-',deri]; %makes the full file name for the vigilance state and derivation output desired
            eval(['load ',pathvs,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);%
            currFro = spectr;
            VS=zeros(7,maxep); 
            VS(1,w)=1; VS(2,w1)=1; VS(3,nr)=1; VS(4,nr2)=1; VS(5,r)=1; VS(6,r3)=1; VS(7,mt)=1;
            clear w nr r w1 nr2 r3 mt ma bastend ; %clears the variables for the next loop 
            
            %%%% wake only
            art=sum(VS([2 3 4 5 6 7],:));  
            %%%% wnr all without artifact
            %%%% art=sum(VS([2 4 6 7],:)); 

            fn1=[mousename,'-',day,'-',ders(2,:)];
            currOcc= load(fullfile(pathvs,[fn1,'.mat']));
            currOcc = currOcc.spectr;
            
            fn2=[mousename,'-',day,'-epochSO'];
            epochs = load(fullfile(pathvs,[fn2,'.mat']));
            epochSO=epochs.epochSO;
            epochSDon=epochs.epochSDon;
            epochSDoff=epochs.epochSDoff;
            
            matsize=epochSDoff-epochSDon;
            eps=epochSDon:epochSDoff-1;
            
            currFrona=currFro; currOccna=currOcc;
            art=logical(art);
            currFrona(art,:)=NaN; currOccna(art,:)=NaN;
            spectrmFro=currFrona(eps,:);  SpectrmOcc=currOccna(eps,:); 
%             SPsFro(anim,:)=nanmean(spectrmFro,1); 
%             SPsOcc(anim,:)=nanmean(SpectrmOcc,1);

        SpecGramFro{anim,dd} = nanmean(spectrmFro,1)';
        SpecGramOcc{anim,dd} = nanmean(SpectrmOcc,1)';
        end
    end
    SpecGramFroG{geno}=SpecGramFro;
    SpecGramOccG{geno}=SpecGramOcc;
end
%%
   
freq = 0.5:0.25:30;
thetaRange = and(freq>=6,freq<=9);
SWARange = and(freq>0.5,freq<4); 
[froSWA,occSWA,froTheta,occTheta,genoID]=deal(cell(3,1));


%% plot average spectrogram normalised to baseline 
figure;

froBlGFP = SpecGramFroG{1}(:,1);
froBlGFP = cat(2,froBlGFP{:});

froStimGFP = SpecGramFroG{1}(:,2);
froStimGFP = cat(2,froStimGFP{:});
froStimGFP = 100*(froStimGFP./froBlGFP);

froBlArch = SpecGramFroG{2}(:,1);
froBlArch = cat(2,froBlArch{:});

froStimArch = SpecGramFroG{2}(:,2);
froStimArch = cat(2,froStimArch{:});
froStimArch = 100*(froStimArch./froBlArch);


aPrismSpectraGFPArchFro=[froStimGFP froStimArch];

occBlGFP = SpecGramOccG{1}(:,1);
occBlGFP = cat(2,occBlGFP{:});

occStimGFP = SpecGramOccG{1}(:,2);
occStimGFP = cat(2,occStimGFP{:});
occStimGFP = 100*(occStimGFP./occBlGFP);

occBlArch = SpecGramOccG{2}(:,1);
occBlArch = cat(2,occBlArch{:});

occStimArch = SpecGramOccG{2}(:,2);
occStimArch = cat(2,occStimArch{:});
occStimArch = 100*(occStimArch./occBlArch);

aPrismSpectraGFPArchOcc=[occStimGFP occStimArch];




