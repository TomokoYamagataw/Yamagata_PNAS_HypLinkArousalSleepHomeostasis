close all

%% v3: compare spontaneous arousals from NR @baseline day to stimulation induced arousals from NR @stimulation day
clear all
% pathToHelpfulFunctions = 'C:\Users\Martin Kahn\Desktop\Scripts_and_Code\Matlab';
% addpath(genpath(pathToHelpfulFunctions));2

path='E:\Optoinh\';
pathStim ='E:\Optoinh\STIMs';
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 230 160 0;240 228 64;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue, orange, yellow, black/white
npgCMap = [npgCMap;npgCMap];
pathFigs = [path,'FiguresWakeSuppression\']; 

%%%%% GFP controls 
mousenames1=[1 2 4 5]; 
days1=['070421 150421';'070421 150421';'050521 150521';'050521 140521'];

%%%%% Arch
mousenames2=[1 2 4 5 8 9];
days2=['070421 160421';'070421 150421';'070421 160421';'050521 140521';'050521 140521';'050521 150521'];


ders=strvcat('fro','occ');der=1;deri=ders(der,:);
groups = {'GFP','Arch'};
maxep=21600;x=1:maxep; x=x./900;epochl=4;
zermat=zeros(1,maxep);

pathvs=[path,'outputVS\'];

vsname=strvcat('Wake','NREM','REM','SWA')
int=2;
numint=24/int;
numh=900; % num epochs in 1h
x=int/2:int:24;
fieldNames = {'wake','nrem','rem','swa','sleep_cumsum','w_cumsum','nr_cumsum','r_cumsum','swa_cumsum','mt_cumsum'};
rez = cell(numel(fieldNames),numel(groups));  % results container
rez=cell2struct(rez,fieldNames);
SWAs=[];
[sandwiches,sandwichTc] = deal({});
[wSpectraFro,wSpectraOcc,wSpecGramFro,wSpecGramOcc] = deal(cell(numel(groups),1));
for geno=1:numel(groups)
    if geno==1
        mousenames=mousenames1;days=days1;gn='Gf';
    elseif geno ==2
        mousenames=mousenames2;days=days2;gn='Ar';

    end
    numanim=length(mousenames);
    [wakeSpectraFro,wakeSpectraOcc,wDurs,wOnsets,wakeSpecGramOcc,wakeSpecGramFro] = deal(cell(numanim,2));
    for anim=1:numanim
        mousename=[gn,num2str(mousenames(anim))];
        daysi=days(anim,:); is=find(isspace(daysi));
        
        [SWA1,W1,N1,R1,SWA1_cs,W1_cs,mt_cs,N1_cs,R1_cs,SWA_cs]=deal([]);
        
        
        wnr=[];
        
        for dd=[2,1]
            
            if dd==1 day=daysi(1:is(1)-1); else day=daysi(is(1):end); end
            day(isspace(day))=[];
            
            % load recording
            fn1=[mousename,'-',day,'-',deri];
            eval(['load ',pathvs,fn1,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);
            currFro = spectr;
            
            
            % get vigilance states
            VS=zeros(7,maxep);
            VS(1,w)=1;VS(2,w1)=1;VS(3,nr)=1;VS(4,nr2)=1;VS(5,r)=1;VS(6,r3)=1;VS(7,mt)=1;
            clear w nr r w1 nr2 r3 mt ma bastend;
            
         
%%%%%%%%%%%%            
% VS(:,10801:end)=[];            
            
            wake=sum(VS(1:2,:));nrem=sum(VS(3:4,:));rems=sum(VS(5:6,:));mt  = VS(7,:);
            
            w_boutOnsets = find(diff([0,wake])==1);
            w_boutOffsets = find(diff([wake,0])==-1); % last epoch that still contains waking
            w_boutDurs = 4.*(w_boutOffsets-w_boutOnsets)./60;
            
            % find waking bouts (that start from NREM and happen in light period)
            minDurInEps = 25;
            wEpochs = find(wake)';
            wBoutStartEp= wEpochs([1;find(diff(wEpochs)~=1)+1]); % find wake bout onsets
            wBoutEndEp=wEpochs([find(diff(wEpochs)~=1);length(wEpochs)]);
            nremEps = find(nrem);
            nrToWOnly=arrayfun(@(x) sum(nremEps==(x-1))~=0,wBoutStartEp);
            lightPeriodIdxs = wBoutEndEp<(12*3600./4);% only light period
            boutsToUse = and(nrToWOnly,lightPeriodIdxs);
            wBoutStartEp=wBoutStartEp(boutsToUse);
            wBoutEndEp=wBoutEndEp(boutsToUse);
            wBoutEps = arrayfun(@(x) wBoutStartEp(x):wBoutEndEp(x),1:length(wBoutStartEp),'un',0);
            
            % if stim day then get episodes where stim starts/ends
            if dd ==2
                fn1=[mousename,'_',day,'_','stim'];
                tmpStim= load(fullfile(pathStim,[fn1,'.mat']));
                stimOffsetEp = tmpStim.startend;
                stimOffsetEp=stimOffsetEp(stimOffsetEp(:,2)<256*(12*3600),:);
                stimOnsetEp = floor(stimOffsetEp(:,1)/(256*4));
                
                stimOffsetEp = floor(stimOffsetEp(:,2)/(256*4));
                allStims = cell2mat(arrayfun(@(x)stimOnsetEp(x):stimOffsetEp(x),1:length(stimOffsetEp),'un',0));
                wBoutStimAtStart = cellfun(@(x) ~isempty(intersect(x(1),allStims)),wBoutEps);
            else % if it's baseline, select same amount of wake bouts and try match first the timing and then the duration
                wBoutStimAtStart = ones(length(wBoutStartEp),1);
            end
            wBoutStimAtStart = logical(wBoutStimAtStart);
            wBoutStartEp=wBoutStartEp(wBoutStimAtStart);
            wBoutEndEp=wBoutEndEp(wBoutStimAtStart);
            wBoutDur = (wBoutEndEp - wBoutStartEp+1)*4./60;
            wBoutEps= wBoutEps(wBoutStimAtStart);
            
            
            % load occipital
            fn1=[mousename,'-',day,'-',ders(2,:)];
            currOcc= load(fullfile(pathvs,[fn1,'.mat']));
            currOcc = currOcc.spectr;
            % remove artefacts
            artFreeWake = find(VS(1,:));
            artFreewBoutEps=cellfun(@(x) intersect(artFreeWake,x),wBoutEps,'un',0); 
            
            wakeSpectraFro{anim,dd} = cell2mat(cellfun(@(x) nanmean(currFro(x,:),1),artFreewBoutEps,'un',0)');
            wakeSpectraOcc{anim,dd} = cell2mat(cellfun(@(x) nanmean(currOcc(x,:),1),artFreewBoutEps,'un',0)');
            wDurs{anim,dd} = wBoutDur;
            wOnsets{anim,dd} = wBoutStartEp;
            
            % make spectrograms 
            longEnough = wBoutDur>2;  
            CC=artFreewBoutEps(longEnough);
            if size(CC,2)>1
                CC=CC(1,1);
                CC1=cell2mat(CC);
    %             lognEnough2 = artFreewBoutEps(longEnough)> 30;
                if min(CC1)<=30
                    longEnough(1,1)=0;
                end
                tmp = cellfun(@(x) currFro(x(1)-30:x(1)+29,:),artFreewBoutEps(longEnough),'un',0);
                tmp = nanmean(cat(3,tmp{:}),3);
                wakeSpecGramFro{anim,dd} =tmp;
                tmp = cellfun(@(x) currOcc(x(1)-30:x(1)+29,:),artFreewBoutEps(longEnough),'un',0);
                tmp = nanmean(cat(3,tmp{:}),3);                
            else
                tmp=NaN*ones(60,121);
            end
            wakeSpecGramOcc{anim,dd} =tmp;

%             pcolor(-120:4:116,freq,tmp'); shading interp

            
            % time course
            te=sum(VS,1); 
            te=sum(reshape(te,numh*int,numint));
            
            
            % SWA
            swa=nanmean(currFro(:,3:17),2);swa=swa'; swaN=swa; swaN(VS(3,:)==0)=NaN;
            swa_sum =nansum(reshape(swaN,numh*int,numint));
            swa=nanmean(reshape(swaN,numh*int,numint));
            
            w=sum(reshape(wake,numh*int,numint));
            n=sum(reshape(nrem,numh*int,numint));
            r=sum(reshape(rems,numh*int,numint));
            mt = sum(reshape(mt,numh*int,numint));
            w = mt+w;
            
            W1_cs=[W1_cs;cumsum(w)];N1_cs=[N1_cs;cumsum(n)];R1_cs=[R1_cs;cumsum(r)]; SWA_cs = [SWA_cs;cumsum(swa_sum)];
            mt_cs = [mt_cs;cumsum(mt)];
            
            te_cs = cumsum(te);
            % normalise to duration of epoch
            w=w./te*100;
            n=n./te*100;
            r=r./te*100;
            W1=[W1;w];N1=[N1;n];R1=[R1;r];
            
            
            SWA1=[SWA1 swa];
            
        end
        
        
        
%         % equalise number and timings of wake bouts to use (try to select similar durs and ZT times)
%         [minNr,minIdx] = min([length(wDurs{anim,1}),length(wDurs{anim,2})]);
%         maxIdx = find((1:2)~=minIdx);
%         longerDurs = wDurs{anim,maxIdx};
%         shorterDurs = wDurs{anim,minIdx};
%         longerOnsets = wOnsets{anim,maxIdx};
%         shorterOnsets = wOnsets{anim,minIdx};
%         newLongerVector = nan(size(shorterDurs));
%         longerVectorOriginalIdxs = 1:numel(longerOnsets);
%         for currEp= 1 :numel(shorterOnsets)
%             similarDur = find(abs(shorterDurs(currEp) - longerDurs)<2); % try to find similar durations
%             if isempty(similarDur)
%                 [~,bestMatchIdx] = min(abs(shorterDurs(currEp) - longerDurs));
%             else
%                 [~,bestMatchIdx] = min(abs(shorterOnsets(currEp) - longerOnsets(similarDur)));
%                 bestMatchIdx = similarDur(bestMatchIdx);
%             end
%             newLongerVector(currEp)=longerVectorOriginalIdxs(bestMatchIdx); % not down which idx of the original vector we just selected
%             longerDurs(bestMatchIdx)=[]; % make sure the same index doesn't get selected twice
%             longerOnsets(bestMatchIdx)=[];
%             longerVectorOriginalIdxs(bestMatchIdx)=[];
%         end
%         % now change th longer vector spectra
%         wakeSpectraOcc{anim,maxIdx}= wakeSpectraOcc{anim,maxIdx}(newLongerVector,:);
%         wakeSpectraFro{anim,maxIdx}= wakeSpectraFro{anim,maxIdx}(newLongerVector,:);
%         
%         
        
        rez(geno).wake(anim,:,:)  =W1;
        rez(geno).nrem(anim,:,:)  =N1;
        rez(geno).rem(anim,:,:)  =R1;
        rez(geno).swa(anim,:,:)  =SWA1;
        
        
        rez(geno).w_cumsum(anim,:,:)=W1_cs.*(4./60);
        rez(geno).w_cumsum(anim,:,:)=W1_cs.*(4./60);
        
        
        rez(geno).nr_cumsum(anim,:,:)=N1_cs.*(4./60);
        rez(geno).r_cumsum(anim,:,:)=R1_cs.*(4./60);
        rez(geno).swa_cumsum(anim,:,:)=SWA_cs;
        rez(geno).mt_cumsum(anim,:,:)=mt_cs.*(4./60);
        rez(geno).sleep_cumsum(anim,:,:)=(N1_cs+ R1_cs).*(4./60);
        
    end
    

    
    aPrism_meanRez(geno)=structfun(@(x)nanmean(x,1),rez(geno),'un',0);
    aPrism_SEMRez(geno)=structfun(@(x)nanstd(x,[],1)./sqrt(anim),rez(geno),'un',0);
    wSpectraFro{geno} = wakeSpectraFro;
    wSpectraOcc{geno} = wakeSpectraOcc;
        wSpecGramFro{geno} = wakeSpecGramFro;
        wSpecGramOcc{geno} = wakeSpecGramOcc;
    
end

%% total duration of NREM during the light period (cause that's what's relevant here)
vsToPlot = {'sleep_cumsum','nr_cumsum','r_cumsum','r_proportion'};
labelling = {'total sleep','NREM','REM','proportion REM'};
[allStats,groups2,etaSqu,LME,LLR] = deal([]);
[postHocTests,rmANOVA] = deal({});
clrArch = [230 160 0]./255;  %define colors according to tomookos plot
clrGFP = [160,160,160]./255;
colorsToUse = {clrGFP,clrArch};

figure; hold on;
for vs = 1:4
    [stats,groups,b,statsRMANOVA] = deal([]);
    for geno = 1:2
        if vs == 4
            currRem = squeeze(rez(geno).r_cumsum(:,:,6));
            currNr = squeeze(rez(geno).nr_cumsum(:,:,6));
            curr = 100.*currRem./(currRem+currNr);
        elseif vs == 5
            currNr = squeeze(rez(geno).nr_cumsum(:,:,6));
            currR = squeeze(rez(geno).r_cumsum(:,:,6));
            currMt = squeeze(rez(geno).mt_cumsum(:,:,6));
            curr = currNr+currR+currMt;           
        else 
            curr = rez(geno).(vsToPlot{vs});
            curr = squeeze(curr(:,:,6));
        end
        relCurr = 100.*curr(:,2)./curr(:,1);
        means = nanmean(relCurr,1); SEM = std(relCurr)./sqrt(numel(relCurr));
        
        [~,postHoc] = ttest(curr(:,1),curr(:,2));
        statsRMANOVA=[statsRMANOVA;curr];
        stats = [stats;postHoc];
        groups= [groups;ones(size(relCurr)).*geno];
        b{geno}=bar(geno+(vs-1)*5,means,'FaceColor',colorsToUse{geno});
        errorbar(geno+(vs-1)*5,means,SEM,'k');
        display(['mean (SEM) ',vsToPlot{vs},num2str(means),'(',num2str(SEM),')'])
    end
    aPrism_rmANOVA{vs} =  statsRMANOVA';
    statsTable = struct;
    statsTable.timeInState = statsRMANOVA(:);
    statsTable.group = repmat(groups,2,1);
    statsTable.day = ones(size(statsRMANOVA)).*(1:2);
    statsTable.day = statsTable.day(:);
    statsTable.mouseID = repmat(1:length(statsRMANOVA),1,2)';
    statsTable = struct2table(statsTable);
    statsTable.group = nominal(statsTable.group);
    statsTable.day = nominal(statsTable.day);
    statsTable.mouseID = nominal(statsTable.mouseID);
%     aPrism_statsTable=statsTable';]
    LMEs{vs} = fitlme(statsTable,'timeInState ~ 1+day*group+(1|mouseID)','Fitmethod','ML');
    LMEcontrol = fitlme(statsTable,'timeInState ~ 1+day+group+(1|mouseID)','Fitmethod','ML');
    LLR{vs} = compare(LMEcontrol,LMEs{vs});
    legend([b{1} b{2}],{'GFP','Arch'},'autoupdate','off'); legend('boxoff');
    postHocTests{vs}=stats;
end
line(get(gca,'XLim'),[100 100],'LineStyle','--','Color','k');
ax = gca; ax.XTick = 2:5:17; ax.XTickLabel = labelling;
ylabel('time in state [% baseline]');
%% plot spectra
% wSpectraOcc(1)=[];wSpectraFro(1)=[];
%% plot spectra
clrArch = [230 160 0]./255;  %define colors according to tomokos plot
clrGFP = 'k';
colorsToUse = {clrGFP,clrArch};
plots = cell(2,2);
figure;
[pwStatsFro,pwStatsOcc,RMStatsFro,RMStatsOcc] = deal([]);
statsTable = struct('Fropower',[],'Occpower',[],'mouseID',[],'freqID',[],'genoType',[]);
freq = 0:0.25:30;
thetaRange = and(freq>=6,freq<=9);
SWARange = and(freq>0.5,freq<4); 
[froSWA,occSWA,froTheta,occTheta,genoID]=deal(cell(2,1));
for geno = 1:numel(wSpectraOcc)
    nonEmptyCells = ~cellfun(@isempty,wSpectraFro{geno});
    
    % stim day 
    meanFro_stim = cell2mat(cellfun(@(x) nanmean(x,1),wSpectraFro{geno}(nonEmptyCells(:,2),2),'un',0));
    meanOcc_stim = cell2mat(cellfun(@(x) nanmean(x,1),wSpectraOcc{geno}(nonEmptyCells(:,2),2),'un',0));
    
    % baseline day
    meanFro_bl = cell2mat(cellfun(@(x) nanmean(x,1),wSpectraFro{geno}(nonEmptyCells(:,1),1),'un',0));
    meanOcc_bl = cell2mat(cellfun(@(x) nanmean(x,1),wSpectraOcc{geno}(nonEmptyCells(:,1),1),'un',0));
    
    froSWA{geno} = [nanmean(meanFro_bl(:,SWARange),2),nanmean(meanFro_stim(:,SWARange),2)];
    occSWA{geno}  = [nanmean(meanOcc_bl(:,SWARange),2),nanmean(meanOcc_stim(:,SWARange),2)];
    froTheta{geno}  = [nanmean(meanFro_bl(:,thetaRange),2),nanmean(meanFro_stim(:,thetaRange),2)];
    occTheta{geno}  = [nanmean(meanOcc_bl(:,thetaRange),2),nanmean(meanOcc_stim(:,thetaRange),2)];
    genoID{geno}  = geno.*ones(size(froSWA{geno}));
    
%     % prepare for mixed anova
%     tmp_stim = meanOcc_stim; 
%     tmp_bl = meanOcc_bl; 
%     tmp_stim(:,end+1)= geno; tmp_stim(:,end+1)= 2;
%     tmp_bl(:,end+1)= geno;tmp_stim(:,end+1)= 1;
%     RMStatsOcc = cat(1,RMStatsOcc,[tmp_bl;tmp_stim]);
     
    % plot relative values 
    meanFro_ratio = 10.*log10(meanFro_stim./meanFro_bl);
    meanOcc_ratio = 10.*log10(meanOcc_stim./meanOcc_bl);
    
    meanFro_ratio = 100.*(meanFro_stim./meanFro_bl);
    meanOcc_ratio = 100.*(meanOcc_stim./meanOcc_bl);
    
    % run multiple pairwise comparisons 
    if geno == 2
    pwStatsFro = meanFro_ratio; 
    pwStatsOcc = meanOcc_ratio; 
    sidaksAlpha=1-(1-0.05).^(1/121); 

    else
        combineGenoStatsFro = cat(1,meanFro_ratio,pwStatsFro);
        combineGenoStatsOcc = cat(1,meanOcc_ratio,pwStatsOcc);
        [combineGenoStatsFro,h] = arrayfun(@(x) signtest(combineGenoStatsFro(:,x),100,'alpha',sidaksAlpha),1:size(pwStatsOcc,2)); % test against 0 
        combineGenoStatsOcc = arrayfun(@(x) signtest(combineGenoStatsOcc(:,x),100,'alpha',sidaksAlpha),1:size(pwStatsOcc,2));
        
        pwStatsFro=arrayfun(@(x) ranksum(pwStatsFro(:,x)',meanFro_ratio(:,x)')<sidaksAlpha,1:size(pwStatsOcc,2));
        pwStatsOcc=arrayfun(@(x) ranksum(pwStatsOcc(:,x)',meanOcc_ratio(:,x)')<sidaksAlpha,1:size(pwStatsOcc,2));
    end


    % prepare LME 
    statsTable.Fropower = [statsTable.Fropower;meanFro_ratio(:)];
    statsTable.Occpower = [statsTable.Occpower;meanOcc_ratio(:)];
    mouseIDs= repmat((1:size(meanFro_ratio,1))',1,121)+(geno-1).*size(froSWA{geno},1);
    statsTable.mouseID = [statsTable.mouseID;mouseIDs(:)];
    %         dayID= ones(size(currFro)).*dd;
    %         statsTable.dayID = [statsTable.dayID;dayID(:)];
    genoType= ones(size(meanFro_ratio)).*geno;
    statsTable.genoType = [statsTable.genoType;genoType(:)];
    freqID= ones(size(meanFro_ratio)).*(1:121);
    statsTable.freqID = [statsTable.freqID;freqID(:)];
    
    
    
    subplot(2,1,1)
    plots{1,geno}=patchMeUp(freq,nanmean(meanFro_ratio),std(meanFro_ratio)./sqrt(size(meanFro_ratio,1)),colorsToUse{geno},0.3); hold on;
%     p1=plot(freq,meanFro_ratio,'Color',colorsToUse{geno}); hold on;for i = 1:numel(p1); p1(i).Color(4)=0.3;end
%     if geno==2
%         bar(freq(pwStatsFro),80.*ones(sum(pwStatsFro),1),'k');
%         bar(freq(combineGenoStatsFro<0.05),50.*ones(sum(combineGenoStatsFro<0.05),1),'Facecolor','r','Edgecolor','none');
%     end
    box off; xlabel('frequency [Hz]'); ylabel('power relative to baseline [%]');title('frontal');
    xlim([0.5,25]); line(get(gca,'XLim'),[100 100],'Color','k','LineStyle','--'); set(gca,'TickDir','out')

    subplot(2,1,2)
    plots{2,geno}=patchMeUp(freq,nanmean(meanOcc_ratio),std(meanOcc_ratio)./sqrt(size(meanOcc_ratio,1)),colorsToUse{geno},0.3); hold on;
%     p1=plot(freq,meanOcc_ratio,'Color',colorsToUse{geno}); hold on;for i = 1:numel(p1); p1(i).Color(4)=0.3;end
%     if geno==2
%         bar(freq(pwStatsOcc),80.*ones(sum(pwStatsOcc),1),'k');
%         bar(freq(combineGenoStatsOcc<0.05),50.*ones(sum(combineGenoStatsOcc<0.05),1),'Facecolor',colorsToUse{geno},'Edgecolor','none');
%     end

    box off; xlabel('frequency [Hz]'); ylabel('power relative to baseline [%]');title('occipital');
    xlim([0.5,25]);line(get(gca,'XLim'),[100 100 ],'Color','k','LineStyle','--'); set(gca,'TickDir','out');
   
end
subplot(2,1,1)
legend([plots{1,1},plots{1,2}],{'GFP','Arch'}); legend('boxoff'); linkaxes(); 
ylim([41.3771  173.1703]);

saveas(gcf,fullfile(pathFigs,'spontVsOptoArousalSuppression_spectra.fig'));
saveas(gcf,fullfile(pathFigs,'spontVsOptoArousalSuppression_spectra.png'));
saveas(gcf,fullfile(pathFigs,'spontVsOptoArousalSuppression_spectra.svg'));
%%
% LME 
statsTable = struct2table(statsTable); 
statsTable.Fropower =statsTable.Fropower; 
statsTable.Occpower =statsTable.Occpower; 
statsTable.mouseID=nominal(statsTable.mouseID);
artefactFreqs = or(and(freq>=9,freq<=11),and(freq>=19,freq<=21)); 
artefactFreqs = find(or(artefactFreqs,and(freq>=29,freq<=31)));
freqsToExclude=cell2mat(arrayfun(@(x)find((statsTable.freqID)==x),artefactFreqs,'un',0)');
statsTable(freqsToExclude,:)=[];
statsTable.freqID=nominal(statsTable.freqID); 
statsTable.genoType=nominal(statsTable.genoType); 
fullLME = fitlme(statsTable,'Occpower~1+freqID*genoType+(1|mouseID)','FitMEthod','ML');
% nullModel = fitlme(statsTable,'Fropower~1+genoType+(1|mouseID)','FitMEthod','ML');
fullLME.anova


% frequency-band-wise interactions
forSPSSFroSWA = cat(1,froSWA{:});
forSPSSOccSWA= cat(1,occSWA{:});
nMiceGFP = size(froSWA{1},1);
nMiceArch = size(froSWA{2},1);

forSPSSFroTheta = cat(1,froTheta{:});
forSPSSOccTheta= cat(1,occTheta{:});
forSPSSgenoID= cat(1,genoID{:});
toPlotALL = cat(3,forSPSSFroSWA,forSPSSOccSWA,forSPSSFroTheta,forSPSSOccTheta); 
toPlotALL = 100.*squeeze(toPlotALL(:,2,:)./toPlotALL(:,1,:));
toPlotArch =toPlotALL(nMiceGFP+1:end,:); 
toPlotGFP = toPlotALL(1:nMiceGFP,:); 
groupNames = {'Fro-SWA','Occ-SWA','Fro-Theta','Occ-Theta'};

figure; 
for pltCnt = 1:size(toPlotGFP,2)
    bar(pltCnt*4,nanmean(toPlotGFP(:,pltCnt),1),'Facecolor',colorsToUse{1},'FaceAlpha',0.3); hold on;
    errorbar(pltCnt*4,nanmean(toPlotGFP(:,pltCnt),1),std(toPlotGFP(:,pltCnt),[],1)./sqrt(nMiceGFP),'-','Color','k')
    
    bar(pltCnt*4+1,nanmean(toPlotArch(:,pltCnt),1),'Facecolor',colorsToUse{2}); hold on;
    errorbar(pltCnt*4+1,nanmean(toPlotArch(:,pltCnt),1),std(toPlotArch(:,pltCnt),[],1)./sqrt(nMiceArch),'-','Color','k')
    
end
ax = gca; ax.TickDir='out'; box off; ax.XTick = 4.5:4:pltCnt*4+0.5; ax.XTickLabel = groupNames;
ylabel('power (% baseline)');line(get(gca,'XLim'),[100 100],'Color','k','LineStyle','--');
%% try correlate increase in theta during W with increase in SWA during sleep 
SWALightStim = cat(1,nanmean(rez(1).swa(:,:,1:6),3),mean(rez(2).swa(:,:,1:6),3)); 
SWALightBl = cat(1,nanmean(rez(1).swa(:,:,13:18),3),mean(rez(2).swa(:,:,13:18),3)); 
diffNRSWA = SWALightStim-SWALightBl;
diffWThetaOcc = forSPSSOccTheta(:,2) - forSPSSOccTheta(:,1);
diffWSWAFro = forSPSSFroSWA(:,2) - forSPSSFroSWA(:,1);
diffWThetaFro = forSPSSFroTheta(:,2) - forSPSSFroTheta(:,1);
diffWSWAOcc = forSPSSOccSWA(:,2) - forSPSSOccSWA(:,1);

toCorrelate = {diffWThetaOcc,diffWSWAFro,diffWThetaFro,diffWSWAOcc};
xLabelling = {'\Delta wake-thetaOcc','\Delta wake-swaFro','\Delta wake-thetaFro','\Delta wake-swaOcc'};
figure;
for pltCnt = 1:numel(toCorrelate)
    subplot(2,2,pltCnt)
    scatter(toCorrelate{pltCnt}(1:nMiceGFP),diffNRSWA(1:nMiceGFP),15,'MarkerFaceColor',colorsToUse{1},'MarkerEdgeColor','none');hold on;
    scatter(toCorrelate{pltCnt}(nMiceGFP+1:end),diffNRSWA(nMiceGFP+1:end),15,'MarkerFaceColor',colorsToUse{2},'MarkerEdgeColor','none');hold on;
    currCorr=corrcoef(toCorrelate{pltCnt},diffNRSWA); ylabel('\Delta NREM SWA');
    currSlope = polyfit(toCorrelate{pltCnt},diffNRSWA,1);
    xlims = get(gca,'XLim');
    regLine = currSlope(2) + (xlims(1):xlims(2)).*currSlope(1);
    plot(xlims,regLine([1,end]),'-k')
    xlabel(xLabelling{pltCnt});
    title(['R=',num2str(currCorr(2))]);
    xlim(get(gca,'XLim')+([-50 50])); 
    ylim(get(gca,'YLim')+([-50 50])); 
    set(gca,'TickDir','out','FontName','Arial','FontSize',9)
end
saveas(gcf,fullfile(pathFigs,'corrSWATheta.fig'));
saveas(gcf,fullfile(pathFigs,'corrSWATheta.png'));
saveas(gcf,fullfile(pathFigs,'corrSWATheta.svg'));


% stats: compare correlation coefficients (Andy Field pp 362), for dependent corrcoefs 
Rxy = corr([diffWThetaOcc,diffNRSWA]);
Rzy = corr([diffWSWAFro,diffNRSWA]);
Rxz = corr([diffWThetaOcc,diffWSWAFro]);
Rxy = Rxy(2); Rzy = Rzy(2); Rxz = Rxz(2); 
nMice = length(diffWThetaOcc);
tVal = (Rxy - Rzy).*sqrt((nMice-3).*(1+Rxy)./(2.*(1-Rxy^2-Rxz^2-Rzy^2+2*Rxy*Rxz*Rzy)));
tinv(0.01,nMice-3)



%% plot average spectrogram normalised to baseline 
occBlGFP=wSpecGramOcc{1}(:,1);
occBlGFP = cat(3,occBlGFP{:});

occStimGFP=wSpecGramOcc{1}(:,2);
occStimGFP = cat(3,occStimGFP{:});
occStimGFP = nanmean(10.*log10(occStimGFP./occBlGFP),3);

occBlArch=wSpecGramOcc{2}(:,1);
occBlArch = cat(3,occBlArch{:});

occStimArch=wSpecGramOcc{2}(:,2);
occStimArch = cat(3,occStimArch{:});
occStimArch = nanmean(10.*log10(occStimArch./occBlArch),3);


freq = 0.5:0.25:30;
timeInS = -116:4:120; %-116:4:120;
% f1=figure('name','fro');f2=figure('name','occ');


figure;

toPlot = {occStimGFP,occStimArch};
minMax = cell2mat(cellfun(@(x)x(:,1:32),toPlot,'un',0)');
 minMax = [min(minMax(:)) max(minMax(:))];
for plotCnt = 1:2
    subplot(2,2,plotCnt)
    curSpec =smoothdata(squeeze(toPlot{plotCnt}(:,3:121)),1,'movmedian',5);
    pcolor(timeInS,freq,curSpec'); shading flat; colorbar();
    xlabel('time to stim [s]'); ylabel('hz');colormap('jet'); caxis(minMax)
end


%% plot example spectrograms 
exampleMice = [4,5]; % GFP and Arch 

occBlGFP=wSpecGramOcc{1}(:,1);
occBlGFP = log10(cat(3,occBlGFP{:}));

occStimGFP=wSpecGramOcc{1}(:,2);
occStimGFP = log10(cat(3,occStimGFP{:}));

occBlArch=wSpecGramOcc{2}(:,1);
occBlArch = log10(cat(3,occBlArch{:}));

occStimArch=wSpecGramOcc{2}(:,2);
occStimArch = log10(cat(3,occStimArch{:}));

toPlot{1} = {occBlGFP,occStimGFP};
toPlot{2} = {occBlArch,occStimArch};

freq = 0.5:0.25:30;
timeInS = -116:4:120;
genos = {'GFP','Arch'};
titleStrings = {'baseline','stim'};
% f1=figure('name','fro');f2=figure('name','occ');
for geno = 1:2
    currToPlot = toPlot{geno};
    for mouseCnt = 1:size(currToPlot{1},3)
        figure;


        minMax = cell2mat(cellfun(@(x)x(:,:,mouseCnt),currToPlot,'un',0)');
        minMax = [min(minMax(:)) max(minMax(:))];
        for plotCnt = 1:2
            subplot(1,2,plotCnt)
            curSpec =smoothdata(squeeze(currToPlot{plotCnt}(:,3:121,mouseCnt)),1,'movmedian',5);
            
            pcolor(timeInS,freq,curSpec'); shading flat; colorbar();
            xlabel('time to stim [s]'); ylabel('hz');colormap('jet'); caxis(minMax)
            title([titleStrings{plotCnt},'-',genos{geno}])

        end
%         if mouseCnt == exampleMice(geno)
%             saveas(gcf,fullfile(pathFigs,['spVsOpto_SpecGram_',genos{geno},'.fig']));
%             saveas(gcf,fullfile(pathFigs,['spVsOpto_SpecGram_',genos{geno},'.png']));
%             saveas(gcf,fullfile(pathFigs,['spVsOpto_SpecGram_',genos{geno},'.svg']));
%         end
    end
    ylim([0 15])
end
