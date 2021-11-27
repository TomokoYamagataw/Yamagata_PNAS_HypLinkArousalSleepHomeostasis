%% compare arousal delay high and low sleep pressure conditions
clear all
addpath(genpath('C:\Users\Martin Kahn\Desktop\Scripts_and_Code\Matlab'))
%% load emg data and spectra etc 
paths.extractedSignals ='F:\Tomoko_OptoStim\EMGAnalysis_matfiles';
LPO=load(fullfile(paths.extractedSignals,'TenHz24h'));
nonLPO=load(fullfile(paths.extractedSignals,'TenHzNonLPO'));

% colormap 
npgCMap = [ 0 114 178; 0 158 115;  213 94 0;  204,121,167;86,180,233; 240 228 66;0 0 0]./255 ; %blue, green, vermillon, purple, skyblue,yellow
npgCMap = [npgCMap;npgCMap];
EMGThreshForFigures = 2; 
%% extract features: arousal delay, preceeding power in delta, theta, sigma, time to last waking

%settings for detecting arousals
EMGThresholds= [1,2,3,4,6,8]'; % EMG threshold
nEpochsSWA= ([1])'; % EMG threshold
LPO.settings.EMGThresholds = EMGThresholds;
LPO.settings.nEpochsSWA = nEpochsSWA;
nonLPO.settings.EMGThresholds = EMGThresholds;
nonLPO.settings.nEpochsSWA = nEpochsSWA;
LPO.settings.minDur = 4;
nonLPO.settings.minDur = 4;

% detect arousals based on thresholds as before
try nonLPO.settings=rmfield(nonLPO.settings,'forceThreshold'); LPO.settings=rmfield(LPO.settings,'forceThreshold'); catch; end
opto_LPO=extract_features_from_opto_evoked_arousals(LPO.EMGs,LPO.stims,LPO.scorings,LPO.EMGVarRatios,LPO.EMGVarBinning,LPO.settings); % runs the script
opto_nonLPO=extract_features_from_opto_evoked_arousals(nonLPO.EMGs,nonLPO.stims,nonLPO.scorings,nonLPO.EMGVarRatios,nonLPO.EMGVarBinning,nonLPO.settings); % runs the script


%% comapre NREM arousal delays LPO vs nonLPO 
f1=figure('name','ignoring mouse ID');
f2=figure('name','using mean of each mouse');
% f3=figure('name','distribution of arousal delay in each mouse');
tArousNR_nonLPO=opto_nonLPO.tArousal_NR(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
tArousNR_LPO=opto_LPO.tArousal_NR(1+((0:numel(EMGThresholds)-1)*numel(nEpochsSWA)),:); % we will ignore the entries for the varying epochs of SWA, as this doesn't change the tArousal
pValsTTest = nan(numel(EMGThresholds),1);
pValsLME = nan(numel(EMGThresholds),1);
nMice = size(tArousNR_nonLPO,2);

for EMGCnt = 1:numel(EMGThresholds) % cycle through various emg thresholds
    %extract the arousalT based on the current EMG threshold
    acrossMiceNR_nonLPO = cell2mat(tArousNR_nonLPO(EMGCnt,:));
    withinMiceNR_nonLPO = cellfun(@nanmean,tArousNR_nonLPO(EMGCnt,:));
    acrossMiceNR_LPO = cell2mat(tArousNR_LPO(EMGCnt,:));
    withinMiceNR_LPO = cellfun(@nanmean,tArousNR_LPO(EMGCnt,:));

    % plot difference LPO nonLPO across all samples, ignoring mouse ID
    figure(f1)
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    combined = [acrossMiceNR_nonLPO,acrossMiceNR_LPO]';
    groups = [ones(size(acrossMiceNR_nonLPO)),2.*ones(size(acrossMiceNR_LPO))]';
    boxplot(combined,groups,'notch','on');
    h = findobj(gca,'Tag','Box');
    h2= findobj(gca,'Tag','Outliers');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),npgCMap(j,:),'FaceAlpha',.4);
        h(j).Color = 'k'; % set box edge color
        h2(j).MarkerEdgeColor=[0.8 0.8 0.8];
    end
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]); 
    ylabel('arousal delay [s]'); ax = gca; ax.XTickLabel ={'nonLPO','LPO'};
    set(gca,'TickDir','out'); box off;set(gcf,'Position',[205 484 1264 415]);
    linkaxes

   % plot difference LPO nonLPO for each mouse sepparately
    figure(f2)
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
%     CIMedianLPO = 1.58*(iqr(withinMiceNR_LPO)./sqrt(numel(withinMiceNR_LPO)));
%     CIMedianNonLPO = 1.58*(iqr(withinMiceNR_nonLPO)./sqrt(numel(withinMiceNR_nonLPO)));
%     errorbar(1,median(withinMiceNR_nonLPO),CIMedianNonLPO,'ok','MarkerFaceColor','k'); hold on;
%     errorbar(2,median(withinMiceNR_LPO),CIMedianLPO,'ok','MarkerFaceColor','k')
    boxplot([withinMiceNR_LPO,withinMiceNR_nonLPO],[ones(size(withinMiceNR_LPO)),2.*ones(size(withinMiceNR_nonLPO))]);hold on;
    plot(2,withinMiceNR_nonLPO,'-ok','MarkerFaceColor','k','MarkerSize',4); hold on;
    plot(1,withinMiceNR_LPO,'-ok','MarkerFaceColor','k','MarkerSize',4);
    h = findobj(gca,'Tag','Box');
    h2= findobj(gca,'Tag','Outliers');
    h3= findobj(gca,'Tag','Median');

    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',.4,'EdgeColor','none');
        h(j).Color = 'none'; % set box edge color
        h2(j).MarkerEdgeColor='none'; h2(j).MarkerFaceColor='none';
        h3(j).LineWidth=1.2;
    end
    
    if EMGCnt == EMGThreshForFigures
        display(['mean LPO: ',num2str(mean(withinMiceNR_LPO))])
        display(['SEM  LPO: ',num2str(std(withinMiceNR_LPO)./sqrt(size(withinMiceNR_LPO,2)))])
        display(['mean  nLPO: ',num2str(mean(withinMiceNR_nonLPO))])
        display(['SEM nLPO: ',num2str(std(withinMiceNR_nonLPO)./sqrt(size(withinMiceNR_nonLPO,2)))])
    end
    
    
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('arousal delay [s]'); ax = gca; ax.XTickLabel ={'LPO','nonLPO'};
    set(gca,'TickDir','out'); box off
    ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415])
    pValsTTest(EMGCnt)=ranksum(withinMiceNR_nonLPO,withinMiceNR_LPO,'tail','both');
    linkaxes
    text(1,max(get(gca,'YLim'))-1,['p= ',num2str(pValsTTest(EMGCnt))])

end
saveas(f1,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\LPOvsnonLPO','tArousal_allDataPoints'))
saveas(f1,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\LPOvsnonLPO','tArousal_allDataPoints.png'))
saveas(f2,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\LPOvsnonLPO','tArousal_1DataPointPerMouse'))
saveas(f2,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\LPOvsnonLPO','tArousal_1DataPointPerMouse.png'))
%% compare arousal delays LPO vs nonLPO from REM 
f1=figure('name','mean arousal delays comparison');
f2=figure('name','difference between NREM REM comparuison');

tArousNR_nonLPO=opto_nonLPO.tArousal_NR; 
tArousNR_LPO=opto_LPO.tArousal_NR;
tArousR_nonLPO=opto_nonLPO.tArousal_R; 
tArousR_LPO=opto_LPO.tArousal_R; 
pValsTTest = nan(numel(EMGThresholds),1);

for EMGCnt = 1:numel(EMGThresholds) % cycle through various emg thresholds
    %extract the arousalT based on the current EMG threshold
    withinMiceNR_nonLPO = cellfun(@nanmean,tArousNR_nonLPO(EMGCnt,:));
    withinMiceNR_LPO = cellfun(@nanmean,tArousNR_LPO(EMGCnt,:));
    withinMiceR_nonLPO = cellfun(@nanmean,tArousR_nonLPO(EMGCnt,:));
    withinMiceR_LPO = cellfun(@nanmean,tArousR_LPO(EMGCnt,:));
    diffNRvREM_LPO = 100.*withinMiceR_LPO./withinMiceNR_LPO;
    diffNRvREM_nonLPO = 100.*withinMiceR_nonLPO./withinMiceNR_nonLPO;

    if EMGCnt == EMGThreshForFigures
        forSPSSLPO = ([withinMiceNR_LPO',withinMiceR_LPO']);
        forSPSSnonLPO =([withinMiceNR_nonLPO',withinMiceR_nonLPO']);
        display(['mean NREM/REM LPO: ',num2str(mean(forSPSSLPO,1))])
        display(['SEM NREM/REM LPO: ',num2str(std(forSPSSLPO,[],1)./sqrt(size(forSPSSLPO,1)))])
        display(['mean NREM/REM nLPO: ',num2str(mean(forSPSSnonLPO,1))])
        display(['SEM NREM/REM nLPO: ',num2str(std(forSPSSnonLPO,[],1)./sqrt(size(forSPSSLPO,1)))])
        signrankTestLPO = signrank(forSPSSLPO(:,1),forSPSSLPO(:,2));
        
        % instead of ANOVA in SPSS, do sepparate rank sum tests 
        
        toPlot = {100.*withinMiceR_LPO./withinMiceNR_LPO,100.*withinMiceR_nonLPO./withinMiceNR_nonLPO};
        titles = {'LPO','nonLPO'};
        figure;
        boxplot([toPlot{1},toPlot{2}],[ones(size(toPlot{1})),2.*ones(size(toPlot{2}))]);hold on;
        plot(1,toPlot{1},'-ok','MarkerFaceColor','k','MarkerSize',4); hold on;
        plot(2,toPlot{2},'-ok','MarkerFaceColor','k','MarkerSize',4);
        h = findobj(gca,'Tag','Box');
        h2= findobj(gca,'Tag','Outliers');
        h3= findobj(gca,'Tag','Median');
        pValsTTest(1)=ranksum(toPlot{1},toPlot{2},'tail','both');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',.4,'EdgeColor','none');
            h(j).Color = 'none'; % set box edge color
            h2(j).MarkerEdgeColor='none'; h2(j).MarkerFaceColor='none';
            h3(j).LineWidth=1.2;
        end
        box off; set(gca,'TickDir','out');
        ax = gca;  ax.XTickLabel ={'LPO','nonLPO'};      
        ylabel('relative REM arousal delay  [% NREM]')
        
        
        % plot difference nonREM/REM for LPO
        toPlot = {withinMiceNR_LPO,withinMiceR_LPO,withinMiceNR_nonLPO,withinMiceR_nonLPO};
        titles = {'LPO','LPO','nonLPO','nonLPO'};
        for plotCnt = [1,3]
            figure;
            boxplot([toPlot{plotCnt},toPlot{plotCnt+1}],[ones(size(toPlot{plotCnt})),2.*ones(size(toPlot{plotCnt+1}))]);hold on;
            plot(1:2,[toPlot{plotCnt};toPlot{plotCnt+1}],'-ok','MarkerFaceColor','k','MarkerSize',4); hold on;
            h = findobj(gca,'Tag','Box');
            h2= findobj(gca,'Tag','Outliers');
            h3= findobj(gca,'Tag','Median');
            pValsTTest(plotCnt+1)=ttest(toPlot{plotCnt},toPlot{plotCnt+1},'tail','both');
            for j=1:length(h)
                patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',.4,'EdgeColor','none');
                h(j).Color = 'none'; % set box edge color
                h2(j).MarkerEdgeColor='none'; h2(j).MarkerFaceColor='none';
                h3(j).LineWidth=1.2;
            end
            box off; set(gca,'TickDir','out');
            ax = gca;  ax.XTickLabel ={'NREM','REM'};
            ylabel('arousal delay  [s]'); title(titles{plotCnt})
        end

    end
    
    % plot difference LPO nonLPO for each mouse sepparately
    figure(f1);
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    plot(1:2,[withinMiceNR_LPO;withinMiceR_LPO],'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3); hold on;
%     plot(2,withinMiceR_LPO,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3);
    plot([3.5,4.5],[withinMiceNR_nonLPO;withinMiceR_nonLPO],'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3); hold on;
%     plot(4.5,withinMiceR_nonLPO,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3);
    
    
    
%     CiNRMedianLPO = 1.58*(iqr(withinMiceNR_LPO)./sqrt(numel(withinMiceNR_LPO)));
%     CiNRMedianNonLPO = 1.58*(iqr(withinMiceNR_nonLPO)./sqrt(numel(withinMiceNR_nonLPO)));
%     CiRMedianLPO = 1.58*(iqr(withinMiceR_LPO)./sqrt(numel(withinMiceNR_LPO)));
%     CiRMedianNonLPO = 1.58*(iqr(withinMiceR_nonLPO)./sqrt(numel(withinMiceNR_nonLPO)));

    CiNRMedianLPO = abs(median(withinMiceNR_LPO)-prctile(withinMiceNR_LPO,[25,75]));
    CiNRMedianNonLPO =  abs(median(withinMiceR_LPO)-prctile(withinMiceR_LPO,[25,75]));
    CiRMedianLPO =  abs(median(withinMiceNR_nonLPO)-prctile(withinMiceNR_nonLPO,[25,75]));
    CiRMedianNonLPO = abs(median(withinMiceR_nonLPO)- prctile(withinMiceR_nonLPO,[25,75]));
    errorbar(1,median(withinMiceNR_LPO),CiNRMedianLPO(1),CiNRMedianLPO(2),'ok','MarkerFaceColor','k')
    errorbar(2,median(withinMiceR_LPO),CiNRMedianNonLPO(1),CiNRMedianNonLPO(2),'ok','MarkerFaceColor','k')
    errorbar(3.5,median(withinMiceNR_nonLPO),CiRMedianLPO(1),CiRMedianLPO(2),'ok','MarkerFaceColor','k')
    errorbar(4.5,median(withinMiceR_nonLPO),CiRMedianNonLPO(1),CiRMedianNonLPO(2),'ok','MarkerFaceColor','k')
 
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('arousal delay [s]'); ax = gca; ax.XTickLabel ={'LPO-NR','LPO-R','nLPO-NR','nLPO-R'};ax.XTickLabelRotation =90;
    set(gca,'TickDir','out'); box off
    ax.XTick = [1,2,3.5,4.5]; ax.XLim = [0,6.5]; set(gcf,'Position',[205 484 1264 415])
%     pValsTTest(EMGCnt)=ranksum(withinMiceNR_nonLPO,withinMiceNR_LPO,'tail','both');
    linkaxes

    
    figure(f2)
    subplot(1,numel(EMGThresholds),EMGCnt); hold on;
    plot(1,diffNRvREM_LPO,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3); hold on;
    plot(2,diffNRvREM_nonLPO,'-o','Color',[.8,.8,.8],'MarkerFaceColor',[.8,.8,.8],'MarkerSize',3);
    CidiffMedianLPO = 1.58*(iqr(diffNRvREM_LPO)./sqrt(numel(diffNRvREM_LPO)));
    CidiffMedianNonLPO = 1.58*(iqr(diffNRvREM_nonLPO)./sqrt(numel(diffNRvREM_nonLPO)));
    errorbar(1,median(diffNRvREM_LPO),CidiffMedianLPO,'ok','MarkerFaceColor','k')
    errorbar(2,median(diffNRvREM_nonLPO),CidiffMedianNonLPO,'ok','MarkerFaceColor','k')
    

    
    
    title(['EMG threshold: ',num2str(EMGThresholds(EMGCnt))]);
    ylabel('REM -NREM arousal delay [s]'); ax = gca; ax.XTickLabel ={'LPO','nLPO'};ax.XTickLabelRotation =90;
    set(gca,'TickDir','out'); box off
    ax.XTick = [1,2]; ax.XLim = [0,3]; set(gcf,'Position',[205 484 1264 415])
    linkaxes
end
%% plot average opto-EMG variance 
figure;
iteration = 2; %which threshold to pklot
conditions = {opto_nonLPO,opto_LPO};
colors = {[0,0,0],npgCMap(1,:)};
for lowHigh = 1:2
    for mouseCnt = 1:nMice
        currCond = conditions{lowHigh};
        subplot(4,2,mouseCnt)
        timeInS = -nonLPO.settings.tPre:0.25:nonLPO.settings.tPost-0.25;
        meanVar= median(currCond.optoEvokedEMG_NR{iteration,mouseCnt},1);
        VarPtiles = prctile(currCond.optoEvokedEMG_NR{iteration,mouseCnt},[25,75],1);
        meanTArous=mean(currCond.tArousal_NR{iteration,mouseCnt}); %mean time to crossing "arousal" threshold
        currThresh = currCond.thresholdsUsed{iteration,mouseCnt};
        plot(timeInS,meanVar,'Color',colors{lowHigh}); hold on; % plot mean evoked variance
        patch([timeInS,fliplr(timeInS)],[VarPtiles(1,:),fliplr(VarPtiles(2,:))],'k','FaceAlpha',0.3,'EdgeColor','none','FaceColor',colors{lowHigh}); % add 25/75%percentiles as shading
        xlabel('time to stim [s]'); ylabel('log EMG var'); box off; set(gca,'TickDir','out','YScale','log');xlim([-10,20]); %cosmetics
        if mouseCnt ==1 && lowHigh==1
            legend('median nonLPO','25/75 %tiles','AutoUpdate','off'); legend boxoff          
        end
%         line([meanTArous,meanTArous],get(gca,'YLim'),'Color',colors{lowHigh})
%         line(get(gca,'XLim'),[currThresh,currThresh],'Color',colors{lowHigh})
    end
end
linkaxes
saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\LPOvsnonLPO','mean_evoked_EMG_Var'))
saveas(gcf,fullfile('F:\Tomoko_OptoStim\figures_EMGVar\LPOvsnonLPO','mean_evoked_EMG_Var.png'))