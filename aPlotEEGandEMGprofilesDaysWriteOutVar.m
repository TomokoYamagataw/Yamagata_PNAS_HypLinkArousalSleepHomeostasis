clear all;
close all;
outp=[];nums=[];

epochl=4;
path='E:\OptoInh\';
pathsig=[path,'OutputSignals\'];
pathV=[path,'OutputSIGvar\'];mkdir(pathV)
pathF=[path,'Figures_Profiles\'];mkdir(pathF)


%%%%%% 'Opto Inhibition'
year='21';
recorddates=strvcat('1704');
mousenames=strvcat('Gf1','Ar1','Ar2','ZZZ','Gf2','ZZZ'); %#ok<DSTRVCT>
mice=[2 3];%[1 3 5];[1:size(mousenames,1)];% numdays=[1:size(recorddates,1)];
numdays=[1:size(recorddates,1)];
artEEG=[1*10^5 1*10^5 1*10^5 1*10^5 1*10^5 1*10^5];
artEMG=[4*10^5 4*10^5 4*10^5 4*10^5 4*10^5 4*10^5];
yaxiss=[7*10^4 7*10^4 7*10^4 7*10^4 8*10^4 7*10^4];
emgamp=[8.0 8.0 8.0 8.0 10.0 8.0];


% %%%%%% 'Opto Inhibition'
% year='21';
% recorddates=strvcat('0704','0804','1204','1504');
% mousenames=strvcat('Gf1','ZZZ','Ar2','ZZZ','Gf2','ZZZ'); %#ok<DSTRVCT>
% mice=[5];%[1 3 5];[1:size(mousenames,1)];% numdays=[1:size(recorddates,1)];
% numdays=[1:size(recorddates,1)];
% artEEG=[1*10^5 1*10^5 1*10^5 1*10^5 1*10^5 1*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5 4*10^5 4*10^5];
% yaxiss=[7*10^4 7*10^4 7*10^4 7*10^4 8*10^4 7*10^4];
% emgamp=[8.0 8.0 8.0 8.0 10.0 8.0];

% %%%%%%%% 'Opto Inhibition'
% year='21';
% recorddates=strvcat('0704','0804','1204','1604');
% mousenames=strvcat('ZZZ','Ar1','ZZZ','Ar3','ZZZ','Ar4');
% mice=[6];%[2 4 6];[1:size(mousenames,1)];% numdays=[1:size(recorddates,1)];
% numdays=[1:size(recorddates,1)];
% artEEG=[1*10^5 1*10^5 1*10^5 1*10^5 1*10^5 1*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5 4*10^5 4*10^5];
% yaxiss=[7*10^4 7*10^4 7*10^4 7*10^4 7*10^4 1*10^5];
% emgamp=[8.0 8.0 8.0 8.0 8.0 16.0];
% 



% %%%%%%% 'Opto Inhibition'
% year='21';
% recorddates=strvcat('0505','0605','1005','1405','1505');
% mousenames=strvcat('Ar8','Gf4','Gf5','Ar9');
% mice=[1 2 3 4];%[1:size(mousenames,1)];% numdays=[1:size(recorddates,1)];
% numdays=[1:size(recorddates,1)];
% artEEG=[1*10^5 1*10^5 1*10^5 1*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5];
% yaxiss=[5.5*10^4 7*10^4 5.5*10^4 8*10^4];
% emgamp=[8.0 5.0 5.0 10.0];


% year='21';
% recorddates=strvcat('0505','0605','0905','1005','1405');
% mousenames=strvcat('Ar5','Gf3');
% mice=[1];%[1:size(mousenames,1)];% numdays=[1:size(recorddates,1)];
% numdays=[1:size(recorddates,1)];
% artEEG=[8*10^4 1*10^5];
% artEMG=[4*10^5 4*10^5];
% yaxiss=[7*10^4 7*10^4];
% emgamp=[5.0 5.0];




fs=256;
numm=60;
numh=12; %numh=24;
maxep=10800;  %21600
x1=1:maxep;x1=x1./900;
zermat=zeros(1,fs*numm*numm*numh);

zermat1=zeros(1,maxep);
zermat2=zeros(1,numh*numm);

% events=strvcat('EEG','estm');
channames=strvcat('fro','emg')
ch=1;

p1=0.5; p2=100; s1=0.1; s2=120;
Wp=[p1 p2]/(fs/2); Ws=[s1 s2]/(fs/2); Rp=3; Rs=20; [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bb1,aa1]=cheby2(n,Rs,Wn);

p1=5; p2=100; s1=4; s2=120;
Wp=[p1 p2]/(fs/2); Ws=[s1 s2]/(fs/2); Rp=3; Rs=20; [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bb2,aa2]=cheby2(n,Rs,Wn);

for mouse=mice
    
    mousename=mousenames(mouse,:);
    mousename(isspace(mousename))=[];
    
    figure
    
    for dd=numdays
        recorddate=[recorddates(dd,:),year,'D4']
        channame=channames(1,:)
        %EEG
        fnout=[mousename,'-EEG',num2str(mouse),'-',recorddate,'-',num2str(channame)];
        eval(['load ',pathsig,fnout,'.mat resampled_sig -mat']);
        if length(resampled_sig)>length(zermat); resampled_sig=resampled_sig(1:length(zermat)); EEG=resampled_sig; 
        else;  EEG=zermat; EEG(1:length(resampled_sig))=resampled_sig; 
        end
        clear sig; EEG=filtfilt(bb1,aa1,EEG); EEG=filtfilt(bb1,aa1,EEG);
        %pause
        EEG(find(abs(EEG)>artEEG(mouse)))=0; EEG(EEG==0)=NaN;
        EEGv=nanvar(reshape(EEG,fs*epochl,maxep)); clear EEG;
        
        %EMG
        channame=channames(2,:)
        fnout=[mousename,'-EEG',num2str(mouse),'-',recorddate,'-',num2str(channame)];
        eval(['load ',pathsig,fnout,'.mat resampled_sig -mat']);
        if length(resampled_sig)>length(zermat); resampled_sig=resampled_sig(1:length(zermat)); EMG=resampled_sig; else;  EMG=zermat; EMG(1:length(resampled_sig))=resampled_sig; end
        clear sig; EMG=filtfilt(bb2,aa2,EMG); EMG=filtfilt(bb2,aa2,EMG);
        
        EMG(find(abs(EMG)>artEMG(mouse)))=0; EMG(EMG==0)=NaN;
        EMGv=nanvar(reshape(EMG,fs*epochl,maxep)); clear EMG;
        
%         fnout1=[mousename,'-EEGfrontal-EMGv-',recorddate];
%         eval(['save ',pathV,fnout1,'.mat mousename EEGv EMGv -mat']);
        
%         EEGv(EEGv>yaEEG(mouse))=NaN;
%         EMGv(EMGv>yaEMG(mouse))=NaN;
        
        subplot ('position',[0.1 0.95-0.16*dd 0.8 0.14])
        plot(x1,EEGv,'LineWidth',1,'Color',[0 0.4 1])
        hold on
        bar(x1,EMGv*emgamp(mouse),'r')
        axis([0 12 0 yaxiss(mouse)])%yaEEG(mouse)%1*10^5
        set(gca,'XTick',[0:2:12]) 
        grid on
%         if dd>1
%             plot([4 4],[0 yaEEG(mouse)],'-k','LineWidth',2);
            %plot([8 8],[0 yaEEG(mouse)],'-k','LineWidth',2);
%         end
%         if dd==1; title([mousename,' OptInhibition']); end
%         if dd==1; title([mousename,' baseline, 10Hz']); end
        text(0.2,1*10^4,recorddate,'FontWeight','bold')
        if dd==max(numdays); xlabel('Hours'); end
        % axis off
        
    end;
%     pause
    orient landscape %tall
    
    % GDCh 18-21
    % figname=[mousename,'-080819-250819-WakeEnhancement']
    % figname=[mousename,'-160819-210819-Baseline1Hz2Hz5Hz10Hz']
%     figname=[mousename,'-130518-160518-EEGfProfile']
    % figname=[mousename,'-100819_220819-Dexmedetomidine'] 
%     figname=[mousename,'-080819-280819-WakeEnhance_Caffeine']    
%     figname=[mousename,'-040919-060919-HighSleepPressure']
%     figname=[mousename,'-300819-080919-Continuous10HzStim']   
%     figname=[mousename,'-150819-170819-Plasticity']   
%     figname=[mousename,'-020919-140919-20Hz']   
%     figname=[mousename,'-120518-160518-10Hz']   
%     figname=[mousename,'-070818-040918-new10Hz']       
%     figname=[mousename,'-070421-160421-optinhSD-24hrstim']     
%     figname=[mousename,'-050521-150521-optinhSD-24hrstim']    
    figname=[mousename,'-170421D4-opt2hrstim']     
   


    saveas(gcf,[pathF,figname],'tiff')
%     close all
    
end

% 
% recorddates=strvcat('0708','0808','0908','1008','1108'); % Baseline
% recorddates=strvcat('0708','1008','1408','1508','1608'); % 10Hz 8sec 2min
% recorddates=strvcat('1308','1408','1508','1608','1708'); % 10Hz 8sec 2min
% recorddates=strvcat('1808','1908','2008','2108','2208'); % Baseline 1Hz 2Hz 5Hz
% recorddates=strvcat('0808','1108','2508','2808'); % WakeEnhancement_caffeine
% recorddates=strvcat('1008','2208'); % Dexmedetomidine
% recorddates=strvcat('0409','0609'); % HighSleepPressure
% recorddates=strvcat('3008','3108','0109','0809'); % Continuous10HzStim
% recorddates=strvcat('1008','1508','1608','1708'); % Plasticity, GDCh18,19,21
% numdays=[1:size(recorddates,1)];
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
% mice=[1:size(mousenames,1)];

% %%%%%%%%%%% GDCh18','GDCh19','GDCh20','GDCh21
% year='19';
% recorddates=strvcat('3008','3108','0109','0809'); % Continuous10HzStim
% numdays=[1:size(recorddates,1)];
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
% mice=[1:size(mousenames,1)];
% 
% artEEG=[4*10^5 4*10^5 4*10^5 4*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5];
% % yaEEG=[4*10^5 4*10^5 4*10^5 4*10^5];
% % yaEMG=[1*10^5 1*10^5 1*10^5 1*10^5];
% yaxiss=[1*10^5 1*10^5 5*10^4 8*10^4];
% emgamp=[10.0 8.0 2.0 4.0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plasticity, GDCh18,19,21
%%%%%%%%%%% GDCh18','GDCh19','GDCh20','GDCh21
% year='19';
% recorddates=strvcat('1008','1508','1608','1708'); % Plasticity, GDCh18,19,21
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
% mice=[1 2 4];%1 2 4 [1:size(mousenames,1)];
% numdays=[1:size(recorddates,1)];
% artEEG=[4*10^5 4*10^5 4*10^5 4*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5];
% yaxiss=[1*10^5 1*10^5 5*10^4 8*10^4];
% emgamp=[10.0 8.0 2.0 4.0];


% %%%%%%%%% 'GDCh20'
% year='19';
% recorddates=strvcat('0708','1008','2708','0209','0309','1409'); % Plasticity, GDCh20
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
% mice=[3];%[1:size(mousenames,1)];% numdays=[1:size(recorddates,1)];
% numdays=[1:size(recorddates,1)];
% artEEG=[4*10^5 4*10^5 4*10^5 4*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5];
% yaxiss=[1*10^5 1*10^5 4*10^4 8*10^4];
% emgamp=[10.0 8.0 0.5 4.0];

% %%%%%%%%%%% 'GDCh14,15,16,17'
% year='18';
% recorddates=strvcat('1812','1912','2012'); % Plasticity, GDCh18,19,21
% mousenames=strvcat('GDCh14','GDCh15','GDCh16','GDCh17');
% mice=[2:4];%[1:size(mousenames,1)];
% numdays=[1:size(recorddates,1)];
% artEEG=[2*10^5 2*10^5 2*10^5 2*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5];
% yaxiss=[1*10^5 1.2*10^5 1.2*10^5 1.2*10^5];
% emgamp=[1.0 20.0 10.0 0.01];


% %%%%%%%%%%% 'GDCh5'
% year='18';
% recorddates=strvcat('0204','0304','0404'); % Plasticity, GDCh18,19,21
% mousenames=strvcat('GDCh','GDCh5');
% mice=[2];%[1:size(mousenames,1)];
% numdays=[1:size(recorddates,1)];
% artEEG=[2*10^5 2*10^5 2*10^5 2*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5];
% yaxiss=[1*10^5 1*10^5 1*10^5 1*10^5];
% emgamp=[1.0 4.0 1.0 1.0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%% 'GDCh7','GDCh8'
% year='18';
% recorddates=strvcat('1205','1305','1405','1505','1605'); % Continuous10HzStim
% numdays=[1:size(recorddates,1)];
% mousenames=strvcat('GDCh7');
% mice=[1:2];
% artEEG=[2*10^5 2*10^5];
% artEMG=[2*10^5 2*10^5];
% yaxiss=[1*10^5 8*10^4];
% emgamp=[10 0.02];

%%%%%%%%%%% GDCh16
% year='19';
% recorddates=strvcat('0801','1501'); % WakeEnhancement
% numdays=[1:size(recorddates,1)];
% mousenames=strvcat('na','na','GDCh16','na');
% numdays=[1:size(recorddates,1)];

% mice=[3];

% artEEG=[4*10^5 4*10^5 4*10^5 4*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5];
% % yaEEG=[4*10^5 4*10^5 4*10^5 4*10^5];
% % yaEMG=[1*10^5 1*10^5 1*10^5 1*10^5];
% yaxiss=[1*10^5 1*10^5 8*10^4 1*10^5];
% emgamp=[10.0 5.0 5.0 1.0];



%%%%%%%%%%% GFP6 & GFP7
% year='18';
% recorddates=strvcat('0810','2210'); % WakeEnhancement
% numdays=[1:size(recorddates,1)];
% mousenames=strvcat('GFP6','GFP7');
% mice=[1:2];


% %%%%%%%%%% GFP5
% year='18';
% recorddates=strvcat('0110','0310','0410'); % WakeEnhancement
% numdays=[1:size(recorddates,1)];
% mousenames=strvcat('GFP','GFP','GFP5');
% mice=[3];
% artEEG=[2*10^5 2*10^5 2*10^5 2*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5];
% yaxiss=[1*10^5 1*10^5 5*10^4 1*10^5];
% emgamp=[1.0 1.0 1.0 1.0];

% %%%%%%%%%% GFP4
% year='18';
% recorddates=strvcat('3006','0107','0207','0307','0407'); % WakeEnhancement
% numdays=[1:size(recorddates,1)];
% mousenames=strvcat('GFP','GFP4');
% mice=[2];
% artEEG=[2*10^5 2*10^5 2*10^5 2*10^5];
% artEMG=[4*10^5 4*10^5 4*10^5 4*10^5];
% yaxiss=[1*10^5 1*10^5 1*10^5 1*10^5];
% emgamp=[1.0 8.0 1.0 1.0];
