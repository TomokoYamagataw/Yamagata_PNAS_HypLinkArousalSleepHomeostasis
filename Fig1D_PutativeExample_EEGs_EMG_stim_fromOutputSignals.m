% path='I:\';      
% pathEDF=[path,'OptoMod\OutputEDFs\']; 
% pathstim=[path,'optogenetics\STIMs\'];%mkdir(pathstim)
% mousename=char('GDCh21');
% recorddate=char('160819');
% fnoutTXT=[mousename,'-EEG-EMG-',recorddate];%fnoutTXT=[mousename,'-EEG-EMG-dFR-dFS',recorddate];       
% eval(['load ',pathEDF,fnoutTXT,'.txt output -ascii']);
% mt=int2str('mousename','_EEG_EMG_','recorddate');


clear all;
close all;


%%%%% COMMENT IN/OUT %%%poa line 104 if the animal has poa LFP

% user input
channames=char('fro','occ','emg')
% channames=strvcat('fro','occ','emg','po1');% channames=strvcat('fro','occ','emg1','405','490'); % channels to extract
% channames=strvcat('fro','occ','emg1','405','490'); % channels to extract
% channames=strvcat('fro','occ','emg','poa','lfr');% channames=strvcat('fro','occ','emg1','405','490'); % channels to extract
% channames=char('fro','occ','emg','po1','pod');
% channames=char('fro','occ','emg','po1','po2');
% channames=strvcat('fro','occ','emg','po1','po2','pod');
% channames=strvcat('fro','occ','emg','po2');
% channames=strvcat('fro','occ','emg','po1','po2','po3');
% channames=char('fro','occ','emg','poa')


path='I:\'; % specify the location of the data tank (named 'tank below') - NOTE: end with \
pathsig=[path,'OptoMod\OutputSignals\']; 
pathstim=[path,'optogenetics\STIMs\'];
pathfig=[path,'optogenetics\Figures_EMG_EEG_STIM\']; mkdir(pathfig);

%enter the number of total channels the EDF will have (i.e. individual '.mat'-files) each mouse has

sname='estm';% sname='est2'; 
numMatFiles = 4; %4 %5 %6 %7

recorddates=char('220819','060919'); %Dex:220819, HSP:110919, 2Hz:200819, 5Hz:210819
% recorddates=char('160819');  %baseline:160819
rds=recorddates;
mousenames=char('GDCh18','GDCh19','GDCh20','GDCh21');
days=[1:1]; %days=size(recorddates,1);
mice=[1 4]; %numanim=size(mousenames,1);mice=1:numanim;
% sname=['est',int2str(mice)];

% with the below two settings exactly 24 hours of data will be used
maxep=21600; %21600 %10800
epochl=4;
fs=256;
% calculate number of samples within an exactly 24 h long recording (not sure why this is done so complicatedly)
fsh24=fs*epochl*maxep;


% this designs a chebyshev type II filter 
p1=0.5; p2=30; s1=0.4; s2=40;
Wp=[p1 p2]/(fs/2); Ws=[s1 s2]/(fs/2); Rp=3; Rs=30; [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bb1,aa1]=cheby2(n,Rs,Wn);

p1=10; p2=100; s1=4; s2=120;
Wp=[p1 p2]/(fs/2); Ws=[s1 s2]/(fs/2); Rp=3; Rs=30; [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bb2,aa2]=cheby2(n,Rs,Wn);


for day=days
    recorddate=[recorddates(day,:)];
    rd=rds(day,:);
    
    for mouse=mice
        
        mousename=mousenames(mouse,:);   mousename(isspace(mousename))=[];
        output=zeros(fsh24,numMatFiles);
        
        % EEG original
        for chanEEG=1:2
            channame=channames(chanEEG,:)
            fnin=[mousename,'-EEG',num2str(mouse),'-',recorddate,'-',num2str(channame)];
            eval(['load ',pathsig,fnin,'.mat resampled_sig -mat']);
            if length(resampled_sig)>fsh24
                signal=resampled_sig(1:fsh24);
            else
                signal=zeros(1,fsh24);
                signal(1:length(resampled_sig))=resampled_sig;
            end

            signal=filtfilt(bb1,aa1,signal); 
            %   signal(abs(signal)>2000)=0;
           
            output(:,chanEEG)=signal';
            clear resampled_sig signal;
        end
    
        
        % EMG
        chanEEG=chanEEG+1;
        channame=channames(chanEEG,:)
        fnin=[mousename,'-EEG',num2str(mouse),'-',recorddate,'-',num2str(channame)];
        
        eval(['load ',pathsig,fnin,'.mat resampled_sig -mat']);
        
        if length(resampled_sig)>fsh24
            signal=resampled_sig(1:fsh24);
        else
            signal=zeros(1,fsh24);
            signal(1:length(resampled_sig))=resampled_sig;
        end
        
        %  what is this? - shouldnt rteally matter though as it just affects the very first datapoint 
        if max(abs(signal))==0
            signal(1)=1;
        end

        signal=filtfilt(bb2,aa2,signal); %  signal=filtfilt(bbb1,aaa1,signal); 
%         hold on
%         plot(signal,'-m')

        output(:,chanEEG)=signal';
        clear resampled_sig signal;
        %         
        
         % Optical stimulation
         chanEEG=chanEEG+1;
         fnin=[mousename,'-',sname,'-',recorddate];%fnin=[mousename,'-dF_rawFitCorrected','-',recorddate,'-ch1'];
         eval(['load ',pathsig,fnin,'.mat resampled_sig -mat']);
         if length(resampled_sig)>fsh24
             signal=resampled_sig(1:fsh24);
         else
             signal=zeros(1,fsh24);
             signal(1:length(resampled_sig))=resampled_sig;
         end
         output(:,chanEEG)=signal';
         clear resampled_sig signal;
         

        % cols='kkkkkkkk'; % colors
        % figure(1)
        % plot(output(1:1024,1:4),'Color',cols(1))      


        fn1=[mousename,'_',recorddate,'_stim'];        
        eval(['load ',pathstim,fn1,'.mat sigstim startend -mat']);

        tfstim=1;
        tnstim=72; 

        before=4;  % [sec] before stimulus onset 
        after=45;   % [sec] after stimulus onset %Baseline 8, Dex 1st stim: 20, 2nd-3rd stim: 45-60

        stimep1=round(startend(:,1)); % stim start
        stimep2=round(startend(:,2)); % stim end
        % stimep1=round(startend(:,1)./(fs*epochl)); % stim start (stim episode 1)
        % stimep2=round(startend(:,2)./(fs*epochl)); % stim end (episode 2) 

        stimep=stimep1; %start of stimulation episodes
        numstim=size(stimep,1)

        % eval(['load ',pathvs,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);%
        % 
        % VS=zeros(7,maxep); 
        % VS(1,w)=1; VS(2,w1)=1; VS(3,nr)=1; VS(4,nr2)=1; VS(5,r)=1; VS(6,r3)=1; VS(7,mt)=1;
        % clear w nr r w1 nr2 r3 mt ma bastend ; %clears the variables for the next loop 
        % 
        % wake=sum(VS(1:2,:));     
        % nrem=sum(VS(3:4,:));     
        % rem=sum(VS(5:6,:)); 
        % move=VS(7,:);

        mats=[];
        xt=0:1:(before+after)*fs-1; xt=xt/256;%xt=1:1:nums*fs1; xt=xt/fs1; % x time axes with steps of 1s/fs1=1/256=0.0039, matrix

%         if day<2; ns=18; else; ns=8; end
        for s=[2]%[2 7 10]%[1 8]%[7 10] %numstim %iterating over all the stimulations
            
            %GDCh18: [7 10],
            
            step=stimep(s); 
            eps=step-before*fs:step+after*fs-1; %returns the desired period

            mat=output(eps,:);
            %%%%%%%%%%%%% plot Frontal EEG %%%%%%%%%%%%%%
            figure(1000*day+100*mouse+s)

            subplot('position',[0.2,0.8,0.7,0.15])
            plot(xt,mat(:,1),'-k','LineWidth',.5);%plot(xt,EEG(4,:),'-k','LineWidth',0.5);
%             text(-4,0,'Frontal')
            hold on;
%             thr=max(abs(mat(:,1)))*1;
            axis([-1 max(xt) -900 700]) %baseline 400, Dex -900, 700
            axis off                   

            subplot('position',[0.2,0.64,0.7,0.15])
            plot(xt,mat(:,2),'-k','LineWidth',.5);%plot(xt,EEG(4,:),'-k','LineWidth',0.5);
%             text(-4,0,'Occipital')
            hold on;
            plot([-0.2 -0.2],[-800 0],'-k','LineWidth',2)
%             text(-nums/10,-400,'200uV')
%             thr=max(abs(mat(:,2)))*1; %1.1
            axis([-1 max(xt) -900 700])%1800
            axis off                   

            %%%%%%%%%%%%% plot EMG %%%%%%%%%%%%%%
            subplot('position',[0.2,0.52,0.7,0.1])
            plot(xt,mat(:,3),'-k','LineWidth',.5);%plot(xt,EEG(4,:),'-k','LineWidth',0.5);
            hold on
%             text(-4,50,'EMG');
            plot([-0.2 -0.2],[-400 0],'-k','LineWidth',2)
            %text(-nums/10,-400,'200uV')
%             thr=max(abs(mat(:,3)))*1;
            axis([-1 max(xt) -400 400])%1800
            axis off  

            %%%%%%%%%%%%% plot Optical stimulation %%%%%%%%%%%%%%
            subplot('position',[0.2,0.40,0.7,0.08])
            plot(xt,mat(:,4),'-k','LineWidth',0.5);%plot(xt,EEG(4,:),'-k','LineWidth',0.5);
%             text(-4,70,'Stimulation');
            %plot([-0.2 -0.20],[-500 -200],'-k','LineWidth',2)
            %text(-nums/10,-400,'200uV')
            thr=max(mat(:,4))*2;
            axis([-1 max(xt) -thr/10 thr])%1800
%             axis off  

%             figname=[mousename,'-',recorddate,'_',num2str(s),'th_stim_Dex_',num2str(floor(step/fs/4))];
            figname=[mousename,'-',recorddate,'_',num2str(s),'th_stim_10Hz_scale_',num2str(floor(step/fs/4))];
%             saveas(gcf,[pathfig,figname],'tif'); %tiffn
            saveas(gcf,[pathfig,figname],'svg'); %tiffn
%             saveas(gcf,[pathfig,figname],'fig'); %tiffn
            
            close all
    %         mats=[mats mat1];
        end
    end

end