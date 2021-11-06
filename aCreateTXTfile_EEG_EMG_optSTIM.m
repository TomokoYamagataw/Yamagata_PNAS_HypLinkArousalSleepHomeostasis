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


% path='I:\OptoMod\'; 
path='E:\OptoInh\';  %mkdir(pathout);
pathsig=[path,'OutputSignals\'];  %mkdir(pathsig);
pathoutEDF=[path,'OutputEDFs\']; mkdir(pathoutEDF)

%enter the number of total channels the EDF will have (i.e. individual '.mat'-files) each mouse has

sname='est1'; %'est1'; %'est4'; %est1 for SD %est2 for 2h continuous stim  %est4 24hr 5min stim
numMatFiles = 4; %4 %5 %6 %7

recorddates=char('060421','070421','080421','110421','120421','150421','160421','170421');
% recorddates=char('050521','060521','100521','140521','150521','170521','180521');
rds=recorddates;
mousenames=char('Gf1','Ar1','Ar2','Ar3','Gf2','Ar4');
% mousenames=char('Ar8','Gf4','Gf5','Ar9'); 
% mousenames=char('Ar5','Gf3'); 
days=[8]; %days=[2]; 
mice=[2 3]; %mice=[2 4 6]; %


% with the below two settings exactly 24 hours of data will be used
maxep=10800; %21600 %10800
epochl=4;
fs=256;
% calculate number of samples within an exactly 24 h long recording (not sure why this is done so complicatedly)
fsh24=fs*epochl*maxep;


% this designs a chebyshev type II filter 
p1=0.5; p2=100; s1=0.1; s2=120;
Wp=[p1 p2]/(fs/2); Ws=[s1 s2]/(fs/2); Rp=3; Rs=30; [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bb1,aa1]=cheby2(n,Rs,Wn);

p1=5; p2=100; s1=4; s2=120;
Wp=[p1 p2]/(fs/2); Ws=[s1 s2]/(fs/2); Rp=3; Rs=30; [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bb2,aa2]=cheby2(n,Rs,Wn);

% p3=10; p4=100; s3=1; s4=120;
% WpM=[p3 p4]/(fs/2); WsM=[s3 s4]/(fs/2); RpM=3; RsM=30; [nM, WnM]=cheb2ord(WpM,WsM,RpM,RsM);
% [bbb1,aaa1]=cheby2(nM,RsM,WnM);


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
        
%         %%%% POAs
%         for chanEEG=4:6 %6
%             channame=channames(chanEEG,:)
%             fnin=[mousename,'-EEG',num2str(mouse),'-',recorddate,'-',num2str(channame)];
%             eval(['load ',pathsig,fnin,'.mat resampled_sig -mat']);
%             if length(resampled_sig)>fsh24
%                 signal=resampled_sig(1:fsh24);
%             else
%                 signal=zeros(1,fsh24);
%                 signal(1:length(resampled_sig))=resampled_sig;
%             end
%             signal=filtfilt(bb1,aa1,signal); 
%            
%             output(:,chanEEG)=signal';
%             clear resampled_sig signal;
%         end
% %         
        
         % 490 FP channel - processed version (the 'N' stands for new....)
         chanEEG=chanEEG+1;
         fnin=[mousename,'-',sname,'-',recorddate];%fnin=[mousename,'-dF_rawFitCorrected','-',recorddate,'-ch1'];
         eval(['load ',pathsig,fnin,'.mat resampled_sig -mat']);
         if length(resampled_sig)>fsh24
             signal=resampled_sig(1:fsh24);
         else
             signal=zeros(1,fsh24);
             signal(1:length(resampled_sig))=resampled_sig;
         end
        % signal=filtfilt(bb1,aa1,signal);
         signal(1,1:2)=100;
         output(:,chanEEG)=signal';
         clear resampled_sig signal;
         
         
         
        fnoutTXT=[mousename,'-EEG-EMG-',recorddate];%fnoutTXT=[mousename,'-EEG-EMG-dFR-dFS',recorddate];       
%         fnoutTXT=[mousename,'-EEG-EMG-POA-ChR2Optstim-',recorddate];%fnoutTXT=[mousename,'-EEG-EMG-dFR-dFS',recorddate];       
%         fnoutTXT=[mousename,'-EEG-EMG-Arch-Optstim-',recorddate];%fnoutTXT=[mousename,'-EEG-EMG-dFR-dFS',recorddate];
        eval(['save ',pathoutEDF,fnoutTXT,'.txt output -ascii']);
        
        
        
        %%%%% check program
%         eval(['save ',pathoutEDF,fnoutTXT,'2.mat output']);
%         clear all; path='I:\OptoMod\'; pathoutEDF=[path,'OutputEDFs\']; fnoutTXT=['Beck-EEG-EMG-POA-OPT4-step10000_251117']; load([pathoutEDF,fnoutTXT,'2.mat']);
        %%%%%
    end
end


% recorddates=char('080819'); %'110819'
% rds=recorddates; % rds=char('080819')
% mousenames=char('GDCh18','','GDCh20','');

% recorddates=char('110819'); 
% rds=recorddates;
% mousenames=char('','GDCh19','','GDCh21');
% 
% recorddates=char('160819','140819'); %'110819'
% rds=recorddates; % rds=char('080819')
% mousenames=char('GDCh18','GDCh19','GDCh20','GDCh21');
% days=[1]; %days=size(recorddates,1);
% mice=[1:4]; %numanim=size(mousenames,1);mice=1:numanim;

% recorddates=char('111018'); %'110819'
% rds=recorddates; % rds=char('080819')
% mousenames=char('GFP6','GFP7','GFP5','GFP5');
% days=[1]; %days=size(recorddates,1);
% mice=[4:4]; %numanim=size(mousenames,1);mice=1:numanim;

% sname='estm';
% sname='est2'; 