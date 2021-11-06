clear all;
close all;

path='I:\optogenetics\';
% 
% mousenames=[5 6 8 15 16 17]; %GDCh1 - event1 %GDCh9 - event4
% event=[2 1 3 2 3 4]; 
% days=['030418';'120418';'150618';'191218';'191218';'191218'];

mousenames=[5 15 16 17 18 19 21]; %GDCh1 - event1 %GDCh9 - event4
event=[2 2 3 4 1 2 4]; 
days=['030418';'191218';'191218';'191218';'160819';'160819';'160819'];

pathsig=['I:\OptoMod\OutputSignals\'];
% pathoutEMG=[path,'OutputEMG\']; mkdir(pathoutEMG)
pathEMGstim=[path,'EMGstim2\']; mkdir(pathEMGstim)
pathstim=[path,'STIMs\']; % pathstim=[path,'outputSTIM24\']; 

maxep=21600; epochl=4;
fs=256;
fs1=100;
fsh24=fs*60*60*24;

before=120;
after=1200;
di=before:after;
period=(before+after);

sr=0.25; %64/256Hz = 0.25

p1=5; p2=100; s1=4; s2=120;
Wp=[p1 p2]/(fs/2); Ws=[s1 s2]/(fs/2); Rp=3; Rs=30; [n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bb1,aa1]=cheby2(n,Rs,Wn);

numanim=length(mousenames);

for n=5:numanim
    mousename=['GDCh',num2str(mousenames(n))];
    day=days(n,:); day(isspace(day))=[];
    count=1;

%     if n==3
%           fnin=[mousename,'-vis',num2str(event(n)),'-',day,'-emg'];
%     else
          fnin=[mousename,'-EEG',num2str(event(n)),'-',day,'-emg'];
%     end
    
    eval(['load ',pathsig,fnin,'.mat resampled_sig -mat']);

    if length(resampled_sig)>fsh24
        signal=resampled_sig(1:fsh24); clear resampled_sig;
    else
        signal=zeros(1,fsh24);
        signal(1:length(resampled_sig))=resampled_sig; clear resampled_sig;
    end
    signal=filtfilt(bb1,aa1,signal);

    fn2=[mousename,'_',day,'_stim'];
    eval(['load ',pathstim,fn2,'.mat startend -mat']);

    nums=size(startend,1);

    EMG250ms=zeros(nums,period*4);
    EMG1s=zeros(nums,period);
    for s=1:nums;
        st=startend(s,1)-before*fs;
        en=startend(s,1)+after*fs;
        if st<1 continue; end
        if en>21600*256*4 continue; end
        emg1=var(reshape(signal(st+1:en),64,[]));%600
        EMG250ms(s,:)=emg1;
        emg2=var(reshape(signal(st+1:en),fs,[]));%150
        EMG1s(s,:)=emg2;
    end

    fnoutEMG=[mousename,'-EMGstim-',day];

    eval(['save ',pathEMGstim,fnoutEMG,'.mat EMG1s EMG250ms startend -mat']);
    clear EMG1s EMG250ms signal;

    %     orient tall
    %     fign=[pathfig,[mousename,'-EMGvar']]
    %     saveas(gcf,fign,'tiff')
    %     close all

end
