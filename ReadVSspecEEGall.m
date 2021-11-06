% 
clear all
close all
% 
exte='.txt'
path='I:\OptoMod\';
discrip=strvcat('all_30Hz'); %% der=strvcat('fro_30Hz'); %%%%%all_30Hz %all_30Hz_bug
ders=strvcat('fro','occ') %% ders=strvcat('fro') %% ders=strvcat('fro','occ','stm')
numder=size(ders,1);


%%%%%%%%%%%%%%%% GDCh8
%%%%%%%%% Baseline/24h 10Hz 2min stim
% pathin=[path,'fft30Hz\Baseline-24h-fft\']
% pathout=[path,'outputVS1\Baseline_10Hz\']; %mkdir(pathout)
% pathfig=[path,'Figure_VSspec\Baseline_10Hz\']; %mkdir(pathfig)
% maxep=21600; %21600 %10800 

% mousenames=strvcat('GDCh8','GDCh8');
% days=strvcat('130518','140518');



%%%%%%%%%%%%%%%% GDCh18-21
%%%%%%%%% Baseline/24h 10Hz 2min stim
% pathin=[path,'fft30Hz\Baseline-24h-fft\'];
% pathout=[path,'outputVS1\Baseline_10Hz\']; %mkdir(pathout)
% pathfig=[path,'Figure_VSspec\Baseline_10Hz\']; %mkdir(pathfig)
% maxep=21600; %21600 %10800 
% % mousenames=strvcat('GDCh5','GDCh18','GDCh18','GDCh19','GDCh19','GDCh20','GDCh20','GDCh21','GDCh21');
% % days=strvcat('030418','160819','100819','160819','070819','160819','100819','160819');
% 
% days=char('','','180919');
% mousenames=char('','','GDCh20');

%%%%%%%%%% Baseline-12h
% pathin=[path,'fft30Hz\Baseline-24h-fft\']
% pathout=[path,'outputVS1\Baseline_10Hz\']; %mkdir(pathout)
% pathfig=[path,'Figure_VSspec\Baseline_10Hz\']; %mkdir(pathfig)
% maxep=10800; 
% mousenames=strvcat('GDCh18','GDCh19');
% days=strvcat('070819','070819');
% discrip=strvcat('12h_30Hz');


% % %%%%%%%% 1Hz2Hz5Hz
% pathin=[path,'fft30Hz\1Hz2Hz5Hz\'];
% pathout=[path,'outputVS1\1Hz2Hz5Hz\']; mkdir(pathout)
% pathfig=[path,'Figure_VSspec\1Hz2Hz5Hz\']; mkdir(pathfig)
% maxep=21600;
% mousenames=strvcat('GDCh18','GDCh18','GDCh18','GDCh19','GDCh19','GDCh19','GDCh20','GDCh20','GDCh20','GDCh21','GDCh21','GDCh21');
% days=strvcat('190819','200819','210819','190819','200819','210819','190819','200819','210819','190819','200819','210819');


%%%%%%% 2h1Hz
% pathin=[path,'fft30Hz\2h1Hz\']
% pathout=[path,'outputVS1\2h1Hz\']; mkdir(pathout)
% pathfig=[path,'Figure_VSspec\2h1Hz\']; mkdir(pathfig)
% maxep=21600; 
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
% days=strvcat('110919','110919','110919','110919');


%%%%%%%% 8sec
% pathin=[path,'fft30Hz\8sec\']
% pathout=[path,'outputVS1\8sec\']; mkdir(pathout)
% pathfig=[path,'Figure_VSspec\8sec\']; mkdir(pathfig)
% maxep=21600; 
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
% days=strvcat('140819','140819','140819','140819');


%%%%%%% Continuous ZT0
% pathin=[path,'fft30Hz\Continuous\']
% pathout=[path,'outputVS1\Continuous\']; mkdir(pathout)
% pathfig=[path,'Figure_VSspec\Continuous\']; mkdir(pathfig)
% maxep=10800; 
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21','GDCh18','GDCh19','GDCh20','GDCh21');
% days=strvcat('010919','010919','010919','010919','080919','080919','080919','080919');


%%%%%%%% Continuous ZT12
% pathin=[path,'fft30Hz\Continuous\']
% pathout=[path,'outputVS1\Continuous\']; mkdir(pathout)
% pathfig=[path,'Figure_VSspec\Continuous\']; mkdir(pathfig)
% maxep=21600; 
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
% days=strvcat('160919','160919','160919','160919');



% %%%%%%%% Caffeine
% pathin=[path,'fft30Hz\Caffeine\'];
% pathout=[path,'outputVS1\Caffeine\']; mkdir(pathout)
% pathfig=[path,'Figure_VSspec\Caffeine\']; mkdir(pathfig)
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21','GDCh18','GDCh19','GDCh20','GDCh21');
% days=strvcat('250819','250819','250819','250819','280819','280819','280819','280819');
% maxep=10800;


% %%%%%%%% SD
% pathin=[path,'fft30Hz\SD\'];
% pathout=[path,'outputVS1\SD\']; mkdir(pathout)
% pathfig=[path,'Figure_VSspec\SD\']; mkdir(pathfig)
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21','GDCh18','GDCh19','GDCh20','GDCh21','GDCh18','GDCh19','GDCh20','GDCh21');
% days=strvcat('040919','040919','040919','040919','060919','060919','060919','060919','100919','100919','100919','100919');
% maxep=10800;


% %%%%%%%% WakeEnhancement
% pathin=[path,'fft30Hz\WakeEnhance\'];
% pathout=[path,'outputVS1\WakeEnhance\']; %mkdir(pathout)
% pathfig=[path,'Figure_VSspec\WakeEnhance\']; %mkdir(pathfig)
% mousenames=strvcat('GDCh18','GDCh18','GDCh19','GDCh19','GDCh20','GDCh20','GDCh21','GDCh21');
% days=strvcat('080819','110819','080819','110819','080819','110819','080819','110819');
% maxep=10800;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% %%%%%%% 20Hz
% pathin=[path,'fft30Hz\20Hz\']
% pathout=[path,'outputVS1\20Hz\']; mkdir(pathout)
% pathfig=[path,'Figure_VSspec\20Hz\']; mkdir(pathfig)
% maxep=21600; 
% mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
% days=strvcat('140919','140919','140919','140919');
% 


%%%%%%%%%% Baseline
% pathin=[path,'fft30Hz\Baseline-24h-fft\'];
% pathout=[path,'outputVS1\Baseline_10Hz\']; %mkdir(pathout)
% pathfig=[path,'Figure_VSspec\Baseline_10Hz\']; %mkdir(pathfig)
% maxep=21600;
% mousenames=strvcat('GFP3');
% days=strvcat('260618');
% mousenames=strvcat('GFP5');
% days=strvcat('051018');
% mousenames=strvcat('GFP1','GFP1','GFP5','GFP5');
% days=strvcat('140518','160518','011018','041018');


%%%%%%% additional Baseline
pathin=[path,'fft30Hz\Baseline-24h-fft\'];
pathout=[path,'outputVS1\Baseline_10Hz\']; %mkdir(pathout)
pathfig=[path,'Figure_VSspec\Baseline_10Hz\']; %mkdir(pathfig)
maxep=21600; 
mousenames=strvcat('GDCh18','GDCh19','GDCh20','GDCh21');
days=strvcat('030919','030919','030919','030919');




numanim=size(mousenames,1);

f=0:0.25:30;
zermat=zeros(1,maxep);
% pathout=[path,'outputVS1\']; mkdir(pathout)


OOO=[];
for n=[3]%1:numanim%1:numanim %1:numanim %11 %[3]
    mouse=mousenames(n,:); mouse(isspace(mouse))=[];
    day=days(n,:);

    figure
    
    fname=[mouse,'-EEG-EMG-',day,'_',discrip]    
%     fname=[mouse,'-EEG-EMG-ChR2Optstim-',day,'_',der]
%   fname=[mouse,'-EEG-EMG-POA-ChR2Optstim-',day,'_',der]

    fnameFFT=[pathin,fname,exte]

    numline=1;
    fidfft=fopen(fnameFFT,'r')

    if fidfft<1 OOO=[OOO; fname]; continue;  end;

    fl=textread(fnameFFT,'%s%*[^\n]');
    fl1=char(fl);
    ep=find(fl1=='E');
    if size(ep,1)>1; numrow=ep(2)-ep(1)-2; else; numrow=maxep; end

    for d=1:numder %1:3

        numskip=ep(d)+d;

        [epoch,state,dateexp,time,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,... % NO DELTA ALPHA THETA IN LINW 38 AND 42
            a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,...
            a42,a43,a44,a45,a46,a47,a48,a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,...
            a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80,a81,a82,...
            a83,a84,a85,a86,a87,a88,a89,a90,a91,a92,a93,a94,a95,a96,a97,a98,a99,a100,a101,...
            a102,a103,a104,a105,a106,a107,a108,a109,a110,a111,a112,a113,a114,a115,a116,a117,a118,a119,a120,a121]...
            = textread(fnameFFT,'%d%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',numrow,'headerlines',numskip);

        spectr=[a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,...
            a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,...
            a42,a43,a44,a45,a46,a47,a48,a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,...
            a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80,a81,...
            a82,a83,a84,a85,a86,a87,a88,a89,a90,a91,a92,a93,a94,a95,a96,a97,a98,a99,a100,a101,...
            a102,a103,a104,a105,a106,a107,a108,a109,a110,a111,a112,a113,a114,a115,a116,a117,a118,a119,a120,a121];

        %spectr=spectr(:,1:81);

        clear a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20;
        clear a21 a22 a23 a24 a25 a26 a27 a28 a29 a30 a31 a32 a33 a34 a35 a36 a37 a38 a39 a40 a41;
        clear a41 a42 a43 a44 a45 a46 a47 a48 a49 a50 a51 a52 a53 a54 a55 a56 a57 a58 a59 a60;
        clear a61 a62 a63 a64 a65 a66 a67 a68 a69 a70 a71 a72 a73 a74 a75 a76 a77 a78 a79 a80 fft;
        clear a81 a82 a83 a84 a85 a86 a87 a88 a89 a90 a91 a92 a93 a94 a95 a96 a97 a98 a99 a100;
        clear a101 a102 a103 a104 a105 a106 a107 a108 a109 a110 a111 a112 a113 a114 a115 a116 a117 a118 a119 a120 a121;

        statestr=char(state);

        nr=find((statestr(:,1)=='N' & statestr(:,2)=='R')); r=find((statestr(:,1)=='R' & statestr(:,2)~='a'));
        w=find(statestr(:,1)=='W' & statestr(:,2)~='a');
        nr2=find((statestr(:,1)=='N' & statestr(:,2)=='a')); r3=find((statestr(:,1)=='R' & statestr(:,2)=='a'));
        w1=find(statestr(:,1)=='W' & statestr(:,2)=='a');
        mt=find(statestr(:,1)=='M' & statestr(:,2)~='a');
        ma=find(statestr(:,1)=='M' & statestr(:,2)=='a');

        nr(nr>maxep)=[];nr2(nr2>maxep)=[];r(r>maxep)=[];r3(r3>maxep)=[];w(w>maxep)=[];w1(w1>maxep)=[];mt(mt>maxep)=[];ma(ma>maxep)=[];

        ww1mt=sort([w;w1;ma;mt]);wake=zermat; wake(ww1mt)=1;[bastend badur]=BriefAwakenings(wake,maxep);
        ba=zermat; for b=1:length(badur) ba(bastend(b,1):bastend(b,2))=1; end; mt=find(ba);
        [x,y]=intersect(w,mt); w(y)=[];[x1,y1]=intersect(w1,mt); w1(y1)=[]; mt=mt';

        deri=ders(d,:);

%         fn=[mouse,'-',day,'-',deri,'-VSspec']
        fn=[mouse,'-',day,'-',deri]
       
        eval(['save ',pathout,fn,'.mat spectr w nr r w1 nr2 r3 mt ma bastend -mat']);

        if d==1 spectr1=spectr; elseif d==2 spectr2=spectr; elseif d==3 spectr3=spectr; else spectr4=spectr; end

        spN=mean(spectr(nr,:));
        spW=mean(spectr(w,:));
        spR=mean(spectr(r,:));

        sp=[spW;spN;spR];
        sp(:,1:2)=NaN;

        
        co = [0 .27 .8;    
              .2 .6 0;
              .86 .1 0;
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840];
            set(groot,'defaultAxesColorOrder',co)



        subplot(1,2,d) %subplot(1,3,d)

        semilogy(f,sp,'LineWidth',2)
        grid on
        legend('W','N','R')
        title (fn)
        ylim([5*10^0 10^4])

    end
%     saveas(gcf,[pathfig,'\',mouse,'_',day,'_VSspec'],'tiff')
            %pause
%             close all
end




% pathin=[path,'fft\']
% pathout=[path,'outputVS\']; mkdir(pathout)







% %%%%%%%%%% all baseline data ChR2

% pathin=[path,'fft30Hz\Baseline-24h-fft\']

% maxep=21600; 

% mousenames=strvcat('GDCh1','GDCh4','GDCh5','GDCh6','GDCh7','GDCh8','GDCh9','GDCh12','GDCh15','GDCh16','GDCh17');
% days=strvcat('130218','020418','020418','290418','120518','100618','100618','100618','181218','181218','181218');
% mousenames=strvcat('GDCh7');
% days=strvcat('300618');
% mousenames=strvcat('GFP7');
% days=strvcat('081018');
% mousenames=strvcat('GDCh8');
% days=strvcat('150618');
% mousenames=strvcat('GFP4');
% days=strvcat('260618');
% mousenames=strvcat('GFP5');
% days=strvcat('051018');
% mousenames=strvcat('GDCh1','GDCh1','GDCh4','GDCh4','GDCh5','GDCh5','GDCh6','GDCh6','GDCh8','GDCh8','GDCh9','GDCh9','GDCh15','GDCh15','GDCh16','GDCh16','GDCh17','GDCh17');
% days=strvcat('180218','190218','300318','310318','300318','310318','140418','150418','060618','070618','060618','070618','201218','211218','201218','211218','201218','211218');

% %%%%%%%%%% all 24h stim data ChR2
% mousenames=strvcat('GDCh1','GDCh4','GDCh5','GDCh6','GDCh7','GDCh8','GDCh9','GDCh12','GDCh15','GDCh16','GDCh17');
% days=strvcat('140218','030418','030418','120418','130518','050618','150618','090618','191218','191218','191218');

% %%%%%%%%%% all baseline data GFP
% mousenames=strvcat('GFP1','GFP2','GFP3','GFP4','GFP5','GFP6','GFP7','GFP8')
% days=strvcat('250218','290418','300618','300618','031018','081018','081018','080119')
% %%%%%%%%%% all 24h stim data GFP
% mousenames=strvcat('GFP1','GFP2','GFP3','GFP4','GFP5','GFP6','GFP7','GFP8')
% days=strvcat('260218','120418','260618','040718','061018','061018','061018','070119')



%%%%%%%% Wake Enhancement ChR2

% pathin=[path,'fft30Hz\WakeEnhance\']

% maxep=10800; 

% mousenames=strvcat('GDCh1','GDCh1','GDCh4','GDCh4','GDCh5','GDCh5','GDCh6','GDCh6','GDCh7','GDCh7','GDCh8','GDCh8','GDCh9','GDCh9','GDCh12','GDCh12','GDCh16','GDCh16')
% days=strvcat('040518','060518','160518','140518','040518','060518','300418','020518','030718','010718','130618','110618','130618','110618','130618','110618','090119','110119')
% mousenames=strvcat('GDCh17','GDCh17') %%occ
% days=strvcat('090119','110119') %%occ
% der=strvcat('occ_30Hz'); 
% ders=strvcat('occ')



%%%%%%%%%% Wake Enhancement GFP

% pathin=[path,'fft30Hz\WakeEnhance\']

% maxep=10800; 

% mousenames=strvcat('GFP1','GFP1','GFP2','GFP2','GFP3','GFP3','GFP4','GFP4','GFP5','GFP5','GFP6','GFP6','GFP7','GFP7','GFP8','GFP8')
% days=strvcat('140518','160518','020518','300418','030718','010718','010718','030718','041018','011018','111018','091018','111018','091018','090119','110119')
 


%%%%%%%%%% 1, 2, 5 Hz stim 

% pathin=[path,'fft30Hz\1Hz2Hz5Hz\']
% maxep=21600; 

% mousenames=strvcat('GDCh15','GDCh15','GDCh15','GDCh16','GDCh16','GDCh16','GDCh17','GDCh17','GDCh17');
% days=strvcat('241218','251218','261218','241218','251218','261218','241218','251218','261218');



%%%%%%%%%% 2h, 1Hz

% pathin=[path,'fft30Hz\2h1Hz\']
% maxep=21600; 

% mousenames=strvcat('GDCh15','GDCh16','GDCh17');
% days=strvcat('231218','231218','231218');




%%%%%%%%%% Sedation - open ReadVSspecEEG_sedation





%%%%%%%%%% 8sec

% pathin=[path,'fft30Hz\8sec\']
% maxep=21600; 

% mousenames=strvcat('GDCh16','GDCh17');
% days=strvcat('060119','060119');
% 
