
%%%%%%%%%%%%%%%
clear all;
close all;
outp=[];nums=[];

% pathout='I:\OptoMod\';
pathout='E:\OptoInh\';  %mkdir(pathout);
pathsig=[pathout,'OutputSignals\'];  %mkdir(pathsig);




% pathtanks=['I:\rawdata\'];  
% tank='May2018yellow' 
% event='estm'  %'est1,2,3 or 4' 
% blocks=char('GDCh7_8_120518','GDCh7_8_130518','GDCh7_8_140518','GDCh7_8_150518','GDCh7_8_160518');
% recorddates=char('120518','130518','140518','150518','160518');
% mousenames=char('GDCh7','GDCh8');
% days=[4:5];
% mice=[2];

% pathtanks=['F:\rawdata\'];  
% tank='August2019' 
% event='estm'  %'est1,2,3 or 4' 
% blocks=char('GDCh18_19_20_21_170919');
% recorddates=char('170919');
% mousenames=char('GDCh18','GDCh19','GDCh20','GDCh21');
% days=[1];
% mice=[3];
% event=['est',int2str(mice)];
% 
% pathtanks=['G:\rawdata\'];       
% tank='September2018b';
% event='est2'
% blocks=char('GFP6_7_5_111018');%
% recorddates=char('111018');
% mousenames=char('GFP6','GFP7','GFP5','GFP5');
% days=[1];
% mice=[3];

% pathtanks=['F:\rawdata\'];  
% tank='August2019' 
% event='estm'  %'est1,2,3 or 4' 
% blocks=char('GDCh18_19_20_21_110819');
% recorddates=char('080819','110819');
% mousenames=char('GDCh18','GDCh19','GDCh20','GDCh21');
% days=[1:2];
% mice=[1:4];

% pathtanks=['G:\rawdata\'];  
% tank='April2021pink';
% event='est4';  %'est1,2 or 4' %est1:SD, est4:24h, est2: blueLED 1h contstim
% blockhead='Gf1_Ar1_2_3_Gf2_Ar4_7';
% blockhead2='2021';
% blocks=char('0408','0412','0415','0416');
% recorddatehead='2021';
% recorddates=char('080421','120421','150421','160421');
% days=[4];
% mousenames=char('Gf1','Ar1','Ar2','Ar3','Gf2','Ar4','Ar7');
% mice=[1];%[1:numanim];


pathtanks=['E:\rawdata\']; 
tank='April2021p';  %April2021pink  %May2021y %May2021b  %May2021y
event='est2';  %'est1,2 or 4' %est1:SD, est4:24h, est2: blueLED 1h contstim
blockhead='Gf1_Ar1_2_3_Gf2_Ar4'; %Gf1_Ar1_2_3_Gf2_Ar4_7 %Ar8_Gf4_5_Ar9 %Arch5_Gf3 %Gf1_Ar1_2_3_Gf2_Ar4
blockhead2='2021';
recorddatehead='2021';
blocks=char('0406','0407','0408','0411','0412','0415','0416','0417');
recorddates=char('060421','070421','080421','110421','120421','150421','160421','170421');
% blocks=char('0506','0510','0514','0515','0517','0518');
% recorddates=char('060521','100521','140521','150521','170521','180521');

days=[8];
mousenames=char('Gf1','Ar1','Ar2','Ar3','Gf2','Ar4','Ar7'); numanim=size(mousenames,1);
% mousenames=char('Ar8','Gf4','Gf5','Ar9'); numanim=size(mousenames,1);
% mousenames=char('Ar5','Gf3'); numanim=size(mousenames,1);
mice=[2];%[1:numanim];

numanim=size(mousenames,1);%1;
numch=1;


num_chunks=4; % for resampling
tail_length=20000;
visualize=0;

for mouse=mice
    mousename=mousenames(mouse,:); mousename(isspace(mousename))=[];

    for ii=days
        block=[blockhead,'_',blockhead2,blocks(ii,:)]; block(isspace(block))=[];                            
        recorddate=recorddates(ii,:); recorddate(isspace(recorddate))=[];

        maxret = 1000000; maxevents = 100000; step=10000;

        TTX = actxcontrol('TTank.X');
        % Then connect to a server.
        invoke(TTX,'ConnectServer', 'Local', 'Me')
        invoke(TTX,'OpenTank', [pathtanks tank], 'R')

        % Select the block to access
        invoke(TTX,'SelectBlock', block)
        invoke(TTX, 'ResetGlobals')
        invoke(TTX, 'SetGlobalV', 'MaxReturn', maxret)

        L = invoke(TTX, 'ReadEventsV', maxret, event, 1, 0, 0, 1, 'ALL');
        SR = invoke(TTX, 'ParseEvInfoV', 0, 1, 9)

        f1=SR;
        f2=256;
        
        t1=0; t2=43200; %t2=43200*2;
        
        ts = [];
        i=1;

        %STIM
        for chan = [1]
%             sig=[];

            %%%%%%%%%%%%%% new!
            invoke(TTX, 'SetGlobalV', 'WavesMemLimit', 1024^3);
            %%%%%%%%%%%%%%%%%%
            
            invoke(TTX, 'SetGlobalV', 'T1', t1);
            invoke(TTX, 'SetGlobalV', 'T2', t2);
            invoke(TTX, 'SetGlobalV', 'Channel', chan);
            y = invoke(TTX,'ReadWavesV',event);

            sig=y';%sig=[sig y'];   
            
            sig1=sig*10^6;
            resampled_sig=resampling_new(sig1,f1,f2,num_chunks,tail_length,visualize); %resampling_rev
            clear sig_or;
            
            %%%%%%%%% Arch - comment out
            %%%%%%%%% ChR2 - comment
            resampled_sig=resampled_sig/1000000; %10000 %10000000 %100000000 %1000 
            %10000000000 for GDCh18-21, 1000000 for Ar1-7, Gf1-2
            
            % make an Xaxis (in seconds)for the resampled signal  
            xAxRSig=0:1/256:(size(resampled_sig,2)/256) -1/256;    

      
            fnout=[mousename,'-',event,'-',recorddate];    
            
            eval(['save ',pathsig,fnout,'.mat event resampled_sig block tank xAxRSig -mat']);

            clear sig resampled_sig xAxRSig;

        end

    end
end

%%%%%%%%% Amplitude Collection - copy and run if necessary
% 
% path='I:\';
% pathout='I:\OptoMod\';
% pathsig=[pathout,'OutputSignals\'];
% fnout2=['GDCh18-estm-080819'];
% load ([pathsig,fnout2,'.mat'],'resampled_sig')
% resampled_sig=resampled_sig/100; %/1000000
% save ([pathsig,fnout2,'.mat'],'resampled_sig')
% clear resampled_sig;




%%%%%%%% log

% blocks=char('GDCh18_19_20_21_080819');
% recorddates=char('080819');
% mousenames=char('GDCh18','','GDCh20','');
% days=[1];
% mice=[1 3];

% blocks=char('GDCh18_19_20_21_110819');
% recorddates=char('110819');
% mousenames=char('','GDCh19','','GDCh21');
% days=[1];
% mice=[2 4];

% blocks=char('GDCh18_19_20_21_140819');
% recorddates=char('140819');
% mousenames=char('GDCh18','GDCh19','GDCh20','GDCh21');
% days=[1];
% mice=[1:4];
% 
% 
