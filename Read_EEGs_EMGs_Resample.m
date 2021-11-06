clear all;
close all;
outp=[];nums=[];

epochl=4;
pathout='E:\OptoInh\';  %mkdir(pathout);
pathsig=[pathout,'OutputSignals\'];  %mkdir(pathsig);


pathtanks=['E:\rawdata\']; 
tank='April2021p';  %April2021p  %May2021y %May2021b %May2021y
blockhead='Gf1_Ar1_2_3_Gf2_Ar4';  %Gf1_Ar1_2_3_Gf2_Ar4_7; %'Ar8_Gf4_5_Ar9'; %'Arch5_Gf3' %Gf1_Ar1_2_3_Gf2_Ar4
blockhead2='2021';

blocks=char('0406','0407','0408','0411','0412','0415','0416','0417','0417D4');
recorddatehead='2021';
recorddates=char('060421','070421','080421','110421','120421','150421','160421','170421','170421D4');
% blocks=char('0505','0506','0509','0510','0513','0514','0515','0517','0518');
% recorddatehead='2021';
% recorddates=char('050521','060521','090521','100521','130521','140521','150521','170521','180521');
days=[9];
mousenames=char('Gf1','Ar1','Ar2','Ar3','Gf2','Ar4','Ar7'); numanim=size(mousenames,1);
% mousenames=char('Ar8','Gf4','Gf5','Ar9'); numanim=size(mousenames,1);
% mousenames=char('Ar5','Gf3'); numanim=size(mousenames,1);
mice=[1 2 3 5];%[1:numanim];
% channames=char('fro','occ','foc','emg'); chans=[1:4]; 
channames=char('fro','occ','emg'); chans=[1:3]; 
events1234=char('EEG');


num_chunks=4; % for resampling
tail_length=20000; % don't know what this does
visualize=0;


for mouse=mice % 1:4
    event=[events1234,num2str(mouse)]; event(isspace(event))=[];
    mousename=mousenames(mouse,:); mousename(isspace(mousename))=[];

    for ii=days
        block=[blockhead,'_',blockhead2,blocks(ii,:)]; block(isspace(block))=[];                            
        recorddate=recorddates(ii,:); recorddate(isspace(recorddate))=[];

        maxret = 10000000;  % maxret = 1000000; 
        TTX = actxcontrol('TTank.X');
        % Then connect to a server.
        invoke(TTX,'ConnectServer', 'Local', 'Me');
        invoke(TTX,'OpenTank', [pathtanks tank], 'R');

        % Select the block to access
        invoke(TTX,'SelectBlock', block);
        
        %set parameters to read the events'
        invoke(TTX, 'ResetGlobals')
        invoke(TTX, 'SetGlobalV', 'MaxReturn', maxret);
        
        L = invoke(TTX, 'ReadEventsV', maxret, event, 1, 0, 0, 1, 'ALL');
        SR = invoke(TTX, 'ParseEvInfoV', 0, 1, 9);
        
        f1=SR;
        f2=256;

        % 21600 epochs 24h
%         t1=0; t2=43200*2;
        % 10800 12hr        
        t1=0; t2=43200;
        
        ts = [];
        i=1;

        %EEGs, EMG
        for chan = chans
            %%%%%%% for May2021yellow and May2021blue %%%
%             if chan==3 
%                 continue 
%             end 
            
            sig=[];
          
            %%%%%%%%%%%%%% new!
            invoke(TTX, 'SetGlobalV', 'WavesMemLimit', 1024^3);
            %%%%%%%%%%%%%%%%%%
            
            invoke(TTX, 'SetGlobalV', 'T1', t1);
            invoke(TTX, 'SetGlobalV', 'T2', t2);
            invoke(TTX, 'SetGlobalV', 'Channel', chan);
            y = invoke(TTX,'ReadWavesV',event);

            sig=[sig y']; 
            
            sig1=sig*10^6;
            resampled_sig=resampling_new(sig1,f1,f2,num_chunks,tail_length,visualize);%resampling_rev
            clear sig_or;

            channame=channames(chan,:);            
            fnout=[mousename,'-',event,'-',recorddate,'-',num2str(channame)]    
            %%%%%% If MUA and LFP            
            %             channame=channames(chan,:);            
            %             fnout=[mousename,'-',event,'-',recorddate,'-',num2str(channame)];            

            eval(['save ',pathsig,fnout,'.mat event resampled_sig block tank -mat']);

            clear sig sig1 resampled_sig y;

        end
    end
end

% log
% recorddates=char('080919','060919','050919','020919','030919','010919','310819','300819','280819','250819','220819','210819','200819','190819','180819','170819','160819','150819','140819','130819','120819','110819','100819','090819','080819','070819');
% blocks=char('GDCh18_19_20_21_080919','GDCh18_19_20_21_060919','GDCh18_19_20_21_050919','GDCh18_19_20_21_020919','GDCh18_19_20_21_030919','GDCh18_19_20_21_040919','GDCh18_19_20_21_010919','GDCh18_19_20_21_310819','GDCh18_19_20_21_300819','GDCh18_19_20_21_280819','GDCh18_19_20_21_250819','GDCh18_19_20_21_220819','GDCh18_19_20_21_210819','GDCh18_19_20_21_200819','GDCh18_19_20_21_190819','GDCh18_19_20_21_180819','GDCh18_19_20_21_170819','GDCh18_19_20_21_160819','GDCh18_19_20_21_150819','GDCh18_19_20_21_140819','GDCh18_19_20_21_130819','GDCh18_19_20_21_120819','GDCh18_19_20_21_110819','GDCh18_19_20_21_100819','GDCh18_19_20_21_090819','GDCh18_19_20_21_080819','GDCh18_19_20_21_070819');

% %%%%%%%%%%%%%%%%%%
% pathtanks=['G:\rawdata\'];       
% tank='September2018b';
% blocks=char('GFP6_7_5_111018');%
% days=[1];
% recorddates=char('111018');
% mousenames=char('GFP6','GFP7','GFP5','GFP5'); numanim=size(mousenames,1);
% mice=[4];%[2 3 4 5];%[1:numanim];
% channames=char('fro','occ','emg'); chans=[1:3];
% events1234=char('EEG');

% %%%%%%%%%%%%%%%%%% with bug
% channames=char('vfr','voc','emg'); chans=[1:3];
% events1234=char('vis');

%%%%%%%%%%%%%%%% no bug
% pathtanks='I:\rawdata\'; 
% tank='May2018yellow' 
% event='estm'  %'est1,2,3 or 4' 
% blocks=char('GDCh7_8_120518','GDCh7_8_130518','GDCh7_8_140518','GDCh7_8_150518','GDCh7_8_160518');
% recorddates=char('120518','130518','140518','150518','160518');
% mousenames=char('GDCh7','GDCh8');
% days=[4:5];
% mice=[2];

% pathtanks=['F:\rawdata\']; 
% tank='August2019'; 
% 
% % blocks=char('GDCh18_19_20_21_170919','GDCh18_19_20_21_190919','GDCh18_19_20_21_180919','GDCh18_19_20_21_220919','GDCh18_19_20_21_200919','GDCh18_19_20_21_160919','GDCh18_19_20_21_140919','GDCh18_19_20_21_120919','GDCh18_19_20_21_110919','GDCh18_19_20_21_100919');
% % recorddates=char('170919','190919','180919','220919','200919','160919','140919','120919','110919','100919');
% blocks=char('GDCh18_19_20_21_270819','GDCh18_19_20_21_020919','GDCh18_19_20_21_030919','GDCh18_19_20_21_070919','GDCh18_19_20_21_150919');
% recorddates=char('270819','020919','030919','070919','150919');
% days=[1];
% mousenames=char('GDCh18','GDCh19','GDCh20','GDCh21'); numanim=size(mousenames,1);
% mice=[3];%[2 3 4 5];%[1:numanim];
% channames=char('fro','occ','emg'); chans=[1:3];
% events1234=char('EEG');
% 
