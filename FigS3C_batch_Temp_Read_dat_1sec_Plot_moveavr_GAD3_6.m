close all;
clear all;


%% Import data from text file.
% Script for importing data from the following text file:
%
%    F:\Optrix_PIX\GFP6_7_191018_0955-201018_0912_23h_Dexmedetomidine.dat
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2019/08/23 23:24:12



%%%%%% Initialize variables.
% filename = 'E:\OptoInh\Calorimetry\Temp\GDVV7-3c_7C_P10_20210412-0941_Dex_postDex_24h.dat'; %GAD34
% camera_batch='GAD3_';
% mouse_batch='GAD3-';
% recorddate='120421';
% fnames=char('GAD3','cold');
% Dst = 12290; % Inj start
% Dst2 = 12546; % Inj end
% Dend = 28584; % Antisedan injection
% 
% filename = 'E:\OptoInh\Calorimetry\Temp\GDVV20-4d_7C_P10_20210403_Dex_Calo2119.dat'; %GAD4
% %%%   GDVV20-4d_7C_P10_20210407-0408_2356_Dex
% %%%   GDVV20-4d_7C_P10_20210403_Dex_Calo2119
% %%%   GDVV20-4d_7C_P10_20210407-0408_2356_Dex
% camera_batch='GAD4_';
% mouse_batch='GAD4-';
% recorddate='030421';
% fnames=char('GAD4','cold');
% Dst = 10371; % Inj start
% Dst2 = 10454; % Inj end
% Dend = 25335; % Antisedan injection
% 
% filename = 'E:\OptoInh\Calorimetry\Temp\GOAR7-3d_7B_P10_20210405_Dex_2359.dat'; %GAD5
% camera_batch='GAD5_';
% mouse_batch='GAD5-';
% recorddate='050421';
% fnames=char('GAD5','cold');
% Dst = 11764; % Inj start
% Dst2 = 12166; % Inj end
% Dend = 25653; % Antisedan injection
% 
% %  
% filename = 'E:\OptoInh\Calorimetry\Temp\GOAR3-1e_7B_P09_20210410-1918_Dex.dat'; %GAD6
% camera_batch='GAD6_';
% mouse_batch='GAD6-';
% recorddate='100421';
% fnames=char('GAD6','cold');
% Dst = 11808; % Inj start
% Dst2 = 12245; % Inj end
% Dend = 26687; % Antisedan injection


% ctemp=22.7;
% numdata=2;


path= 'E:\OptoInh\';
pathin = [path,'Dex_dat\']; %mkdir(pathout)
pathout = [path,'temp\';];%mkdir(pathout)
pathout2 = [path,'Calorimetry\'];
pathfig = [path,'Figure_Temp\']; %mkdir(pathfig)
exte='.dat';
delimiter = '\t';
startRow = 8;
% endRow = 81769;

for n=1  %1:size(recorddates)
%     recorddate=recorddates(n,:);   
%     filename = [pathin,camera_batch,recorddate,exte];

%% Format for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
% GFP671910180955201018091223hDexmedetomidine = table(dataArray{1:end-1}, 'VariableNames', {'Time','VarName2','VarName3','ambient','sm','sm1'});

    %% Allocate imported array to column variable names
    Time = dataArray{:, 1}; 
    L1=length(Time)-2; 
%     s1=size(dataArray,2)-2; %exclude the time and the last blank
    Time = Time(1:L1,1);
   
%      
%      %%%% baseline correction %
%      temp = dataArray{:, 3}; temp = temp(1:L1,1);
%      resampled_temp=temp';
%      out=find(resampled_temp>24);resampled_temp(out)=[];
%      basetemp=nanmean(resampled_temp);
%      corrtemp=ctemp-basetemp;   
       
     %%%%%%%%%%%%%%%% 36sec max and Moveavr       
        ii=1;
        temp = dataArray{:, ii+1}; temp = temp(1:L1,1);
        out0=find(temp>40);temp(out0)=NaN;
%         out2=find(temp<29);temp(out2)=NaN;

        resampled_temp=temp';
%         corrected_temp=resampled_temp+corrtemp;
        fname=fnames(ii,:); fname(isspace(fname))=[];

        maxep=2400; % maxep=10800; %21600
        rectemp=NaN*ones(1,86400);
        Lm=length(resampled_temp);
        rectemp(1,1:Lm)=resampled_temp;
    %         tempZsec=reshape(rectemp,4,21600); 
        tempZsec=reshape(rectemp,36,[]); 
        Temp=max(tempZsec); 
        %%%%%%%%% moving average %%%%%%%%%%
        lag = 50; %lag = 450;
        Tempmoveavr = movmean(Temp,lag,'omitnan');
    %         Tempmoveavr=Tempmoveavr(1:maxep);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot
    x1=1:maxep; x1=x1./100;
    axistemp=[0.01 12.01 20 40]; %axistemp=[0 24 20 35]; 
    options.color_line4 = [255 0 0]./255; % Temp2

    s1=Dst/4/900; s2=Dst2/4/900;
    % Temperature
    figure
    fill([s1 s2 s2 s1], [20 20 40 40],'g')
    hold on
    plot(x1,Temp,'LineWidth',2,'Color',options.color_line4) %Tempmoveavr  %Temp
    axis(axistemp)
    set(gca,'YTick',[20:5:40],'TickDir','out','FontSize',14)% set(gca, 'YTick', []); 
    set(gca, 'XTick', []); set(gca,'XTick',[0:2:12],'TickDir','out')
    ylabel('Temperature[C]','FontSize',16)

    fn0=[fname,'_',recorddate,'_36secMax'];
%         eval(['save ',pathout,fn0,'.mat Temp Tempmoveavr -mat']);
%         saveas(gcf,[pathfig,fn0,'_tmp36'],'tiff')%'_moveavrtemp'
%     saveas(gcf,[pathfig,fn0,'_moveavrtemp'],'tiff')%'_moveavrtemp'



    %%%%%%%% aligned to Dex/atipamezole injection
    T0=Temp';
    p1=300;
    nanmat1=NaN*ones(p1,1);
    nanmat2=NaN*ones(p1,1);
    
    a1=floor(Dst/36); a2=floor(Dst2/36); a3=floor(Dend/36);
    if a1<p1
    T1=nanmat1;
    T1(p1-a1+1:p1,:)=T0(1:a1,1);
    else
    T1=T0(a1-(p1-1):a1);
    end
    out1=find(T1<29);T1(out1)=NaN;
    
    T2=T0(a2:a2+(p1-1));
    
    if Lm<a3+(p1-1)
        T3=nanmat1; 
        T3(1:Lm-a3+1,:)=T0(a3:Lm,1);
    else
     T3=T0(a3:a3+(p1-1));
    end
%     out3=find(T3<29);T3(out3)=NaN;
    
    aPrism_Temp2=[T1;T2;T3];
    fn5=[camera_batch,'DexTempAligned2'];
    eval(['save ',pathout2,fn5,'.mat aPrism_Temp2 -mat']);

    figure
    plot(aPrism_Temp2)
   

    %% Clear temporary variables
    clearvars filename fileID dataArray ans;
end
clearvars delimiter startRow formatSpec ;