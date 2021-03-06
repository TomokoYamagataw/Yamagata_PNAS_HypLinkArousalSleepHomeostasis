
%% Import data from text file.
% Script for importing data from the following text file:
%
%    I:\optogenetics\Calorimetry\Calo_2020_01_09.CSV
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2020/01/23 09:22:50

%% Initialize variables.
close all;
clear all;

path='E:\Optoinh\'; 
pathout = [path,'Calorimetry\';]; mkdir(pathout)
exte='.dat';
delimiter = ',';



%% Chose file
% mouse='GAD1'
% recorddate='090120'jo
% filename = 'I:\optogenetics\Calorimetry\Calo_2020_01_09.CSV';
% startRow = 499;
% endRow = 9504; %infinit

% mouse='GAD2'
% recorddate='100120'
% filename = 'I:\optogenetics\Calorimetry\Calo_2020_01_10.CSV';
% startRow = 3; %5792; %3
% endRow = Inf; %infinit

%%%% GDVV7-3c GAD3
%%%% GDVV20-4d GAD4
%%%% GOAR7-3d GAD5
%%%% GOAR3-1e GAD6

% pathin = [path,'Calorimetry\Calorimetry_C\']; %Calorimetry_B %Calorimetry
% mouse='GAD3' %GDVV7-3GAD3c_7C_P10_20210411-1820_Dex
% recorddate='110421' %120421
% filename =  [pathin,'Calo_2021_04_11_pre.CSV']; %Calo_2021_04_12
% startRow = 3; %5792; %3
% endRow = Inf; %infinit

% % % % pathin = [path,'Calorimetry\Calorimetry_C\']; %Calorimetry_B %Calorimetry
% % % % mouse='GAD4' %GDVV20-4d_7C_P10_20210407-0408_2356_Dex
% % % % recorddate='100421'
% % % % filename =  [pathin,'Calo_2021_04_10_Dex.CSV'];
% % % % startRow = 3; %5792; %3
% % % % endRow = Inf; %infinit
pathin = [path,'Calorimetry\Calorimetry_C\']; %Calorimetry_B %Calorimetry
mouse='GAD4' %GDVV20-4d_7C_P10_20210403_Dex_Calo2119
recorddate='030421'
filename =  [pathin,'Calo_2021_04_03_Dex.CSV'];
startRow = 3; %5792; %3
endRow = Inf; %infinit

% pathin = [path,'Calorimetry\Calorimetry_B\']; %Calorimetry_B %Calorimetry
% mouse='GAD5' %GOAR7-3d_7B_P10_20210405_Dex_2359
% recorddate='050421'
% filename =  [pathin,'Calo_2021_04_05.CSV'];
% startRow = 3; %5792; %3
% endRow = Inf; %infinit

% pathin = [path,'Calorimetry\Calorimetry_B\']; %Calorimetry_B %Calorimetry
% mouse='GAD6' %GOAR3-1e_7B_P09_20210410-1918_Dex
% recorddate='100421'
% filename =  [pathin,'Calo_2021_04_10_Dex.CSV'];
% startRow = 3; %5792; %3
% endRow = Inf; %infinit


%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');


% %
% for block=2:length(startRow)
%     frewind(fileID);
%     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     for col=1:length(dataArray)
%         dataArray{col} = [dataArray{col};dataArrayBlock{col}];
%     end
% end
% %


%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[3,4,5,6,7,8]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

% Convert the contents of columns with dates to MATLAB datetimes using the
% specified date format.
try
    dates{2} = datetime(dataArray{2}, 'Format', 'HH:mm:ss', 'InputFormat', 'HH:mm:ss');
catch
    try
        % Handle dates surrounded by quotes
        dataArray{2} = cellfun(@(x) x(2:end-1), dataArray{2}, 'UniformOutput', false);
        dates{2} = datetime(dataArray{2}, 'Format', 'HH:mm:ss', 'InputFormat', 'HH:mm:ss');
    catch
        dates{2} = repmat(datetime([NaN NaN NaN]), size(dataArray{2}));
    end
end

dates = dates(:,2);

%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [3,4,5,6,7,8]);
rawStringColumns = string(raw(:, [1,9]));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

%% Create output variable
Calo = table;
Calo.Function = categorical(rawStringColumns(:, 1));
Calo.Time = dates{:, 1};

Calo.TempC = cell2mat(rawNumericColumns(:, 1));
TempC = cell2mat(rawNumericColumns(:, 1));

Calo.FlowLmin = cell2mat(rawNumericColumns(:, 2));
FlowL = cell2mat(rawNumericColumns(:, 2));

Calo.O2mL = cell2mat(rawNumericColumns(:, 3));
O2mL = cell2mat(rawNumericColumns(:, 3));

Calo.CO2mLmin = cell2mat(rawNumericColumns(:, 4));
CO2mL = cell2mat(rawNumericColumns(:, 4));

Calo.RER = cell2mat(rawNumericColumns(:, 5));
RER = cell2mat(rawNumericColumns(:, 5));

Calo.HPmW = cell2mat(rawNumericColumns(:, 6));
HPmW = cell2mat(rawNumericColumns(:, 6));

Calo.Comment = rawStringColumns(:, 2);

%%%%%%%%%% O2 calibration
figure
plot(O2mL)
hold on

% out1=find(O2mL<=0);
O2mLc=O2mL;
for ii=2:length(O2mLc)-20
    if O2mLc(ii,1)<=0
        O2mLc(ii-1:ii+20,1)=NaN; 
    end
    ii=ii+15;
end
out1=isnan(O2mLc);

figure
plot(O2mLc)

%%%%%%%%%% CO2 calibration
figure
plot(O2mL)
hold on

CO2mLc=CO2mL;
CO2mLc(out1)=NaN;

figure
plot(CO2mLc)

%%%%%%%%%% HPmWc calibration
figure
plot(HPmW)
hold on

HPmWc=HPmW;
HPmWc(out1)=NaN;

figure
plot(HPmWc)


%%%%%%%%%% RER calibration
figure
plot(RER)
hold on

RERc=RER;
RERc(out1)=NaN;

figure
plot(RERc)

%%%%%%% remove noisy RER %%%%%
stdRER=nanstd(RERc);
medRER=nanmedian(RERc);
out2=find(RERc>medRER+4*stdRER);
RERc(out2)=NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(RERc)

pause % check RER noise removal working well


fname2=['calo_',mouse,'_',recorddate,'_Dex']; %=Dex %_pre


eval(['save ',pathout,fname2,'.mat Calo TempC FlowL O2mL O2mLc CO2mL CO2mLc RER RERc HPmW HPmWc -mat']);

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% Calo20200109.Time=datenum(Calo20200109.Time);

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp dates blankDates anyBlankDates invalidDates anyInvalidDates rawNumericColumns rawStringColumns R idx;