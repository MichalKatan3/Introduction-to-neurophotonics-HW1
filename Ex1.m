%% Introduction to neurophotonics - Ex. 1
% Michal Katan (206799793)
% Channa Shapira (314762006)

%% Reset
clear all;
close all;
clc;

%% Load files
extinctionCoefficientsFile = readtable('ExtinctionCoefficientsData.csv');
relDPFfile = readtable('RelativeDPFCoefficients.csv');
DPFperTissueFile = readtable('DPFperTissue.txt');

%% Variables the user should choose

prompt = {'File: 1 or 2?','Tisue type: 1 - adult forehead, 2 - baby head, 3 - adut head, 4 - adult leg','Channels:'};
dlgtitle = 'Enter your selected';
fieldsize = [1 50; 1 50; 1 50];
definput = {'1','3','1,2'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

file_selct = str2num(answer{1});
tissueType = str2num(answer{2});
plotChannelIdx = str2num(answer{3});

%% Check validity of inputs
if  isempty(file_selct)
    errordlg('The first field must be an integer! The code will use the default value of 1.') 
    file_selct = 1;
elseif file_selct == 1
    dataFile = load('FN_032_V1_Postdose1_Nback.mat');
elseif file_selct == 2
    dataFile = load('FN_031_V2_Postdose2_Nback.mat');
elseif file_selct ~= 1 && file_selct ~= 2
    errordlg('The first field should be a number of 1 or 2! The code will use the default value of 1.')
    file_selct = 1;
end

% In case tissueType is not 1, 2, 3 or 4.
if isempty(tissueType)
    errordlg('The second field must be an integer! The code will use the default value of 3.')
    tissueType = 3;
elseif tissueType~=1 &&  tissueType~=2 && tissueType~=3 && tissueType~=4
    errordlg('tissueType must be a number between 1-4! The code will use the default value of 3.')
    tissueType = 3;
end

% In case the user choses more than 20 channels or a channel bigger than 20
if isempty(plotChannelIdx)
   errordlg('The third field must be an integer! The code will use the default value of 1 and 2.') 
   plotChannelIdx = [1,2];
elseif length(plotChannelIdx) > 20 || max(plotChannelIdx) > 20
   errordlg('There are only 20 channels! The code will use the default value of 1 and 2.') 
   plotChannelIdx = [1,2];
elseif  min(plotChannelIdx) <= 0 % In case the input is not a positive integer
   errordlg('The third field must be a positive integer! The code will use the default value of 1 and 2.') 
   plotChannelIdx = [1,2];
end

%% Extracting dHbR and dHbO
SDS = 3; %cm
[dHbR , dHbO, fig] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile );


