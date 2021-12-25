function [Ca] = getCaFromPtScan(filename, rationize, point, channel)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

st = 7; en = 35;
dtInterp = 0.01;

%% Default parameters

if ~exist('rationize', 'var')
    rationize = false;
end

if ~exist('point', 'var')
    point = 0;
end

if ~exist('channel', 'var')
    channel = 0;
end

%% Extracting the informations
%Fluorescence
imageArray = imread([filename '-point' num2str(point) ' - Channel' ...
    num2str(channel) '.png']);
fluoArray = mean(imageArray(:, st:en), 2);

if rationize
    ratioingImage = imread([filename '-point' num2str(point) ...
        ' - Channel' num2str(~channel) '.png']);
    fluoArray = fluoArray ./mean(ratioingImage(:, st:en), 2);  
end

%Time info
ini = split(fileread([filename '.ini']), {newline, char(13)});
ini = ini(~cellfun('isempty', ini));

nbSamplePerDecay = split(ini{4}, '=');
nbSamplePerDecay = str2num(nbSamplePerDecay{2});

acqRate = split(ini{19}, '=');
acqRate = str2num(acqRate{2})*1000;

interPts = split(ini{13}, '=');
interPts = str2num(interPts{2})/1000;

dt = nbSamplePerDecay/acqRate + interPts;

%% Treating the signal
%SGolayFilt
Ca(1, :) = dt:dt:length(fluoArray)*dt;
Ca(2, :) = fluoArray;

end

