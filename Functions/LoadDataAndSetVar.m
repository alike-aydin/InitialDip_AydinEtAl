%% LoadDataAndSetVar.m
% LOADDATAANDSETVAR Load dataset and set some useful variables for the scripts
%
%   Author : Ali-Kemal Aydin, PhD student
%   Date : July 4th, 2021
%   Mail: ali-kemal.aydin@inserm.fr
%   Affiliation : U968, Institut de la Vision, Paris
%   License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
%       See LICENSE.txt or <a href="matlab:web('https://creativecommons.org/licenses/by/4.0/')">here</a> 
%       for a human-readable version.
%
%__________________________________________________________________________

colors = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; ...
    [0.9290, 0.6940, 0.1250]; [0.4940, 0.1840, 0.5560]; ...
    [0.4660, 0.6740, 0.1880]; [0.3010, 0.7450, 0.9330]; ...
    [0.6350, 0.0780, 0.1840]};
colorStim = [0.7 0.7 0.7];

load('data.mat');
d = data;