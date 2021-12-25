function cf = zscoreData(structure)
% ZSCOREDATA Turns the raw data acquired in the awake animal into Z-score.
%
% function cf = zscoreData(structure)
%
%   Author: Ali-Kemal Aydin, PhD student & Camille Verdier, Master's intern
%   Mail: ali-kemal.aydin@inserm.fr
%   Affiliations: 
%       * INSERM, CNRS, Institut de la Vision, Sorbonne Université, Paris, France
%   License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
%       See LICENSE.txt or <a href="matlab:web('https://creativecommons.org/licenses/by/4.0/')">here</a>
%       for a human-readable version.
%
%   DESCRIPTION: Turns the raw data acquired an analyzed by Camille into
%   Z-Score. Adapted to the data structure and to the baseline we use.
%__________________________________________________________________________
%   PARAMETERS:
%       structure (struct): structure containing each acquisitions, sorted
%       by date, capillaries and so on.
%__________________________________________________________________________
%   RETURN:
%       cf (struct): structure organized as the input, with data turned
%       into Z-scores.


f1 = fieldnames(structure); %contient les souris
cf = structure;

for k = 1:length(f1)
    f = fieldnames(structure.(f1{k}));
    for j = 1:length(f) % capi
        f2 = fieldnames(structure.(f1{k}).(f{j})); %contient les jours
        for i = 1:length(f2)
            f3 = fieldnames(structure.(f1{k}).(f{j}).(f2{i})); %contient les concentrations
            for a = 1:length(f3)
                f4 = fieldnames(structure.(f1{k}).(f{j}).(f2{i}).(f3{a})); %contient les acquisitions
                for b = 1:length(f4)
                    time = cf.(f1{k}).(f{j}).(f2{i}).(f3{a}).(f4{b})(:,1);
                    data = cf.(f1{k}).(f{j}).(f2{i}).(f3{a}).(f4{b})(:,2);
                    
                    bsl = data(time < 10);
                    data = (data - mean(bsl))/mean(std(bsl));
                    cf.(f1{k}).(f{j}).(f2{i}).(f3{a}).(f4{b})(:,2) = data;
                end
            end
        end
    end
end

end

