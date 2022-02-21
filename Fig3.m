function Fig3()
% FIG3 Generating Figure 3 from Aydin et al.
%
% function Fig3() = []
%
%   Author : Ali-Kemal Aydin, PhD student
%   Date : April 19th, 2021
%   Mail: ali-kemal.aydin@inserm.fr
%   Affiliation : U968, Institut de la Vision, Paris
%   License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
%       See LICENSE.txt or <a href="matlab:web('https://creativecommons.org/licenses/by/4.0/')">here</a> 
%       for a human-readable version.
%
%   DESCRIPTION : Generates the panels from Figure 3 in Aydin et al.
%
%__________________________________________________________________________
% PARAMETERS:
%
%__________________________________________________________________________
% RETURN:
%
%__________________________________________________________________________

%% Loading dataset, to be run before any other individual panel
LoadDataAndSetVar;

%% Data treatment
m = {'M1951_110321'; 'M1873_080321'; 'M1854_160321'; 'M1952_070421'};
g = {'Glom1'; 'Glom2'; 'Glom3'; 'Glom4'; 'Glom5'; 'Glom6'};

typeCa = 'Delta';
typeRBC = 'Relative';
typePO2 = 'Delta';

Ca = struct();
RBC = struct();
PO2 = struct();
avgsZScore = [];
xL = [0 29];

tIndexCa = d.M1951_110321.(typeCa).Glom1.Ca(1, :) >= xL(1) ...
    & d.M1951_110321.(typeCa).Glom1.Ca(1, :) <= xL(2);
timeCa = d.M1951_110321.(typeCa).Glom1.Ca(1, tIndexCa);

tIndexVasc = d.M1951_110321.(typePO2).Glom1.RBC(1, :) >= xL(1) ...
    & d.M1951_110321.(typePO2).Glom1.RBC(1, :) <= xL(2);
timeVasc = d.M1951_110321.(typePO2).Glom1.RBC(1, tIndexVasc);

for i=1:length(m)
    Ca.(m{i}) = [];
    RBC.(m{i}) = [];
    PO2.(m{i}) = [];
    for j=1:length(g)
       if isfield(d.(m{i}).(typeCa), g{j})
           Ca.(m{i})(:, j) = sgolayfilt(d.(m{i}).(typeCa).(g{j}).Ca(end-1, tIndexCa), 3, 9);
           RBC.(m{i})(:, j) = d.(m{i}).(typeRBC).(g{j}).RBC(end-1, tIndexVasc);
           PO2.(m{i})(:, j) = d.(m{i}).(typePO2).(g{j}).PO2All(end-1, tIndexVasc);
           avgsZScore(end+1, :) = d.(m{i}).ZScore.(g{j}).PO2All(end-1, tIndexVasc);
       end
    end
end

%% Panel A
% This panel is only an image and is not reproduced here.

%% Panel B

mShow = 1;
xL = [0 29];

figure;
subplot(131); hold on;
yL = [-0.5 1.2];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
for i=1:size(Ca.(m{mShow}), 2)
    plot(timeCa, Ca.(m{mShow})(:, i)/max(Ca.(m{mShow})(:, i)), ...
        'Color', colors{i}, 'LineWidth', 2);
end
xlabel('Time (s)'); ylabel('Normalized Ca^{2+} (a.u.)'); 
xlim(xL); ylim(yL)
set(gca, 'FontSize', 13);

subplot(132); hold on;
yL = [-100 250];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
for i=1:size(RBC.(m{mShow}), 2)
    plot(timeVasc, RBC.(m{mShow})(:, i)*100, 'Color', colors{i}, 'LineWidth', 2);
end
xlabel('Time (s)'); ylabel('\Delta RBC Velocity (%)'); 
xlim(xL); ylim(yL)
set(gca, 'FontSize', 13);

subplot(133); hold on;
yL = [-25 45];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
for i=1:size(PO2.(m{mShow}), 2)
    plot(timeVasc, PO2.(m{mShow})(:, i), 'Color', colors{i}, 'LineWidth', 2);
end
xlabel('Time (s)'); ylabel('\Delta pO_{2} Mean (mmHg)'); 
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

suptitle('FIGURE 5B');

%% Panel C

yL = [-4 8];
bslFrame = [1 10];
dipFrame = [10.5 11.5];
nboot = 10000;
alpha = 0.05;
rng('default');

figure; subplot(211); hold on;
patch('XData', [dipFrame(1) dipFrame(2) dipFrame(2) dipFrame(1)], ...
    'YData', [yL(1) yL(1) yL(2) yL(2)], 'FaceColor', ...
    colors{5}, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', ...
    'Initial Dip Period');
patch('XData', [bslFrame(1) bslFrame(2) bslFrame(2) bslFrame(1)], ...
    'YData', [yL(1) yL(1) yL(2) yL(2)], 'FaceColor', ...
    colors{3}, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', ...
    'Baseline Period');
legend('AutoUpdate', 'off');
plot(timeVasc, avgsZScore(1:end-1, :), 'Color', [0.5 0.5 0.5], ...
    'LineWidth', 0.5);
legend('AutoUpdate', 'on');
plot(timeVasc, avgsZScore(end, :), 'Color', [0.5 0.5 0.5], ...
    'LineWidth', 0.5, 'DisplayName', 'Single Glomerulus');
plot(timeVasc, mean(avgsZScore, 1), 'Color', 'k', 'LineWidth', 2, ...
    'DisplayName', 'Average');
xlabel('Time (s)'); ylabel('\Delta pO_{2} Mean (SD)');
xticks([0 10 20]); yticks([-4:2:8]);
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);



bt_bsl = matrixBootstrap(nboot, @mean, ...
    avgsZScore(:, timeVasc > bslFrame(1) & timeVasc < bslFrame(2)));
m = mean(mean(avgsZScore(:, timeVasc > dipFrame(1) & timeVasc < dipFrame(2)), 2));
thresh = prctile(bt_bsl(:, end), alpha/2*100);
threshUp = prctile(bt_bsl(:, end), 100*(1-alpha/2));
pval = sum(bt_bsl(:, end) < m)/nboot;
disp('### Figure 2A ###');
disp(['No dip in the average PO2 response with p = ' ...
    num2str(pval) '.']);

subplot(212); hold on;
patch('XData', [thresh threshUp threshUp thresh], 'YData', ...
    [0 0 1000 1000], 'FaceColor', colorStim, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.5, 'DisplayName', 'CI 95%');
line([m m], ylim, 'Color', colors{5}, 'LineWidth', 3, 'DisplayName', ...
    'Initial Dip Period');
histogram(bt_bsl(:, end), [-1:0.01:1], 'DisplayName', ...
    'Baseline Distribution', 'LineStyle', 'none', 'FaceColor', colors{3});
legend();
xticks([-0.2 0 0.2]); yticks([0 500 1000]);
xlim([-0.25 0.25]); ylim([0 1000]);
ylabel('Count'); xlabel('Average Z-score (SD)');
set(gca, 'FontSize', 13);


end

