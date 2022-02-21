function Fig2()
% FIG2 Generating Figure 2 from Aydin et al.
%
% function Fig2() = []
%
%   Author : Ali-Kemal Aydin, PhD student
%   Date : January 5th, 2021
%   Mail: ali-kemal.aydin@inserm.fr
%   Affiliation : U968, Institut de la Vision, Paris
%   License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
%       See LICENSE.txt or <a href="matlab:web('https://creativecommons.org/licenses/by/4.0/')">here</a> 
%       for a human-readable version.
%
%   DESCRIPTION : Generates the panels from Figure 2 in Aydin et al.
%
%__________________________________________________________________________
% PARAMETERS:
%
%__________________________________________________________________________
% RETURN:
%
%__________________________________________________________________________

%% Loading dataset, to be run before any other individual panel
LoadDataAndSetVar

%% Data Treatment
% More details on the bootstrapping in the Methods of the paper.

type = 'ZScore';
xL_TS = [0 29]; yL_TS = [-3 15]; xL_BS = [-1.1 1.1]; yL_BS = [0 1000];


tIndexPO2 = d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, :) >= xL_TS(1) ...
    & d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, :) <= xL_TS(2);
timePO2 = d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2);


exps = [d.M1393_021020.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1509_261020.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1511_211020.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    mean([d.M1514_051120.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1514_191120.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2)], 1); ...
    d.M1875_160221.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1512_281020.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1594_261120.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1804_080221.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1854_090321.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1951_180321.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1873_220321.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    mean([d.M1952_240321.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1952_310321.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2)], 1); ...
    mean([d.M1954_240321.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1954_010421.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2)], 1)];

all = {d.M1393_021020.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1511_211020.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1509_261020.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1512_281020.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    [d.M1514_051120.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1514_191120.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2)]; ...
    d.M1594_261120.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1804_080221.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1875_160221.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1854_090321.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1951_180321.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1873_220321.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    [d.M1952_240321.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1952_310321.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2)]; ...
    [d.M1954_240321.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1954_010421.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2)]; ...
    exps}; % This last line is for computing the bootstrap distribution on all mice at once.
 
bslFrame = [1 10];
dipFrame = [10.5 11.5];
nboot = 10000;
bt_bsl = [];
rng('default'); % For reproducibility

% Generation of the bootstrap distribution for all mice and the average
for i=1:length(all)
    mat = all{i};
    bt_bsl(:, i) = matrixBootstrap(nboot, @mean, mat(:, timePO2 > bslFrame(1) & timePO2 < bslFrame(2)));
end


alpha = 0.05;

figure;
subplot(211); hold on;
m = mean(mean(all{end}(:, timePO2 > dipFrame(1) & timePO2 < dipFrame(2)), 2));
thresh = prctile(bt_bsl(:, end), alpha/2*100);
threshUp = prctile(bt_bsl(:, end), 100*(1-alpha/2));
patch('XData', [thresh threshUp threshUp thresh], 'YData', ...
    [0 0 1000 1000], 'FaceColor', colorStim, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.5, 'DisplayName', 'CI 95%');
line([m m], ylim, 'Color', colors{5}, 'LineWidth', 3, 'DisplayName', ...
    'Initial Dip Period');
histogram(bt_bsl(:, end), [-1:0.01:1], 'DisplayName', ...
    'Baseline Distribution', 'LineStyle', 'none', 'FaceColor', colors{3});
legend();
xlim(xL_BS); ylim(yL_BS);
ylabel('Count'); xlabel('Average Z-score (SD)');
set(gca, 'FontSize', 13);

pval = sum(bt_bsl(:, end) < m)/nboot;
disp('### Figure 2A ###');
disp(['Significant dip in the average PO2 response with p = ' ...
    num2str(pval) '.']);

subplot(212); hold on;
patch('XData', [dipFrame(1) dipFrame(2) dipFrame(2) dipFrame(1)], ...
    'YData', [yL_TS(1) yL_TS(1) yL_TS(2) yL_TS(2)], 'FaceColor', ...
    colors{5}, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', ...
    'Initial Dip Period');
patch('XData', [bslFrame(1) bslFrame(2) bslFrame(2) bslFrame(1)], ...
    'YData', [yL_TS(1) yL_TS(1) yL_TS(2) yL_TS(2)], 'FaceColor', ...
    colors{3}, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', ...
    'Baseline Period');
plot(timePO2, exps, ...
    'Color', [0.5 0.5 0.5], 'DisplayName', 'Single Mouse');
plot(timePO2, mean(exps, 1), ...
    'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Average');
xlim(xL_TS); ylim(yL_TS); 
ylabel('\Delta pO_{2} Mean (SD)'); xlabel('Time (s)');
set(gca, 'FontSize', 13);

suptitle('Figure 2A');
%% PANEL B

% Sorting dipped mice (s) and non-dipped (ns)
ns = [3 13 2 6 8 1 14]; ss = [4 5 7 9 10 11 12];
bt_bsl = [bt_bsl(:, ns) bt_bsl(:, ss)];
all = [all(ns); all(ss)];
bonf = alpha/14; % Bonferroni correction for the 14 mice
yL_BS = [0 600]; xL_BS = [-1.1 1.1];
figure;
disp('### Figure 2B ###');

n = 1;
for i=1:3
    for j=1:5
        if n < 15
            if n <= 7
                subplot(4, 4, n); hold on;
            else
                subplot(4, 4, n+2); hold on;
            end
            m = mean(mean(all{n}(:, timePO2 > dipFrame(1) & timePO2 < dipFrame(2)), 2));
            thresh = prctile(bt_bsl(:, n), (bonf/2)*100);
            threshUp = prctile(bt_bsl(:, n), (1-bonf/2)*100);
            patch('XData', [thresh threshUp threshUp thresh], 'YData', ...
                [yL_BS(1) yL_BS(1) yL_BS(2) yL_BS(2)], 'FaceColor', ...
                colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
                'DisplayName', 'CI 95%');
            histogram(bt_bsl(:, n), [-1:0.01:1], 'LineStyle', 'none', ...
                'FaceColor', colors{3});
            line([m m], ylim, 'Color', colors{5}, 'LineWidth', 3 );
            xlim(xL_BS); set(gca, 'FontSize', 13);
            
            pval = sum(bt_bsl(:, n) < m)/nboot;
            if pval < bonf/2
                disp(['Plot n°' num2str(n) ': Significant dip in the ' ...
                    'average PO2 response with p = ' num2str(pval) '.']);
            else
                disp(['Plot n°' num2str(n) ': Non-significant dip in the ' ...
                    'average PO2 response with p = ' num2str(pval) '.']);
            end
            n = n+1;
        end
    end
end

suptitle('FIGURE 2B');

%% PANEL C

flowsDelay = [];
PO2Delay = [];
type = 'Delta';
thresh = 0.2; % 20% threshold for the onset computation
tIndexPO2 = d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, :) >= 1 & d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, :) <= 29;
timePO2 = d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2);

% 3 mice from the previous plot are not in this computation due to the
% inability to correctly measure their onset in PO2 or Flow: M1511_211020, 
% M1804_080221, M1952_240321.

% Computation of the onsets for the flow
flowsDelay(1) = getThreshold_fitSig(squeeze(d.M1393_021020.(type).Oxygen_ET_1p5V_2sec.Flow(1, tIndexPO2)), squeeze(d.M1393_021020.Delta.Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2))', [5 10], [5 15], thresh, false);
flowsDelay(2) = getThreshold_fitSig(squeeze(d.M1514_051120.(type).Oxygen_ET_1p5V_2sec.Flow(1, tIndexPO2)), squeeze(mean([d.M1514_051120.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); d.M1514_191120.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2)], 1))', [5 10], [8 18], thresh, false);
flowsDelay(3) = getThreshold_fitSig(squeeze(d.M1854_090321.(type).Oxygen_ET_1p5V_2sec.Flow(1, tIndexPO2)), squeeze(d.M1854_090321.Delta.Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2))', [5 10], [8 18], thresh, false);
flowsDelay(4) = getThreshold_fitSig(squeeze(d.M1873_220321.(type).Oxygen_ET_1p5V_2sec.Flow(1, tIndexPO2)), squeeze(d.M1873_220321.Delta.Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2))', [5 10], [8 18], thresh, false);
flowsDelay(5) = getThreshold_fitSig(squeeze(d.M1954_240321.(type).Oxygen_ET_1p5V_2sec.Flow(1, tIndexPO2)), squeeze(mean([d.M1954_240321.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); d.M1954_010421.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2)], 1))', [5 10], [8 18], thresh, false);
flowsDelay(6) = getThreshold_fitSig(squeeze(d.M1512_281020.(type).Oxygen_ET_1p5V_2sec.Flow(1, tIndexPO2)), squeeze(d.M1512_281020.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2))', [5 10], [5 18], thresh, false);
flowsDelay(7) = getThreshold_fitSig(squeeze(d.M1511_211020.(type).Oxygen_ET_1p5V_2sec.Flow(1, tIndexPO2)), squeeze(d.M1511_211020.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2))', [5 10], [5 15], 0.2, false);
flowsDelay(8) = getThreshold_fitSig(squeeze(d.M1804_080221.(type).Oxygen_ET_1p5V_2sec.Flow(1, tIndexPO2)), squeeze(d.M1804_080221.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2))', [5 10], [8 18], 0.2, false);
flowsDelay(9) = getThreshold_fitSig(squeeze(d.M1952_240321.(type).Oxygen_ET_1p5V_2sec.Flow(1, tIndexPO2)), squeeze(d.M1952_240321.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2))', [5 10], [8 18], 0.2, false);

% Computation of the onsets for the PO2
PO2Delay(1) = getThreshold_fitSig(squeeze(d.M1393_021020.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2)), squeeze(d.M1393_021020.Delta.Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2))', [5 10], [5 15], thresh, false);
PO2Delay(2) = getThreshold_fitSig(squeeze(d.M1514_051120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2)), squeeze(mean([d.M1514_051120.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); d.M1514_191120.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2)], 1))', [5 10], [5 15], thresh, false);
PO2Delay(3) = getThreshold_fitSig(squeeze(d.M1854_090321.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2)), squeeze(d.M1854_090321.Delta.Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2))', [5 10], [5 18], thresh, false);
PO2Delay(4) = getThreshold_fitSig(squeeze(d.M1873_220321.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2)), squeeze(d.M1873_220321.Delta.Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2))', [5 10], [5 18], thresh, false);
PO2Delay(5) = getThreshold_fitSig(squeeze(d.M1954_240321.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2)), squeeze(mean([d.M1954_240321.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); d.M1954_010421.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2)], 1))', [5 10], [5 15], thresh, false);
PO2Delay(6) = getThreshold_fitSig(squeeze(d.M1512_281020.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2)), squeeze(d.M1512_281020.Delta.Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2))', [5 10], [5 18], thresh, false);

% Calcium onset is at 10 s, so we substract 10 to normalize.
flowsDelay = flowsDelay - 10;
PO2Delay = PO2Delay - 10;

yL = [-0.6 2]; xL = [5 20];
m = 'M1512_281020'; 

figure;
subplot(221); hold on;
x = d.(m).(type).Oxygen_ET_1p5V_2sec.Flow(1, :);
y = d.(m).(type).Oxygen_ET_1p5V_2sec.Flow(end-1, :);
[tFlow, oFlow] = getThreshold_fitSig(x, y', [5 10], [1 17], thresh, false);

patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'LineWidth', 0.5, ...
    'FaceAlpha', 0.5);
plot(x, y/(max(oFlow.fitresult(x))), 'Color', 'k', 'LineWidth', 2);
plot(x, oFlow.fitresult(x)/(max(oFlow.fitresult(x))), 'Color', ...
    colors{1}, 'LineWidth', 3);
line([x(1) x(end)], [thresh thresh], 'Color', colors{1}, 'LineWidth', 3, ...
    'LineStyle', ':');
line([tFlow tFlow], [-1 thresh], 'Color', colors{1}, 'LineWidth', 3, ...
    'LineStyle', ':');
scatter(tFlow, yL(1), 60, colors{1}, 'v', 'filled');
yticks([0 1 2]);
ylim(yL); xlim(xL);
xlabel('Time (s)'); ylabel('Normalized Flow');
set(gca, 'FontSize', 13);

subplot(223); hold on;
x = d.(m).(type).Oxygen_ET_1p5V_2sec.PO2All(1, :);
y = d.(m).(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, :);
[tPO2, oPO2] = getThreshold_fitSig(x, y', [5 10], [1 20], thresh, false);

patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');
plot(x, y/(max(oPO2.fitresult(x))), 'Color', 'k', 'LineWidth', 2, ...
    'DisplayName', 'Data');
plot(x, oPO2.fitresult(x)/(max(oPO2.fitresult(x))), 'Color', colors{2}, ...
    'LineWidth', 3, 'DisplayName', 'Fit');
legend('AutoUpdate', 'off');
line([x(1) x(end)], [thresh thresh], 'Color', colors{2}, 'LineWidth', 3, ...
    'LineStyle', ':');
line([tPO2 tPO2], [-1 thresh], 'Color', colors{2}, 'LineWidth', 3, ...
    'LineStyle', ':');
scatter(tPO2, yL(1), 60, colors{2}, 'v', 'filled');
yticks([0 1 2]);
ylim(yL); xlim(xL); legend();
xlabel('Time (s)'); ylabel('Normalized pO_{2}');
set(gca, 'FontSize', 13);

subplot(222); hold on;
b = bar(1, mean(flowsDelay(1:6)), 'FaceColor', 'flat');
b.CData = colors{1};
c = bar(2, mean(PO2Delay), 'FaceColor', 'flat');
c.CData = colors{2};
scatter(repmat(1, 1, length(PO2Delay)), flowsDelay(1:6), 40, 'k', ...
    'filled');
scatter(repmat(2, 1, length(PO2Delay)), PO2Delay, 40, 'k', 'filled');
h = line([repmat(1, 1, length(PO2Delay)); ...
    repmat(2, 1, length(PO2Delay))], [flowsDelay(1:6); PO2Delay], 'Color', ...
    'k', 'LineWidth', 0.5);
set(h(end), 'LineStyle', '--');
xticks([1 2]); yticks([0 2 4]);
ylim([0 4.5]); 
xticklabels({'Flow'; 'pO_{2} Mean'}); 
xtickangle(45);
ylabel('Delay from Ca^{2+} onset (s)');
set(gca, 'FontSize', 13);


subplot(224); hold on;
dip = [3 4 6]; ndip = [1 2 5 7 8 9];
b = bar(1, mean(flowsDelay(ndip)), 'FaceColor', 'flat');
b.CData = colors{1};
c = bar(2, mean(flowsDelay(dip)), 'FaceColor', 'flat');
c.CData = colors{2};
scatter(repmat(1, 1, length(flowsDelay(ndip))), flowsDelay(ndip), 40, 'k', ...
    'filled');
scatter(repmat(2, 1, length(flowsDelay(dip))), flowsDelay(dip), 40, 'k', ...
    'filled');
xticks([1 2]); yticks([0 1 2 3]);
ylim([0 3.5]); ylim([0 4]);
xticklabels({'Without dip'; 'With dip'}); 
xtickangle(45);
ylabel('Delay between Ca^{2+} - Flow (s)');
set(gca, 'FontSize', 13);

suptitle('FIGURE 2C');

disp('### Figure 2C ###');
[p, h] = signrank(flowsDelay(1:6), PO2Delay, 'tail', 'left');
disp('Wilcoxon signed rank test - One-Sided for Flow onset < PO2 onset');
disp(['Significant difference between Flow and PO2 onset, p = ' ...
    num2str(p)])
[h, p] = ttest2(flowsDelay(ndip),flowsDelay(dip), 'tail', 'left');
disp('Two-sampled t-test - One-Sided for Flow onset without dip < with dip');
disp(['Significant difference between without and with dip, p = ' ...
    num2str(p)])

end

