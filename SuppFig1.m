function SuppFig1()
% FIG3 Generating Supp Fig 1 from Aydin et al.
%
% function SuppFig1() = []
%
%   Author : Ali-Kemal Aydin, PhD student
%   Date : January 5th, 2021
%   Mail: ali-kemal.aydin@inserm.fr
%   Affiliation : U968, Institut de la Vision, Paris
%   License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
%       See LICENSE.txt or <a href="matlab:web('https://creativecommons.org/licenses/by/4.0/')">here</a> 
%       for a human-readable version.
%
%   DESCRIPTION : Generates the panels from Supp Fig 1 in Aydin et al.
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

%% Panel A
% This panel is not reproduced here as it is only an image.

%% Panel B

m = 'M1514_191120';
type = 'Interp';
xL = [0 29]; yL = [20 80];

% n = 3 acquisitions for this mouse
figure;
% Individual mouse
subplot(221); hold on;
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(d.(m).(type).Oxygen_ET_1p5V_2sec.PO2All(1, :), ...
    d.(m).(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, :), 'Color', ...
    [0.5 0.5 0.5]);
plot(d.(m).(type).Oxygen_ET_1p5V_2sec.PO2All(1, :), ...
    d.(m).(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, :), 'Color', 'k', ...
    'LineWidth', 2);
ylabel('\Delta Po_{2} Mean (mmHg)'); 
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

subplot(222); hold on;
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(d.(m).(type).Oxygen_ET_1p5V_2sec.PO2RBC(1, :), ...
    d.(m).(type).Oxygen_ET_1p5V_2sec.PO2RBC(2:end-2, :), 'Color', ...
    [0.5 0.5 0.5]);
plot(d.(m).(type).Oxygen_ET_1p5V_2sec.PO2RBC(1, :), ...
    d.(m).(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, :), 'Color', 'k', ...
    'LineWidth', 2);
ylabel(' \Delta Po_{2} RBC (mmHg)'); 
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

subplot(223); hold on;
patch('XData', [10 12 12 10], 'YData', [0 0 35 35], 'FaceColor', ...
    colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(d.(m).(type).Oxygen_ET_1p5V_2sec.Flow(1, :), ...
    d.(m).(type).Oxygen_ET_1p5V_2sec.Flow(2:end-2, :), 'Color', ...
    [0.5 0.5 0.5]);
plot(d.(m).(type).Oxygen_ET_1p5V_2sec.Flow(1, :), ...
    d.(m).(type).Oxygen_ET_1p5V_2sec.Flow(end-1, :), 'Color', 'k', ...
    'LineWidth', 2);
ylabel('\Delta Flow (RBC.s^{-1})'); xlabel('Time (s)'); 
xlim(xL); ylim([0 35]);
set(gca, 'FontSize', 13);

subplot(224); hold on;
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');
legend('AutoUpdate', 'off');
plot(d.(m).(type).Oxygen_ET_1p5V_2sec.PO2Inter(1, :), ...
    d.(m).(type).Oxygen_ET_1p5V_2sec.PO2Inter(2:end-3, :), 'Color', ...
    [0.5 0.5 0.5]);
legend('AutoUpdate', 'on');
plot(d.(m).(type).Oxygen_ET_1p5V_2sec.PO2Inter(1, :), ...
    d.(m).(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-2, :), 'Color', ...
    [0.5 0.5 0.5], 'DisplayName', 'Single acquisition');
plot(d.(m).(type).Oxygen_ET_1p5V_2sec.PO2Inter(1, :), ...
    d.(m).(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, :), 'Color', 'k', ...
    'LineWidth', 2, 'DisplayName', 'Average');
ylabel('\Delta Po_{2} Inter (mmHg)'); xlabel('Time (s)');
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

suptitle('SUPP FIG 1B');

%% 2C. Mean Responses in ZScore

type = 'ZScore';
tIndexPO2 = d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, :) >= 1 & d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, :) <= 29;
timePO2 = d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2);

exps.PO2RBC =  [d.M1393_021020.(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, tIndexPO2); ...
    d.M1511_211020.(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, tIndexPO2); ...
    d.M1512_281020.(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, tIndexPO2); ...
    d.M1514_191120.(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, tIndexPO2); ...
    d.M1804_080221.(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, tIndexPO2); ...
    d.M1854_090321.(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, tIndexPO2); ...
    d.M1873_220321.(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, tIndexPO2); ...
    d.M1952_240321.(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, tIndexPO2); ...
    mean([d.M1954_240321.(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, tIndexPO2); ...
    d.M1954_010421.(type).Oxygen_ET_1p5V_2sec.PO2RBC(end-1, tIndexPO2)], 1)];

exps.PO2Inter =  [d.M1393_021020.(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, tIndexPO2); ...
    d.M1511_211020.(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, tIndexPO2); ...
    d.M1512_281020.(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, tIndexPO2); ...
    d.M1514_191120.(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, tIndexPO2); ...
    d.M1804_080221.(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, tIndexPO2); ... % P_Inter values are out of line
    d.M1854_090321.(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, tIndexPO2); ...
    d.M1873_220321.(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, tIndexPO2); ...
    d.M1952_240321.(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, tIndexPO2); ...
    mean([d.M1954_240321.(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, tIndexPO2); ...
    d.M1954_010421.(type).Oxygen_ET_1p5V_2sec.PO2Inter(end-1, tIndexPO2)], 1)];

exps.Flow =  [d.M1393_021020.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); ...
    d.M1511_211020.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); ...
    d.M1512_281020.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); ...
    d.M1514_191120.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); ...
    d.M1804_080221.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); ... 
    d.M1854_090321.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); ...
    d.M1873_220321.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); ...
    d.M1952_240321.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); ...
    mean([d.M1954_240321.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2); ...
    d.M1954_010421.(type).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexPO2)], 1)];

xL = [0 29]; yL = [-3 12];

figure;
subplot(131); hold on;
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(timePO2, exps.Flow, 'Color', [0.5 0.5 0.5]);
plot(timePO2, mean(exps.Flow, 1), 'Color', 'k', 'LineWidth', 2);
xlim(xL); ylim(yL);
ylabel('\Delta Flow (SD)'); xlabel('Time (s)');
set(gca, 'FontSize', 13);

yL = [-3 7];
subplot(132); hold on;
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(timePO2, exps.PO2RBC, 'Color', [0.5 0.5 0.5]);
plot(timePO2, mean(exps.PO2RBC, 1), 'Color', 'k', 'LineWidth', 2);
ylabel('\Delta Po_{2} RBC (SD)'); xlabel('Time (s)');
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

subplot(133); hold on;
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'LineWidth', 0.5, ...
    'FaceAlpha', 0.5, 'DisplayName', 'Stimulus');
legend('AutoUpdate', 'off');
plot(timePO2, exps.PO2Inter(1:end-1, :), 'Color', [0.5 0.5 0.5]);
legend('AutoUpdate', 'on');
plot(timePO2, exps.PO2Inter(end, :), 'Color', [0.5 0.5 0.5], ...
    'DisplayName', 'Single mouse');
plot(timePO2, mean(exps.PO2Inter, 1), 'Color', 'k', 'LineWidth', 2, ...
    'DisplayName', 'Average');
ylabel('\Delta Po_{2} Inter (SD)'); xlabel('Time (s)');
xlim(xL); ylim(yL); legend();
set(gca, 'FontSize', 13);

suptitle('SUPP FIG 1C');


end