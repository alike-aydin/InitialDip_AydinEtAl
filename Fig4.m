function Fig4()
% FIG4 Generating Figure 4 from Aydin et al.
%
% function Fig4() = []
%
%   Author : Ali-Kemal Aydin, PhD student
%   Date : April 19th, 2021
%   Mail: ali-kemal.aydin@inserm.fr
%   Affiliation : U968, Institut de la Vision, Paris
%   License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
%       See LICENSE.txt or <a href="matlab:web('https://creativecommons.org/licenses/by/4.0/')">here</a>
%       for a human-readable version.
%
%   DESCRIPTION : Generates the panels from Figure 4 in Aydin et al.
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
% Data in this panel is from a single mouse, but merges two experiments
m1 = 'M1954_240321'; m2 = 'M1954_010421';
typeVasc = 'Interp'; typeCa = 'Delta';
xL = [0 29];

figure;
% Ca2+ responses
subplot(321); hold on;
yL = [-50 450];
% Oxygenated data
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
x = d.(m1).(typeCa).Oxygen_ET_1p5V_2sec.Ca(1, :);
y = mean([d.(m1).(typeCa).Oxygen_ET_1p5V_2sec.Ca(2:end-2, :); ...
    d.(m2).(typeCa).Oxygen_ET_1p5V_2sec.Ca(2:end-2, :)]);
ystd = std([d.(m1).(typeCa).Oxygen_ET_1p5V_2sec.Ca(2:end-2, :); ...
    d.(m2).(typeCa).Oxygen_ET_1p5V_2sec.Ca(2:end-2, :)], 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{1}, ... % Thanks to Paul Baudin, from the Charpier team at the Paris Brain Institute for this piece of code
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{1}, 'LineWidth', 2);

% Reduced O2 data
x = d.(m1).(typeCa).NoOxy_ET_1p5V_2sec.Ca(1, :);
y = mean([d.(m1).(typeCa).NoOxy_ET_1p5V_2sec.Ca(2:end-2, :); ...
    d.(m2).(typeCa).NoOxy_ET_1p5V_2sec.Ca(2:end-2, :)]);
ystd = std([d.(m1).(typeCa).NoOxy_ET_1p5V_2sec.Ca(2:end-2, :); ...
    d.(m2).(typeCa).NoOxy_ET_1p5V_2sec.Ca(2:end-2, :)], 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{2}, 'LineWidth', 2);
ylabel('\Delta Ca^{2+} (a.u.)');
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13);

% RBC velocity data
subplot(323); hold on;
yL = [0 1.5];
% Oxygenated data
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
x = d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, :);
y = mean([d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.RBC(2:end-2, :); ...
    d.(m2).(typeVasc).Oxygen_ET_1p5V_2sec.RBC(2:end-2, :)]);
ystd = std([d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.RBC(2:end-2, :); ...
    d.(m2).(typeVasc).Oxygen_ET_1p5V_2sec.RBC(2:end-2, :)], 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{1}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{1}, 'LineWidth', 2);

% Reduced O2 data
x = d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.RBC(1, :);
y = mean([d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.RBC(2:end-2, :); ...
    d.(m2).(typeVasc).NoOxy_ET_1p5V_2sec.RBC(2:end-2, :)]);
ystd = std([d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.RBC(2:end-2, :); ...
    d.(m2).(typeVasc).NoOxy_ET_1p5V_2sec.RBC(2:end-2, :)], 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{2}, 'LineWidth', 2);
ylabel('RBC Velocity (mm.s^{-1})');
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13);

% Flow data
subplot(325); hold on;
yL = [0 150];
% Oxygenated data
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
x = d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.Flow(1, :);
y = mean([d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.Flow(2:end-2, :); ...
    d.(m2).(typeVasc).Oxygen_ET_1p5V_2sec.Flow(2:end-2, :)]);
ystd = std([d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.Flow(2:end-2, :); ...
    d.(m2).(typeVasc).Oxygen_ET_1p5V_2sec.Flow(2:end-2, :)], 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{1}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{1}, 'LineWidth', 2);

% Reduced O2 data
x = d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.Flow(1, :);
y = mean([d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.Flow(2:end-2, :); ...
    d.(m2).(typeVasc).NoOxy_ET_1p5V_2sec.Flow(2:end-2, :)]);
ystd = std([d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.Flow(2:end-2, :); ...
    d.(m2).(typeVasc).NoOxy_ET_1p5V_2sec.Flow(2:end-2, :)], 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{2}, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('RBC Flow (RBC.s^{-1})');
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

% PO2 Mean data
subplot(322); hold on;
yL= [20 60];
% Oxygenated data
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
x = d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(1, :);
y = mean([d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, :); ...
    d.(m2).(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, :)]);
ystd = std([d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, :); ...
    d.(m2).(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, :)], 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{1}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{1}, 'LineWidth', 2);

% Reduced O2 data
x = d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(1, :);
y = mean([d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, :); ...
    d.(m2).(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, :)]);
ystd = std([d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, :); ...
    d.(m2).(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, :)], 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{2}, 'LineWidth', 2);
ylabel('Po_{2} Mean (mmHg)');
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13);

% PO2 RBC Data
subplot(324); hold on;
yL= [30 110];
% Oxygenated data
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
x = d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.PO2RBC(1, :);
y = mean([d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.PO2RBC(2:end-2, :); ...
    d.(m2).(typeVasc).Oxygen_ET_1p5V_2sec.PO2RBC(2:end-2, :)]);
ystd = std([d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.PO2RBC(2:end-2, :); ...
    d.(m2).(typeVasc).Oxygen_ET_1p5V_2sec.PO2RBC(2:end-2, :)], 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{1}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{1}, 'LineWidth', 2);

% Reduced O2 data
x = d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.PO2RBC(1, :);
y = mean([d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.PO2RBC(2:end-2, :); ...
    d.(m2).(typeVasc).NoOxy_ET_1p5V_2sec.PO2RBC(2:end-2, :)]);
ystd = std([d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.PO2RBC(2:end-2, :); ...
    d.(m2).(typeVasc).NoOxy_ET_1p5V_2sec.PO2RBC(2:end-2, :)], 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{2}, 'LineWidth', 2);
ylabel('Po_{2} RBC (mmHg)');
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13);

% PO2 Inter Data
subplot(326); hold on;
yL= [0 50];
% Oxygenated data
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');
x = d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.PO2Inter(1, :);
y = mean([d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.PO2Inter(2:end-2, :); ...
    d.(m2).(typeVasc).Oxygen_ET_1p5V_2sec.PO2Inter(2:end-2, :)]);
ystd = std([d.(m1).(typeVasc).Oxygen_ET_1p5V_2sec.PO2Inter(2:end-2, :); ...
    d.(m2).(typeVasc).Oxygen_ET_1p5V_2sec.PO2Inter(2:end-2, :)], 0, 1);
legend('AutoUpdate', 'off');
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{1}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
legend('AutoUpdate', 'on');
plot(x, y, 'Color', colors{1}, 'LineWidth', 2, 'DisplayName', 'Oxygenated');

x = d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.PO2Inter(1, :);
y = mean([d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.PO2Inter(2:end-2, :); ...
    d.(m2).(typeVasc).NoOxy_ET_1p5V_2sec.PO2Inter(2:end-2, :)]);
ystd = std([d.(m1).(typeVasc).NoOxy_ET_1p5V_2sec.PO2Inter(2:end-2, :); ...
    d.(m2).(typeVasc).NoOxy_ET_1p5V_2sec.PO2Inter(2:end-2, :)], 0, 1);
legend('AutoUpdate', 'off');
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
legend('AutoUpdate', 'on');
plot(x, y, 'Color', colors{2}, 'LineWidth', 2, 'DisplayName', ...
    'Reduced O_{2}');
legend(); xlabel('Time (s)'); ylabel('PO_{2} Inter (mmHg)');
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13);

suptitle('Figure 4A');

%% Panel B & C

mAll = {'M1804_080221'; 'M1854_090321'; 'M1873_220321'};
mDouble = {{'M1514_051120'; 'M1514_191120'}; {'M1952_240321'; 'M1952_310321'}; {'M1954_240321'; 'M1954_010421'}};
mNoFlow = {'M1397_041120'; 'M1594_261120'; 'M1951_180321'};
mKineticBsl = {'M1594_161220'; 'M1594_211220'; 'M1514_181220'};

means = struct();
means.RBC.Oxy = []; means.RBC.NoOxy = [];
means.Flow.Oxy = []; means.Flow.NoOxy = [];
means.PO2All.Oxy = []; means.PO2All.NoOxy = [];
means.PO2RBC.Oxy = []; means.PO2RBC.NoOxy = [];
means.PO2Inter.Oxy = []; means.PO2Inter.NoOxy = [];

tBSL = data.(mAll{1}).Interp.Oxygen_ET_1p5V_2sec.RBC(1, :) <= 10;
conds = {'Oxy'; 'NoOxy'};

ca = struct(); rbc = struct(); flow = struct(); PO2All = struct(); PO2RBC = struct(); PO2Inter = struct();


for i=1:length(mAll)
    d = data.(mAll{i});

    rbc.oxy = mean(d.Interp.Oxygen_ET_1p5V_2sec.RBC(2:end-2, tBSL), 2);
    rbc.noOxy = mean(d.Interp.NoOxy_ET_1p5V_2sec.RBC(2:end-2, tBSL), 2);
    means.RBC.Oxy = [means.RBC.Oxy mean(rbc.oxy)];
    means.RBC.NoOxy = [means.RBC.NoOxy mean(rbc.noOxy)];

    flow.oxy = mean(d.Interp.Oxygen_ET_1p5V_2sec.Flow(2:end-2, tBSL), 2);
    flow.noOxy = mean(d.Interp.NoOxy_ET_1p5V_2sec.Flow(2:end-2, tBSL), 2);
    means.Flow.Oxy = [means.Flow.Oxy mean(flow.oxy)];
    means.Flow.NoOxy = [means.Flow.NoOxy mean(flow.noOxy)];

    PO2All.oxy = mean(d.Interp.Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tBSL), 2);
    PO2All.noOxy = mean(d.Interp.NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tBSL), 2);
    means.PO2All.Oxy = [means.PO2All.Oxy mean(PO2All.oxy)];
    means.PO2All.NoOxy = [means.PO2All.NoOxy mean(PO2All.noOxy)];

    PO2RBC.oxy = mean(d.Interp.Oxygen_ET_1p5V_2sec.PO2RBC(2:end-2, tBSL), 2);
    PO2RBC.noOxy = mean(d.Interp.NoOxy_ET_1p5V_2sec.PO2RBC(2:end-2, tBSL), 2);
    means.PO2RBC.Oxy = [means.PO2RBC.Oxy mean(PO2RBC.oxy)];
    means.PO2RBC.NoOxy = [means.PO2RBC.NoOxy mean(PO2RBC.noOxy)];

    PO2Inter.oxy = mean(d.Interp.Oxygen_ET_1p5V_2sec.PO2Inter(2:end-2, tBSL), 2);
    PO2Inter.noOxy = mean(d.Interp.NoOxy_ET_1p5V_2sec.PO2Inter(2:end-2, tBSL), 2);
    means.PO2Inter.Oxy = [means.PO2Inter.Oxy mean(PO2Inter.oxy)];
    means.PO2Inter.NoOxy = [means.PO2Inter.NoOxy mean(PO2Inter.noOxy)];

end

% Averaging the multiple experiments on the same mice
for i=1:length(mDouble)
    d1 = data.(mDouble{i}{1});
    d2 = data.(mDouble{i}{2});

    rbc.oxy = [mean(mean(d1.Interp.Oxygen_ET_1p5V_2sec.RBC(2:end-2, tBSL), 2)) mean(mean(d2.Interp.Oxygen_ET_1p5V_2sec.RBC(2:end-2, tBSL), 2))];
    rbc.noOxy = [mean(mean(d1.Interp.NoOxy_ET_1p5V_2sec.RBC(2:end-2, tBSL), 2)) mean(mean(d2.Interp.NoOxy_ET_1p5V_2sec.RBC(2:end-2, tBSL), 2))];
    means.RBC.Oxy = [means.RBC.Oxy mean(rbc.oxy)];
    means.RBC.NoOxy = [means.RBC.NoOxy mean(rbc.noOxy)];

    if i ~= 2 % No Flow in the 2nd exp of 1952
        flow.oxy = [mean(mean(d1.Interp.Oxygen_ET_1p5V_2sec.Flow(2:end-2, tBSL), 2)) mean(mean(d2.Interp.Oxygen_ET_1p5V_2sec.Flow(2:end-2, tBSL), 2))];
        flow.noOxy = [mean(mean(d1.Interp.NoOxy_ET_1p5V_2sec.Flow(2:end-2, tBSL), 2)) mean(mean(d2.Interp.NoOxy_ET_1p5V_2sec.Flow(2:end-2, tBSL), 2))];
    else
        flow.oxy = mean(d1.Interp.Oxygen_ET_1p5V_2sec.Flow(2:end-2, tBSL), 2);
        flow.noOxy = mean(d1.Interp.NoOxy_ET_1p5V_2sec.Flow(2:end-2, tBSL), 2);
    end
    means.Flow.Oxy = [means.Flow.Oxy mean(flow.oxy)];
    means.Flow.NoOxy = [means.Flow.NoOxy mean(flow.noOxy)];

    PO2All.oxy = [mean(mean(d1.Interp.Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tBSL), 2)) mean(mean(d2.Interp.Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tBSL), 2))];
    PO2All.noOxy = [mean(mean(d1.Interp.NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tBSL), 2)) mean(mean(d2.Interp.NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tBSL), 2))];
    means.PO2All.Oxy = [means.PO2All.Oxy mean(PO2All.oxy)];
    means.PO2All.NoOxy = [means.PO2All.NoOxy mean(PO2All.noOxy)];

    if i ~= 2 % No Flow in the 2nd exp of 1952
        PO2RBC.oxy = [mean(mean(d1.Interp.Oxygen_ET_1p5V_2sec.PO2RBC(2:end-2, tBSL), 2)) mean(mean(d2.Interp.Oxygen_ET_1p5V_2sec.PO2RBC(2:end-2, tBSL), 2))];
        PO2RBC.noOxy = [mean(mean(d1.Interp.NoOxy_ET_1p5V_2sec.PO2RBC(2:end-2, tBSL), 2)) mean(mean(d2.Interp.NoOxy_ET_1p5V_2sec.PO2RBC(2:end-2, tBSL), 2))];
    else
        PO2RBC.oxy = mean(d1.Interp.Oxygen_ET_1p5V_2sec.PO2RBC(2:end-2, tBSL), 2);
        PO2RBC.noOxy = mean(d1.Interp.NoOxy_ET_1p5V_2sec.PO2RBC(2:end-2, tBSL), 2);
    end
    means.PO2RBC.Oxy = [means.PO2RBC.Oxy mean(PO2RBC.oxy)];
    means.PO2RBC.NoOxy = [means.PO2RBC.NoOxy mean(PO2RBC.noOxy)];

    if i ~= 2
        PO2Inter.oxy = [mean(mean(d1.Interp.Oxygen_ET_1p5V_2sec.PO2Inter(2:end-2, tBSL), 2)) mean(mean(d2.Interp.Oxygen_ET_1p5V_2sec.PO2Inter(2:end-2, tBSL), 2))];
        PO2Inter.noOxy = [mean(mean(d1.Interp.NoOxy_ET_1p5V_2sec.PO2Inter(2:end-2, tBSL), 2)) mean(mean(d2.Interp.NoOxy_ET_1p5V_2sec.PO2Inter(2:end-2, tBSL), 2))];
    else
        PO2Inter.oxy = mean(d1.Interp.Oxygen_ET_1p5V_2sec.PO2Inter(2:end-2, tBSL), 2);
        PO2Inter.noOxy = mean(d1.Interp.NoOxy_ET_1p5V_2sec.PO2Inter(2:end-2, tBSL), 2);
    end
    means.PO2Inter.Oxy = [means.PO2Inter.Oxy mean(PO2Inter.oxy)];
    means.PO2Inter.NoOxy = [means.PO2Inter.NoOxy mean(PO2Inter.noOxy)];
end

for i=1:length(mNoFlow)
    d = data.(mNoFlow{i});

    rbc.oxy = mean(d.Interp.Oxygen_ET_1p5V_2sec.RBC(2:end-2, tBSL), 2);
    rbc.noOxy = mean(d.Interp.NoOxy_ET_1p5V_2sec.RBC(2:end-2, tBSL), 2);
    means.RBC.Oxy = [means.RBC.Oxy mean(rbc.oxy)];
    means.RBC.NoOxy = [means.RBC.NoOxy mean(rbc.noOxy)];

    PO2All.oxy = mean(d.Interp.Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tBSL), 2);
    PO2All.noOxy = mean(d.Interp.NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tBSL), 2);
    means.PO2All.Oxy = [means.PO2All.Oxy mean(PO2All.oxy)];
    means.PO2All.NoOxy = [means.PO2All.NoOxy mean(PO2All.noOxy)];
end


for i=1:length(mKineticBsl)
    d = data.(mKineticBsl{i});

    PO2RBC.Oxy = d.Raw.Oxygen_CapStart.EAT.RBC;
    PO2Inter.Oxy = d.Raw.Oxygen_CapStart.EAT.Inter;
    PO2All.Oxy = d.Raw.Oxygen_CapStart.EAT.Mean;

    rbc.Oxy = mean(d.Interp.Oxygen.KPtSpeed(2:end-2, :), 2);
    density.Oxy = mean(d.Interp.Oxygen.KPtDensity(2:end-2, :), 2);
    flow.Oxy = mean(d.Interp.Oxygen.KPtDensity(2:end-2, :) .* d.Interp.Oxygen.KPtSpeed(2:end-2, :), 2);

    PO2RBC.NoOxy = d.Raw.NoOxy_CapStart.EAT.RBC;
    PO2Inter.NoOxy = d.Raw.NoOxy_CapStart.EAT.Inter;
    PO2All.NoOxy = d.Raw.NoOxy_CapStart.EAT.Mean;

    rbc.NoOxy = mean(d.Interp.NoOxy.KPtSpeed(2:end-2, :), 2);
    density.NoOxy = mean(d.Interp.NoOxy.KPtDensity(2:end-2, :), 2);
    flow.NoOxy = mean(d.Interp.NoOxy.KPtDensity(2:end-2, :) .* d.Interp.NoOxy.KPtSpeed(2:end-2, :), 2);


    means.RBC.Oxy = [means.RBC.Oxy mean(rbc.Oxy)];
    means.Flow.Oxy = [means.Flow.Oxy mean(flow.Oxy)];
    means.PO2All.Oxy = [means.PO2All.Oxy mean(PO2All.Oxy)];
    means.PO2RBC.Oxy = [means.PO2RBC.Oxy mean(PO2RBC.Oxy)];
    means.PO2Inter.Oxy = [means.PO2Inter.Oxy mean(PO2Inter.Oxy)];

    means.RBC.NoOxy = [means.RBC.NoOxy mean(rbc.NoOxy)];
    means.Flow.NoOxy = [means.Flow.NoOxy mean(flow.NoOxy)];
    means.PO2All.NoOxy = [means.PO2All.NoOxy mean(PO2All.NoOxy)];
    means.PO2RBC.NoOxy = [means.PO2RBC.NoOxy mean(PO2RBC.NoOxy)];
    means.PO2Inter.NoOxy = [means.PO2Inter.NoOxy mean(PO2Inter.NoOxy)];
end

yL = [0 85];
figure;
subplot(131); hold on;
b = bar([mean(means.PO2All.Oxy) mean(means.PO2All.NoOxy)], ...
    'FaceColor', 'flat');
b.CData(2, :) = colors{2};
scatter(repelem(1, 1, length(means.PO2All.Oxy)), ...
    means.PO2All.Oxy, 20, 'k', 'filled');
scatter(repelem(2, 1, length(means.PO2All.NoOxy)), ...
    means.PO2All.NoOxy, 20, 'k', 'filled');
h = line([repmat(1, 1, length(means.PO2All.Oxy)); ...
    repmat(2, 1, length(means.PO2All.NoOxy))], ...
    [means.PO2All.Oxy; means.PO2All.NoOxy], 'Color', 'k', 'LineWidth', 0.5);
set(h(6), 'LineStyle', '--'); % This is the mouse from Panel A
xticks([1:2]); xlim([0.5 2.5]); xticklabels({}); ylim(yL);
ylabel('pO_{2} Mean (mmHg)');
set(gca, 'FontSize', 13);

subplot(132); hold on;
b = bar([mean(means.PO2RBC.Oxy) mean(means.PO2RBC.NoOxy)], ...
    'FaceColor', 'flat');
b.CData(2, :) = colors{2};
scatter(repelem(1, 1, length(means.PO2RBC.Oxy)), ...
    means.PO2RBC.Oxy, 20, 'k', 'filled');
scatter(repelem(2, 1, length(means.PO2RBC.NoOxy)), ...
    means.PO2RBC.NoOxy, 20, 'k', 'filled');
h = line([repmat(1, 1, length(means.PO2RBC.Oxy)); ...
    repmat(2, 1, length(means.PO2RBC.NoOxy))], ...
    [means.PO2RBC.Oxy; means.PO2RBC.NoOxy], 'Color', 'k', 'LineWidth', 0.5);
set(h(6), 'LineStyle', '--'); % This is the mouse from Panel A
xticks([1:2]); xlim([0.5 2.5]); ylim(yL);
xticklabels({'Oxygenated'; 'Deoxygenated'}); xtickangle(45);
ylabel('pO_{2} RBC (mmHg)');
set(gca, 'FontSize', 13);

subplot(133); hold on;
b = bar([mean(means.PO2Inter.Oxy) mean(means.PO2Inter.NoOxy)], ...
    'FaceColor', 'flat');
b.CData(2, :) = colors{2};
scatter(repelem(1, 1, length(means.PO2Inter.Oxy)), ...
    means.PO2Inter.Oxy, 20, 'k', 'filled');
scatter(repelem(2, 1, length(means.PO2Inter.NoOxy)), ...
    means.PO2Inter.NoOxy, 20, 'k', 'filled');
h = line([repmat(1, 1, length(means.PO2Inter.Oxy)); ...
    repmat(2, 1, length(means.PO2Inter.NoOxy))], ...
    [means.PO2Inter.Oxy; means.PO2Inter.NoOxy], 'Color', 'k', 'LineWidth', 0.5);
set(h(6), 'LineStyle', '--'); % This is the mouse from Panel A
xticks([1:2]); xlim([0.5 2.5]); xticklabels({});
 ylim(yL);
ylabel('pO_{2} Inter RBC (mmHg)');
set(gca, 'FontSize', 13);

suptitle('FIGURE 4B')

p = struct();
[p.PO2All, h] = signrank(means.PO2All.Oxy, means.PO2All.NoOxy, 'tail', 'right');
[p.PO2RBC, h] = signrank(means.PO2RBC.Oxy, means.PO2RBC.NoOxy, 'tail', 'right');
[p.PO2Inter, h] = signrank(means.PO2Inter.Oxy, means.PO2Inter.NoOxy, 'tail', 'right');

disp('### Figure 4B ###')
disp('Wilcoxon signed rank test - One-Sided for Baseline values of PO2 Oxy > NoOxy');
disp(['Significant differences between baseline values of PO2 Mean, RBC' ...
    ' and Inter when mouse is oxygenated or with reduced O2, p = ' ...
    num2str(p.PO2All) ', ' num2str(p.PO2RBC) ' and ' num2str(p.PO2Inter) ...
    ' respectively.']);

figure;

subplot(211); hold on;
b = bar([mean(means.RBC.Oxy) mean(means.RBC.NoOxy)], 'FaceColor', 'flat');
b.CData(2, :) = colors{2};
scatter(repelem(1, 1, length(means.RBC.Oxy)), means.RBC.Oxy, 20, ...
    'k', 'filled');
scatter(repelem(2, 1, length(means.RBC.NoOxy)), means.RBC.NoOxy, 20, ...
    'k', 'filled');
h = line([repmat(1, 1, length(means.RBC.Oxy)); ...
    repmat(2, 1, length(means.RBC.NoOxy))], ...
    [means.RBC.Oxy; means.RBC.NoOxy], 'Color', 'k', 'LineWidth', 0.5);
set(h(6), 'LineStyle', '--'); % This is the mouse from Panel A
xticks([1:2]); xlim([0.5 2.5]); xticklabels({});
ylabel('RBC velocity (mm.s^{-1})');
set(gca, 'FontSize', 13);

subplot(212); hold on;
b = bar([mean(means.Flow.Oxy) mean(means.Flow.NoOxy)], ...
    'FaceColor', 'flat');
b.CData(2, :) = colors{2};
scatter(repelem(1, 1, length(means.Flow.Oxy)), means.Flow.Oxy, 20, ...
    'k', 'filled');
scatter(repelem(2, 1, length(means.Flow.NoOxy)), means.Flow.NoOxy, 20, ...
    'k', 'filled');
h = line([repmat(1, 1, length(means.Flow.Oxy)); ...
    repmat(2, 1, length(means.Flow.NoOxy))], ...
    [means.Flow.Oxy; means.Flow.NoOxy], 'Color', 'k', 'LineWidth', 0.5);
set(h(6), 'LineStyle', '--'); % This is the mouse from Panel A
xticks([1:2]); xlim([0.5 2.5]);
ylabel('Flow (RBC.s^{-1})');
set(gca, 'FontSize', 13);

suptitle('FIGURE 4C')

[p.RBC, h] = signrank(means.RBC.Oxy, means.RBC.NoOxy, 'tail', 'left');
[p.Flow, h] = signrank(means.Flow.Oxy, means.Flow.NoOxy, 'tail', 'left');


disp('### Figure 4C ###')
disp('Wilcoxon signed rank test - One-Sided for Baseline values of PO2 Oxy < NoOxy');
disp(['Significant differences between baseline values of RBC velocity' ...
    ' and Flow when mouse is oxygenated or with reduced O2, p = ' ...
    num2str(p.RBC) ' and ' num2str(p.Flow) ' respectively.']);
%% Panel D

caTTP = struct(); caTTP.Oxy = []; caTTP.NoOxy = [];
maxCa = struct(); maxCa.Oxy = []; maxCa.NoOxy = [];

typeCa = 'Delta';
xL = [0 29];
d = data;

%%%% Calcium data - Missing Exp from 1875 - 9 mice & 12 exp %%%%
tIndexCa = d.M1397_041120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(1, :) >= xL(1) ...
    & d.M1397_041120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(1, :) <= xL(2);
timeCa = d.M1397_041120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(1, tIndexCa);

Paired.Ca.Oxy = [d.M1397_041120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    mean([d.M1514_051120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1514_191120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1); ...
    d.M1594_261120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1804_080221.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1854_090321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1951_180321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1873_220321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    mean([d.M1952_240321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1952_310321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1); ...
    mean([d.M1954_240321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1954_010421.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1);];

Paired.Ca.NoOxy = [d.M1397_041120.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    mean([d.M1514_051120.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1514_191120.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1); ...
    d.M1594_261120.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1804_080221.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1854_090321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1951_180321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1873_220321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    mean([d.M1952_240321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1952_310321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1); ...
    mean([d.M1954_240321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1954_010421.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1);];

for i=1:size(Paired.Ca.NoOxy, 1)
    [M, I] = max(Paired.Ca.Oxy(i, timeCa < 12));
    caTTP.Oxy(i) = timeCa(I)-10;
    maxCa.Oxy(i) = M;

    [M, I] = max(Paired.Ca.NoOxy(i, timeCa < 12));
    caTTP.NoOxy(i) = timeCa(I)-10;
    maxCa.NoOxy(i) = M;
end

figure;
subplot(211); hold on;
b = bar([mean(caTTP.Oxy) mean(caTTP.NoOxy)], 'FaceColor', 'flat');
b.CData(2, :) = colors{2};
scatter(repmat(1, 1, length(caTTP.Oxy)), caTTP.Oxy, 20, ...
    'k', 'filled');
scatter(repmat(2, 1, length(caTTP.NoOxy)), caTTP.NoOxy, 20, ...
    'k', 'filled');
h = line([repmat(1, 1, length(caTTP.Oxy)); ...
    repmat(2, 1, length(caTTP.Oxy))], ...
    [caTTP.Oxy; caTTP.NoOxy], 'Color', 'k');
set(h(end), 'LineStyle', '--'); % This is the mouse from Panel A
xticklabels({}); xticks([1:2]); xtickangle(45);
xlim([0.5 2.5]); ylim([0 2]);
ylabel('Ca^{2+} time-to-peak (s)');
set(gca, 'FontSize', 13);

subplot(212); hold on;
b = bar([mean(maxCa.Oxy) mean(maxCa.NoOxy)], 'FaceColor', 'flat');
b.CData(2, :) = colors{2};
scatter(repmat(1, 1, length(maxCa.Oxy)), maxCa.Oxy, 20, ...
    'k', 'filled');
scatter(repmat(2, 1, length(maxCa.NoOxy)), maxCa.NoOxy, 20, ...
    'k', 'filled');
h = line([repmat(1, 1, length(maxCa.Oxy)); ...
    repmat(2, 1, length(maxCa.Oxy))], ...
    [maxCa.Oxy; maxCa.NoOxy], 'Color', 'k');
set(h(end), 'LineStyle', '--'); % This is the mouse from Panel A
ylabel('Ca^{2+} Amplitude (a.u.)');
xticklabels({}); xticks([1:2]); xlim([0.5 2.5]); ylim([0 650]);
set(gca, 'FontSize', 13);

suptitle('FIGURE 4D');

p = struct();
[p.caTTP, h] = signrank(caTTP.Oxy, caTTP.NoOxy);
[p.maxCa, h] = signrank(maxCa.Oxy, maxCa.NoOxy);

disp('### Figure 4D ###')
disp('Wilcoxon signed rank test - Two-Sided for onset of responses Oxy =/= NoOxy');
disp(['No differences between Ca2+ time to peak' ...
    ' and amplitude when mouse is oxygenated or with reduced O2, p = ' ...
    num2str(p.caTTP) ' and ' num2str(p.maxCa) ' respectively.']);

%% Panel E
flowsDelay = struct();
flowsDelay.Oxy = [];
flowsDelay.NoOxy = [];
rbcDelay.Oxy = [];
rbcDelay.NoOxy = [];
d = data;
xL = [0 29];

typeVasc = 'Delta';
tIndexVasc = d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, :) >= xL(1) ...
    & d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, :) <= xL(2);
timeVasc = d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc);

Paired.RBC.Oxy = [d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    mean([d.M1514_051120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1514_191120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1); ...
    d.M1594_261120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1804_080221.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1854_090321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1951_180321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1873_220321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    mean([d.M1952_240321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1952_310321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1); ...
    mean([d.M1954_240321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1954_010421.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1);];

flowsDelay.Oxy(1) = getThreshold_fitSig(squeeze(d.M1512_281020.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(d.M1512_281020.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc))', [5 10], [5 18], 0.2, false);
flowsDelay.NoOxy(1) = getThreshold_fitSig(squeeze(d.M1512_281020.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(d.M1512_281020.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc))', [5 10], [5 18], 0.2, false);
flowsDelay.Oxy(2) = getThreshold_fitSig(squeeze(d.M1514_051120.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(mean([d.M1514_051120.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc); d.M1514_191120.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc)], 1))', [5 10], [8 18], 0.2, false);
flowsDelay.NoOxy(2) = getThreshold_fitSig(squeeze(d.M1514_051120.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(mean([d.M1514_051120.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc); d.M1514_191120.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc)], 1))', [5 10], [8 18], 0.2, false);
flowsDelay.Oxy(3) = getThreshold_fitSig(squeeze(d.M1804_080221.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(d.M1804_080221.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc))', [5 10], [8 18], 0.2, false);
flowsDelay.NoOxy(3) = getThreshold_fitSig(squeeze(d.M1804_080221.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(d.M1804_080221.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc))', [5 10], [8 18], 0.2, false);
flowsDelay.Oxy(4) = getThreshold_fitSig(squeeze(d.M1854_090321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(d.M1854_090321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc))', [5 10], [8 18], 0.2, false);
flowsDelay.NoOxy(4) = getThreshold_fitSig(squeeze(d.M1854_090321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(d.M1854_090321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc))', [5 10], [8 18], 0.2, false);
flowsDelay.Oxy(5) = getThreshold_fitSig(squeeze(d.M1873_220321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(d.M1873_220321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc))', [5 10], [8 18], 0.2, false);
flowsDelay.NoOxy(5) = getThreshold_fitSig(squeeze(d.M1873_220321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(d.M1873_220321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc))', [5 10], [8 18], 0.2, false);
flowsDelay.Oxy(6) = getThreshold_fitSig(squeeze(d.M1952_240321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(d.M1952_240321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc))', [5 10], [8 18], 0.2, false);
flowsDelay.NoOxy(6) = getThreshold_fitSig(squeeze(d.M1952_240321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(d.M1952_240321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc))', [5 10], [1 18], 0.2, false);
flowsDelay.Oxy(7) = getThreshold_fitSig(squeeze(d.M1954_240321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(mean([d.M1954_240321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc); d.M1954_010421.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc)], 1))', [5 10], [8 18], 0.2, false);
flowsDelay.NoOxy(7) = getThreshold_fitSig(squeeze(d.M1954_240321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(1, tIndexVasc)), squeeze(mean([d.M1954_240321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc); d.M1954_010421.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc)], 1))', [5 10], [8 18], 0.2, false);

rbcDelay.Oxy(1) = getThreshold_fitSig(squeeze(d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, :)), squeeze(d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, :))', [5 10], [1 20], 0.2, false);
rbcDelay.NoOxy(1) = getThreshold_fitSig(squeeze(d.M1397_041120.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(1, :)), squeeze(d.M1397_041120.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, :))', [5 10], [1 30], 0.2, false);
rbcDelay.Oxy(2) = getThreshold_fitSig(squeeze(d.M1514_051120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(mean([d.M1514_051120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); d.M1514_191120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1))', [5 10], [8 18], 0.2, false);
rbcDelay.NoOxy(2) = getThreshold_fitSig(squeeze(d.M1514_051120.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(mean([d.M1514_051120.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); d.M1514_191120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1))', [5 10], [8 18], 0.2, false);
rbcDelay.Oxy(3) = getThreshold_fitSig(squeeze(d.M1594_261120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, :)), squeeze(d.M1594_261120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, :))', [5 10], [1 15], 0.2, false);
rbcDelay.NoOxy(3) = getThreshold_fitSig(squeeze(d.M1594_261120.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(1, :)), squeeze(d.M1594_261120.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, :))', [5 10], [1 17], 0.2, false);
rbcDelay.Oxy(4) = getThreshold_fitSig(squeeze(d.M1804_080221.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(d.M1804_080221.Delta.Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc))', [5 10], [8 15], 0.2, false);
rbcDelay.NoOxy(4) = getThreshold_fitSig(squeeze(d.M1804_080221.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(d.M1804_080221.Delta.NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc))', [5 10], [8 15], 0.2, false);
rbcDelay.Oxy(5) = getThreshold_fitSig(squeeze(d.M1854_090321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(d.M1854_090321.Delta.Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc))', [5 10], [8 15], 0.2, false);
rbcDelay.NoOxy(5) = getThreshold_fitSig(squeeze(d.M1854_090321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(d.M1854_090321.Delta.NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc))', [5 10], [8 18], 0.2, false);
rbcDelay.Oxy(6) = getThreshold_fitSig(squeeze(d.M1873_220321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(d.M1873_220321.Delta.Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc))', [5 10], [8 15], 0.2, false);
rbcDelay.NoOxy(6) = getThreshold_fitSig(squeeze(d.M1873_220321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(d.M1873_220321.Delta.NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc))', [5 10], [8 15], 0.2, false);
rbcDelay.Oxy(7) = getThreshold_fitSig(squeeze(d.M1951_180321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(d.M1951_180321.Delta.Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc))', [5 10], [8 17], 0.2, false);
rbcDelay.NoOxy(7) = getThreshold_fitSig(squeeze(d.M1951_180321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(d.M1951_180321.Delta.NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc))', [5 10], [8 17], 0.2, false);
rbcDelay.Oxy(8) = getThreshold_fitSig(squeeze(d.M1952_240321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(mean([d.M1952_240321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); d.M1952_310321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1))', [5 10], [8 15], 0.2, false);
rbcDelay.NoOxy(8) = getThreshold_fitSig(squeeze(d.M1952_240321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(mean([d.M1952_240321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); d.M1952_310321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1))', [5 10], [8 18], 0.2, false);
rbcDelay.Oxy(9) = getThreshold_fitSig(squeeze(d.M1954_240321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(mean([d.M1954_240321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); d.M1954_010421.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1))', [5 10], [8 15], 0.2, false);
rbcDelay.NoOxy(9) = getThreshold_fitSig(squeeze(d.M1954_240321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(1, tIndexVasc)), squeeze(mean([d.M1954_240321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); d.M1954_010421.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1))', [5 10], [1 20], 0.2, false);

% Normalizing for the onset of Ca2+ at 10 s
flowsDelay.Oxy = flowsDelay.Oxy - 10;
flowsDelay.NoOxy = flowsDelay.NoOxy-10;
rbcDelay.Oxy = rbcDelay.Oxy-10;
rbcDelay.NoOxy = rbcDelay.NoOxy-10;

figure;
subplot(211); hold on;
b = bar([mean(rbcDelay.Oxy) mean(rbcDelay.NoOxy)], 'FaceColor', 'flat');
b.CData(2, :) = colors{2};
scatter(repmat(1, 1, length(rbcDelay.Oxy)), rbcDelay.Oxy, 20, ...
    'k', 'filled');
scatter(repmat(2, 1, length(rbcDelay.NoOxy)), rbcDelay.NoOxy, 20, ...
    'k', 'filled');
h = line([repmat(1, 1, length(rbcDelay.Oxy)); ...
    repmat(2, 1, length(rbcDelay.Oxy))], ...
    [rbcDelay.Oxy; rbcDelay.NoOxy], 'Color', 'k');
set(h(end), 'LineStyle', '--'); % This is the mouse from Panel A
xticklabels({'Oxygenated'; 'Deoxygenated'}); xticks([1:2]); xtickangle(45);
xlim([0.5 2.5]); ylim([0 3]);
ylabel('Delay Velocity onset (s)');
set(gca, 'FontSize', 13);

subplot(212); hold on;
b = bar([mean(flowsDelay.Oxy) mean(flowsDelay.NoOxy)], 'FaceColor', 'flat');
b.CData(2, :) = colors{2};
scatter(repmat(1, 1, length(flowsDelay.Oxy)), flowsDelay.Oxy, 20, ...
    'k', 'filled');
scatter(repmat(2, 1, length(flowsDelay.NoOxy)), flowsDelay.NoOxy, 20, ...
    'k', 'filled');
h = line([repmat(1, 1, length(flowsDelay.Oxy)); ...
    repmat(2, 1, length(flowsDelay.Oxy))], ...
    [flowsDelay.Oxy; flowsDelay.NoOxy], 'Color', 'k');
set(h(end), 'LineStyle', '--'); % This is the mouse from Panel A
ylabel('Delay Flow onset (s)');
xticklabels({}); xticks([1:2]); xlim([0.5 2.5]); ylim([0 3]);
set(gca, 'FontSize', 13);

suptitle('FIGURE 4E');

p = struct();
[p.vel, h] = signrank(rbcDelay.Oxy, rbcDelay.NoOxy);
[p.flow, h] = signrank(flowsDelay.Oxy, flowsDelay.NoOxy);

disp('### Figure 4D ###')
disp('Wilcoxon signed rank test - Two-Sided for onset of responses Oxy =/= NoOxy');
disp(['No differences between velocity and flow onset' ...
    ' when mouse is oxygenated or with reduced O2, p = ' ...
    num2str(p.vel) ' and ' num2str(p.flow) ' respectively.']);

%% Panel F
typeCa = 'Delta';
typeVasc = 'Delta';
xL = [0 29];
d = data;

%%%% Calcium data - Missing Exp from 1875 - 9 mice & 12 exp %%%%
tIndexCa = d.M1397_041120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(1, :) >= xL(1) ...
    & d.M1397_041120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(1, :) <= xL(2);
timeCa = d.M1397_041120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(1, tIndexCa);

Paired.Ca.Oxy = [d.M1397_041120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    mean([d.M1514_051120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1514_191120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1); ...
    d.M1594_261120.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1804_080221.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1854_090321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1951_180321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1873_220321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    mean([d.M1952_240321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1952_310321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1); ...
    mean([d.M1954_240321.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1954_010421.(typeCa).Oxygen_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1);];

Paired.Ca.NoOxy = [d.M1397_041120.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    mean([d.M1514_051120.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1514_191120.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1); ...
    d.M1594_261120.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1804_080221.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1854_090321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1951_180321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1873_220321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    mean([d.M1952_240321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1952_310321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1); ...
    mean([d.M1954_240321.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa); ...
    d.M1954_010421.(typeCa).NoOxy_ET_1p5V_2sec.Ca(end-1, tIndexCa)], 1);];

%%%% RBC data - Missing Exp from 1875 - 9 mice & 12 exp %%%%
tIndexVasc = d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, :) >= xL(1) ...
    & d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, :) <= xL(2);
timeVasc = d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc);

Paired.RBC.Oxy = [d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    mean([d.M1514_051120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1514_191120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1); ...
    d.M1594_261120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1804_080221.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1854_090321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1951_180321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1873_220321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    mean([d.M1952_240321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1952_310321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1); ...
    mean([d.M1954_240321.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1954_010421.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1);];

Paired.RBC.NoOxy = [d.M1397_041120.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    mean([d.M1514_051120.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1514_191120.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1); ...
    d.M1594_261120.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1804_080221.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1854_090321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1951_180321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1873_220321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    mean([d.M1952_240321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1952_310321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1); ...
    mean([d.M1954_240321.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc); ...
    d.M1954_010421.(typeVasc).NoOxy_ET_1p5V_2sec.RBC(end-1, tIndexVasc)], 1);];

%%%% PO2All data - 10 mice & 13 exp %%%%
Paired.PO2All.Oxy = [d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    mean([d.M1514_051120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1514_191120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1); ...
    d.M1594_261120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1804_080221.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1875_160221.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1854_090321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1951_180321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1873_220321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    mean([d.M1952_240321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1952_310321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1); ...
    mean([d.M1954_240321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1954_010421.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1);];

Paired.PO2All.NoOxy = [d.M1397_041120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    mean([d.M1514_051120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1514_191120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1); ...
    d.M1594_261120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1804_080221.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1875_160221.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1854_090321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1951_180321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1873_220321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    mean([d.M1952_240321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1952_310321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1); ...
    mean([d.M1954_240321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1954_010421.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1);];

%%%% Flow data - 6 mice & 8 exp %%%%
Paired.Flow.Oxy = [mean([d.M1514_051120.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    d.M1514_191120.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc)], 1); ...
    d.M1804_080221.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    d.M1854_090321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    d.M1873_220321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    d.M1952_240321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    mean([d.M1954_240321.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    d.M1954_010421.(typeVasc).Oxygen_ET_1p5V_2sec.Flow(end-1, tIndexVasc)], 1);];

Paired.Flow.NoOxy = [mean([d.M1514_051120.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    d.M1514_191120.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc)], 1); ...
    d.M1804_080221.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    d.M1854_090321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    d.M1873_220321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    d.M1952_240321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    mean([d.M1954_240321.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc); ...
    d.M1954_010421.(typeVasc).NoOxy_ET_1p5V_2sec.Flow(end-1, tIndexVasc)], 1);];

figure;
% Calcium data
subplot(221); hold on;
yL = [-100 500];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
x = timeCa;
y = mean(Paired.Ca.Oxy, 1);
ystd = std(Paired.Ca.Oxy, 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{1}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{1}, 'LineWidth', 2);

y = mean(Paired.Ca.NoOxy, 1);
ystd = std(Paired.Ca.NoOxy, 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{2}, 'LineWidth', 2);
set(gca, 'FontSize', 13);
ylabel('\Delta Ca^{2+} (a.u.)');
xlim(xL); ylim(yL);

% RBC velocity data
subplot(222); hold on;
yL = [-0.2 0.8];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
x = timeVasc;
y = mean(Paired.RBC.Oxy, 1);
ystd = std(Paired.RBC.Oxy, 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{1}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{1}, 'LineWidth', 2);

y = mean(Paired.RBC.NoOxy, 1);
ystd = std(Paired.RBC.NoOxy, 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{2}, 'LineWidth', 2);
ylabel('\Delta RBC Velocity (mm.s^{-1})');
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

% Flow data
subplot(223); hold on;
yL = [-10 50];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
y = mean(Paired.Flow.Oxy, 1);
ystd = std(Paired.Flow.Oxy, 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{1}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{1}, 'LineWidth', 2);

y = mean(Paired.Flow.NoOxy, 1);
ystd = std(Paired.Flow.NoOxy, 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{2}, 'LineWidth', 2);
ylabel('\Delta Flow (RBC.s^{-1})'); xlabel('Time (s)');
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

% PO2 data
subplot(224); hold on;
yL = [-12 30];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
y = mean(Paired.PO2All.Oxy, 1);
ystd = std(Paired.PO2All.Oxy, 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{1}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{1}, 'LineWidth', 2);

y = mean(Paired.PO2All.NoOxy, 1);
ystd = std(Paired.PO2All.NoOxy, 0, 1);
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], colors{2}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, y, 'Color', colors{2}, 'LineWidth', 2);
ylabel('\Delta pO_{2} Mean (mmHg)'); xlabel('Time (s)');
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

suptitle('FIGURE 4F');

%% Panel G
d = data;
xL = [0 29];

tIndexVasc = d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, :) >= xL(1) ...
    & d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, :) <= xL(2);
timeVasc = d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.RBC(1, tIndexVasc);
typeVasc = 'ZScore';

%%%% PO2All data - 10 mice & 13 exp %%%%
Paired.PO2All.Oxy = [d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    mean([d.M1514_051120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1514_191120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1); ...
    d.M1594_261120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1804_080221.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1875_160221.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1854_090321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1951_180321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1873_220321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    mean([d.M1952_240321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1952_310321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1); ...
    mean([d.M1954_240321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1954_010421.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1);];

Paired.PO2All.NoOxy = [d.M1397_041120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    mean([d.M1514_051120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1514_191120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1); ...
    d.M1594_261120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1804_080221.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1875_160221.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1854_090321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1951_180321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1873_220321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    mean([d.M1952_240321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1952_310321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1); ...
    mean([d.M1954_240321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc); ...
    d.M1954_010421.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(end-1, tIndexVasc)], 1);];

Oxy = {d.M1397_041120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    [d.M1514_051120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1514_191120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc)]; ...
    d.M1594_261120.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1804_080221.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1875_160221.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1854_090321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1951_180321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1873_220321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    [d.M1952_240321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1952_310321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc)]; ...
    [d.M1954_240321.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1954_010421.(typeVasc).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc)]; ...
    Paired.PO2All.Oxy};

NoOxy = {d.M1397_041120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    [d.M1514_051120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1514_191120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc)]; ...
    d.M1594_261120.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1804_080221.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1875_160221.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1854_090321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1951_180321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1873_220321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    [d.M1952_240321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1952_310321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc)]; ...
    [d.M1954_240321.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc); ...
    d.M1954_010421.(typeVasc).NoOxy_ET_1p5V_2sec.PO2All(2:end-2, tIndexVasc)]; ...
    Paired.PO2All.NoOxy};


% Bootstrap test
bslFrame = [1 10];
dipFrame = [10.5 11.5];
nboot = 10000;
bt_bsl = struct(); pval = struct();
rng('default'); % For reproduciblity

for i=1:length(Oxy)
    mat.Oxy = Oxy{i}; mat.NoOxy = NoOxy{i};

    bt_bsl.Oxy(:, i) = matrixBootstrap(nboot, @mean, mat.Oxy(:, timeVasc > bslFrame(1) & timeVasc < bslFrame(2)));
    bt_bsl.NoOxy(:, i) = matrixBootstrap(nboot, @mean, mat.NoOxy(:, timeVasc > bslFrame(1) & timeVasc < bslFrame(2)));
end

alpha = 0.05;

figure;
subplot(221); hold on;
yL = [0 1000];
m = mean(mean(Oxy{end}(:, timeVasc > dipFrame(1) & timeVasc < dipFrame(2)), 2));
thresh = prctile(bt_bsl.Oxy(:, end), alpha/2*100);
threshUp = prctile(bt_bsl.Oxy(:, end), 100*(1-alpha/2));
patch('XData', [thresh threshUp threshUp thresh], 'YData', ...
    [yL(1) yL(1) yL(2) yL(2)], 'FaceColor', colorStim, 'EdgeColor', ...
    'none', 'FaceAlpha', 0.5);
line([m m], ylim, 'Color', colors{5}, 'LineWidth', 5);
histogram(bt_bsl.Oxy(:, end), [-1:0.01:1], 'LineStyle', 'none', ...
    'FaceColor', colors{3});
ylabel('Count'); xlabel('Average Z-score (SD)');
xlim([-0.5 0.5]); ylim([0 1000]);
set(gca, 'FontSize', 13);

pval.Oxy = sum(bt_bsl.Oxy(:, end) < m)/nboot;


subplot(222); hold on;
m = mean(mean(NoOxy{end}(:, timeVasc > dipFrame(1) & timeVasc < dipFrame(2)), 2));
thresh = prctile(bt_bsl.NoOxy(:, end), alpha/2*100);
threshUp = prctile(bt_bsl.NoOxy(:, end), 100*(1-alpha/2));
patch('XData', [thresh threshUp threshUp thresh], 'YData', ...
    [yL(1) yL(1) yL(2) yL(2)], 'FaceColor', colorStim, 'EdgeColor', ...
    'none', 'FaceAlpha', 0.5, 'DisplayName', 'CI 95%');
line([m m], ylim, 'Color', colors{5}, 'LineWidth', 5, 'DisplayName', ...
    'Initial Dip Period', 'Color', colors{5});
histogram(bt_bsl.NoOxy(:, end), [-1:0.01:1], 'DisplayName', ...
    'Baseline Distribution', 'LineStyle', 'none', 'FaceColor', colors{3});
xlim([-0.5 0.5]); ylim([0 1000]);
xlabel('Average Z-score (SD)');
set(gca, 'FontSize', 13);

pval.NoOxy = sum(bt_bsl.NoOxy(:, end) < m)/nboot;
disp('### Figure 4F ###');
disp(['Significant dip in both Oxygenated and reduced O2 mice with p = ' ...
    num2str(pval.Oxy) ' and ' num2str(pval.NoOxy) ' respectively.']);

xL = [0 29]; yL = [-3 12];
subplot(223); hold on;
patch('XData', [dipFrame(1) dipFrame(2) dipFrame(2) dipFrame(1)], ...
    'YData', [yL(1) yL(1) yL(2) yL(2)], 'FaceColor', colors{5}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
patch('XData', [bslFrame(1) bslFrame(2) bslFrame(2) bslFrame(1)], ...
    'YData', [yL(1) yL(1) yL(2) yL(2)], 'FaceColor', colors{3}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(timeVasc, Paired.PO2All.Oxy, 'Color', [0.5 0.5 0.5]);
plot(timeVasc, mean(Paired.PO2All.Oxy, 1), 'Color', colors{1}, ...
    'LineWidth', 2);
xlim(xL); ylim(yL);
ylabel('\Delta Po_{2} Mean (SD)'); xlabel('Time (s)');
set(gca, 'FontSize', 13);

subplot(224); hold on;
patch('XData', [dipFrame(1) dipFrame(2) dipFrame(2) dipFrame(1)], ...
    'YData', [yL(1) yL(1) yL(2) yL(2)], 'FaceColor', colors{5}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
patch('XData', [bslFrame(1) bslFrame(2) bslFrame(2) bslFrame(1)], ...
    'YData', [yL(1) yL(1) yL(2) yL(2)], 'FaceColor', colors{3}, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(timeVasc, Paired.PO2All.NoOxy, 'Color', [0.5 0.5 0.5]);
plot(timeVasc, mean(Paired.PO2All.NoOxy, 1), 'Color', colors{2}, ...
    'LineWidth', 2);
xlim(xL); ylim(yL);
xlabel('Time (s)');
set(gca, 'FontSize', 13);


suptitle('FIGURE 4G');


end
