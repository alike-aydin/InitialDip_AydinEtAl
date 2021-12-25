function Fig1()
% FIG1 Generating Figure 1 from Aydin et al.
%
% function Fig1() = []
%
%   Author : Ali-Kemal Aydin, PhD student
%   Date : December 30th, 2020
%   Mail: ali-kemal.aydin@inserm.fr
%   Affiliation : U968, Institut de la Vision, Paris
%   License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
%       See LICENSE.txt or <a href="matlab:web('https://creativecommons.org/licenses/by/4.0/')">here</a> 
%       for a human-readable version.
%
%   DESCRIPTION : Generates the panels from Figure 1 in Aydin et al.
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

%% Panel A
% Image on the right is not shown by this script

cd = 'Data_Movie_Fig1A/';

ET200mV = dlmread(fullfile(cd, 'ET_200mV_2sec_deltaCa_MSG8.txt'), '\t', 1, 0);
ET200mV = ET200mV(:, [1:10 13:20]);% Columns 11 & 12 are an empty ROI (full of 0s)
MSG_ROI.ET200mV = 4; % Most sensitive glomerulus is in the 4th ROI of this analysis

ET4mV = dlmread(fullfile(cd, 'ET_4mV_2sec_deltaCa_MSG8.txt'), '\t', 1, 0);
ET4mV = ET4mV(:, [1:8 11:20]);% Columns 9 & 10 are an empty ROI (full of 0s)
MSG_ROI.ET4mV = 3; % Most sensitive glomerulus is in the 3rd ROI of this analysis

yL = [-50 100]; xL = [2 10];

figure;
subplot(211); title('ET 1% 2 s'); hold on;
patch('XData', [5 7 7 5], 'YData', ...
    [yL(1) yL(1) yL(2) yL(2)], 'FaceColor', colorStim, 'EdgeColor', ...
    'none', 'FaceAlpha', 0.5, 'DisplayName', 'Stimulus');

for i=1:9
    if i~= MSG_ROI.ET200mV
        plot(ET200mV(:, i*2-1), ET200mV(:, i*2), 'Color', 'k', ...
            'LineWidth', 1);
    end
end

plot(ET200mV(:, MSG_ROI.ET200mV*2-1), ET200mV(:, MSG_ROI.ET200mV*2), ...
    'Color', colors{1}, 'LineWidth', 2);
ylabel('\Delta Ca^{2+} (a.u.)'); 
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13);

subplot(212); title('ET 0.02% 2 s'); hold on;
patch('XData', [5 7 7 5], 'YData', ...
    [yL(1) yL(1) yL(2) yL(2)], 'FaceColor', ...
    colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', ...
    'Stimulus');

for i=1:9
    if i~= MSG_ROI.ET4mV
        plot(ET4mV(:, i*2-1), ET4mV(:, i*2), 'Color', 'k', ...
            'LineWidth', 1);
    end
end
plot(ET4mV(:, MSG_ROI.ET4mV*2-1), ET4mV(:, MSG_ROI.ET4mV*2), ...
    'Color', colors{1}, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta Ca^{2+} (a.u.)'); 
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13);

suptitle('FIGURE 1A');

%% Panel B

m = 'M1393_021020';
type = 'Delta';
xL = [0 29];
yLCa = [-100 500]; yLRBC = [-0.2 0.35];

CaTS = struct(); RBCTS = struct();
CaTS.ET_1p5V_2sec = d.(m).(type).Oxygen_ET_1p5V_2sec.Ca; % n = 2 acq
CaTS.ET_200mV_2sec = d.(m).(type).Oxygen_ET_200mV_2sec.Ca; % n = 2 acq
CaTS.ET_20mV_2sec = d.(m).(type).Oxygen_ET_20mV_2sec.Ca; % n = 2 acq

RBCTS.ET_1p5V_2sec = d.(m).(type).Oxygen_ET_1p5V_2sec.RBC;
RBCTS.ET_200mV_2sec = d.(m).(type).Oxygen_ET_200mV_2sec.RBC;
RBCTS.ET_20mV_2sec = d.(m).(type).Oxygen_ET_20mV_2sec.RBC;


figure;
subplot(211); hold on;
patch('XData', [10 12 12 10], 'YData', ...
    [yLCa(1) yLCa(1) yLCa(2) yLCa(2)], 'FaceColor', ...
    colorStim, 'EdgeColor', 'none','FaceAlpha', 0.5, 'DisplayName', ...
    'Stimulus');

legend('AutoUpdate', 'off'); % To avoid having one legend per acquisition
plot(CaTS.ET_20mV_2sec(1, :), CaTS.ET_20mV_2sec(2:end-2, :), ...
    'Color', colors{3}*0.8); legend('AutoUpdate', 'on');
plot(CaTS.ET_20mV_2sec(1, :), CaTS.ET_20mV_2sec(end-1, :), ...
    'Color', colors{3}, 'DisplayName', '0.1%', 'LineWidth', 2);

legend('AutoUpdate', 'off');
plot(CaTS.ET_200mV_2sec(1, :), CaTS.ET_200mV_2sec(2:end-2, :), ...
    'Color', colors{2}*0.8); legend('AutoUpdate', 'on');
plot(CaTS.ET_200mV_2sec(1, :), CaTS.ET_200mV_2sec(end-1, :), ...
    'Color', colors{2}, 'DisplayName', '1%', 'LineWidth', 2);

legend('AutoUpdate', 'off');
plot(CaTS.ET_1p5V_2sec(1, :), CaTS.ET_1p5V_2sec(2:end-2, :), ...
    'Color', colors{1}*0.8); legend('AutoUpdate', 'on');
plot(CaTS.ET_1p5V_2sec(1, :), CaTS.ET_1p5V_2sec(end-1, :), ...
    'Color', colors{1}, 'DisplayName', '6%', 'LineWidth', 2);

xlim(xL); ylim(yLCa); ylabel('\Delta Ca^{2+} (a.u.)');
xlabel('Time (s)'); legend(); set(gca, 'FontSize', 13);
hold off;

subplot(212); hold on;
patch('XData', [10 12 12 10], 'YData', ...
    [yLRBC(1) yLRBC(1) yLRBC(2) yLRBC(2)], 'FaceColor', ...
    colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

plot(RBCTS.ET_20mV_2sec(1, :), RBCTS.ET_20mV_2sec(2:end-2, :), ...
    'Color', colors{3}*0.8);
plot(RBCTS.ET_20mV_2sec(1, :), RBCTS.ET_20mV_2sec(end-1, :), ...
    'Color', colors{3}, 'LineWidth', 2);
plot(RBCTS.ET_200mV_2sec(1, :), RBCTS.ET_200mV_2sec(2:end-2, :), ...
    'Color', colors{2}*0.8);
plot(RBCTS.ET_200mV_2sec(1, :), RBCTS.ET_200mV_2sec(end-1, :), ...
    'Color', colors{2}, 'LineWidth', 2);
plot(RBCTS.ET_1p5V_2sec(1, :), RBCTS.ET_1p5V_2sec(2:end-2, :), ...
    'Color', colors{1}*0.8);
plot(RBCTS.ET_1p5V_2sec(1, :), RBCTS.ET_1p5V_2sec(end-1, :), ...
    'Color', colors{1}, 'LineWidth', 2);

xlim(xL); ylim(yLRBC); xlabel('Time (s)');
ylabel('\Delta RBC velocity (mm.s^{-1})');
set(gca, 'FontSize', 13);
hold off;

suptitle('FIGURE 1B');

%% Panel C
% This panel reproduces the treatment done on PointScan to extract a Ca2+
% signal from the green channel, and thus be able to align pointscan to
% neuronal activation. For more details on the process, see the Methods.
% The cartoon on the left is not reproduced here.

% Data from M1804_220221
fldr = 'Data_PointScan_1C\';

acqs = {'d07_P1Pt_ET_1p5V_2sec_1'; ...
    'd08_P1Pt_ET_1p5V_2sec_2'; ...
    'd09_P1Pt_ET_1p5V_2sec_3'; ...
    'd10_P1Pt_ET_1p5V_2sec_4'};

alignTime = [.3282 0.3080 0.4742 0.4056]; 


dt = 312/1250e3; thresh = 0.3;
winCenter = round(30e-3/dt);

Ca = {};
CaCh2 = {};
Ca_RBCPart = {};
Ca_RBCPartCh2 = {};
InterpCa = []; InterpCaCh2 = [];
DeltaCa = []; DeltaCaCh2 = [];
NormCa = []; NormCaCh2 = [];
AlignedCaRBCPart = {}; AlignedCaRBCPartCh2 = {};

timeInterp = [0:0.05:29];

% For each acquisition
for i=1:length(acqs)
    % Getting the raw signal from both channels during the ON period of
    % the linescan
    Ca{i} = getCaFromPtScan([fldr acqs{i}]);
    CaCh2{i} = getCaFromPtScan([fldr acqs{i}], false, 0, 1);
    Ca_RBCPart{i} = []; Ca_RBCPart_Ch2{i} = [];
    
    % Each chunk of signal is treated as to get the mean bottom 30%, to
    % focus on the part where there is RBC passing by, so as to not get
    % polluted by plasmatic fluorescence from FITC
    for j=1:floor(length(Ca{i}(2, :))/winCenter)
        dataSample = Ca{i}(2, (j-1)*winCenter+1:j*winCenter);
        timeSample = Ca{i}(1, (j-1)*winCenter+1:j*winCenter);
        
        range = max(dataSample)- min(dataSample);
        belowMinThreshIdx = find(dataSample < min(dataSample) + thresh * range);
        
        Ca_RBCPart{i} = [Ca_RBCPart{i} [timeSample(belowMinThreshIdx); ...
            dataSample(belowMinThreshIdx)]];
    end
    
    % Same treatment for the red channel, for comparison purposes
    for j=1:floor(length(CaCh2{i}(2, :))/winCenter)
        dataSample = CaCh2{i}(2, (j-1)*winCenter+1:j*winCenter);
        timeSample = CaCh2{i}(1, (j-1)*winCenter+1:j*winCenter);
        
        range = max(dataSample)- min(dataSample);
        belowMinThreshIdx = find(dataSample < min(dataSample) + thresh * range);
        
        Ca_RBCPart_Ch2{i} = [Ca_RBCPart_Ch2{i} [timeSample(belowMinThreshIdx); ...
            dataSample(belowMinThreshIdx)]];
    end
    
    % Aligning and normalizing
    Ca_RBCPart{i}(2, :) = sgolayfilt(Ca_RBCPart{i}(2, :), 3, 599);
    AlignedCaRBCPart{i} = [Ca_RBCPart{i}(1, :) - alignTime(i); Ca_RBCPart{i}(2, :)];
    InterpCa(i, :) = interp1(AlignedCaRBCPart{i}(1, :), AlignedCaRBCPart{i}(2, :), timeInterp, 'pchip');
    DeltaCa(i, :) = InterpCa(i, :) - mean(InterpCa(i, timeInterp < 10));
    NormCa(i, :) = DeltaCa(i, :)/max(DeltaCa(i, :));
    
    Ca_RBCPart_Ch2{i}(2, :) = sgolayfilt(Ca_RBCPart_Ch2{i}(2, :), 3, 599);
    AlignedCaCh2{i} = [Ca_RBCPart_Ch2{i}(1, :) - alignTime(i); Ca_RBCPart_Ch2{i}(2, :)];
    InterpCaCh2(i, :) = interp1(AlignedCaCh2{i}(1, :), AlignedCaCh2{i}(2, :), timeInterp, 'pchip');
    DeltaCaCh2(i, :) = InterpCaCh2(i, :) - mean(InterpCaCh2(i, timeInterp < 10));
    NormCaCh2(i, :) = DeltaCaCh2(i, :)/max(DeltaCaCh2(i, :));
end


figure; 
subplot(131); title('Green Channel'); hold on;
patch('XData', [10 12 12 10], 'YData', [-1 -1 2 2], 'FaceColor', ...
    colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(Ca{3}(1, :), Ca{3}(2, :), 'Color', colors{1});
plot(Ca_RBCPart{3}(1, :), Ca_RBCPart{3}(2, :), 'Color', colors{2}, ...
    'LineWidth', 2);
ylim([0.3 1.3]); xlim([5 15]); yticks({});
xlabel('Time (s)'); ylabel('Photon Count (a.u.)');
set(gca, 'FontSize', 13);

subplot(132); title('Red Channel'); hold on;
patch('XData', [10 12 12 10], 'YData', ...
    [3.1635e4 3.1635e4 3.4635e4 3.4635e4], 'FaceColor', colorStim, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(CaCh2{3}(1, :), CaCh2{3}(2, :), 'Color', colors{1});
plot(Ca_RBCPart_Ch2{3}(1, :), Ca_RBCPart_Ch2{3}(2, :), 'Color', ...
    colors{4}, 'LineWidth', 2);
xlim([5 15]); ylim([3.1635e4 3.4635e4]); yticks({});
xlabel('Time (s)'); 
set(gca, 'FontSize', 13);

subplot(133); title({'Aligned average'; 'neural response'}); hold on;
xL = [5 15];
patch('XData', [10 12 12 10], 'YData', [-0.5 -0.5 1.5 1.5], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');

x = d.M1804_220221.Delta.Oxygen_ET_1p5V_2sec.Ca(1, :);
y = d.M1804_220221.Delta.Oxygen_ET_1p5V_2sec.Ca(end-1, :); 
ystd = d.M1804_220221.Delta.Oxygen_ET_1p5V_2sec.Ca(end, :)/max(y);
y = y/max(y); Y = y(x >= xL(1) & x <= xL(2));
legend('AutoUpdate', 'off');
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], 'k', ...
    'edgecolor', 'none', 'facealpha', 0.5);
legend('AutoUpdate', 'on');
plot(x, y, 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Linescan Ca^{2+}');

x = timeInterp;
y = mean(NormCa, 1); ystd = std(NormCa, 0, 1); 
legend('AutoUpdate', 'off');
patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], ...
    colors{2}, 'edgecolor', 'none', 'facealpha', 0.5);
legend('AutoUpdate', 'on');
plot(x, y, 'Color', colors{2}, 'LineWidth', 2, ...
    'DisplayName', 'PointScan Ca^{2+}');

xlim(xL); ylim([-0.5 1.5]); legend();
xlabel('Time (s)'); ylabel('Normalized Signal');
set(gca, 'FontSize', 13);

suptitle('FIGURE 1C');

y = y(x >= xL(1) & x <= xL(2));
ssresid = sum((Y - y).^2);
sstotal = (length(Y)-1) * var(Y);
rsq = 1 - ssresid/sstotal;
disp('### Figure 1C ###');
disp(['R² between normalized average linescan and normalized average ' ...
    'pointscan traces: ' num2str(rsq) '.']);

%% Panel D

m = 'M1393_021020';
type = 'Interp';
yL = [0 80]; xL = [0 29];

PO2TS = struct();
PO2TS.ET_1p5V_2sec = d.(m).(type).Oxygen_ET_1p5V_2sec.PO2All;
PO2TS.ET_200mV_2sec = d.(m).(type).Oxygen_ET_200mV_2sec.PO2All;
PO2TS.ET_20mV_2sec = d.(m).(type).Oxygen_ET_20mV_2sec.PO2All;


figure;
% There is no legend in this figure, this is the same color code as 1B.

subplot(131); hold on; title('ET 0.1% 2 s');
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(PO2TS.ET_20mV_2sec(1, :), PO2TS.ET_20mV_2sec(2:end-2, :), ...
    'Color', [0.5 0.5 0.5]);
plot(PO2TS.ET_20mV_2sec(1, :), PO2TS.ET_20mV_2sec(end-1, :), ...
    'Color', colors{3}, 'LineWidth', 2);
ylabel('Po_{2} (mmHg)'); xlabel('Time (s)');
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13);

subplot(132); hold on; title('ET 1% 2 s');
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(PO2TS.ET_200mV_2sec(1, :), PO2TS.ET_200mV_2sec(2:end-2, :), ...
    'Color', [0.5 0.5 0.5]);
plot(PO2TS.ET_200mV_2sec(1, :), PO2TS.ET_200mV_2sec(end-1, :), ...
    'Color', colors{2}, 'LineWidth', 2);
xlabel('Time (s)');
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13);

subplot(133); hold on; title('ET 6% 2 s');
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(PO2TS.ET_1p5V_2sec(1, :), PO2TS.ET_1p5V_2sec(2:end-2, :), ...
    'Color', [0.5 0.5 0.5]);
plot(PO2TS.ET_1p5V_2sec(1, :), PO2TS.ET_1p5V_2sec(end-1, :), ...
    'Color', colors{1}, 'LineWidth', 2);
xlim(xL); ylim(yL);
xlabel('Time (s)');
set(gca, 'FontSize', 13);
hold off;

suptitle('FIGURE 1D');

%% Panel E

type = 'Delta';
yL = [-15 60]; xL = [0 29];

tIndexPO2 = d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, :) >= xL(1) ...
    & d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, :) <= xL(2);
timePO2 = d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2);
exps = struct();

% Gathering the average response for each mouse at each concentration

% n = 14 mice, for 17 experiments
exps.ET1p5V = [d.M1393_021020.(type).Oxygen_ET_1p5V_2sec.PO2All(end-1, tIndexPO2); ...
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

% n = 4 mice, for 4 experiments
exps.ET200mV = [d.M1393_021020.(type).Oxygen_ET_200mV_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1397_131020.(type).Oxygen_ET_200mV_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1509_261020.(type).Oxygen_ET_200mV_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1511_211020.(type).Oxygen_ET_200mV_2sec.PO2All(end-1, tIndexPO2)];

% n = 6 mice, for 6 experiments
exps.ET20mV = [d.M1393_021020.(type).Oxygen_ET_20mV_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1397_131020.(type).Oxygen_ET_20mV_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1509_261020.(type).Oxygen_ET_20mV_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1511_211020.(type).Oxygen_ET_20mV_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1514_051120.(type).Oxygen_ET_20mV_2sec.PO2All(end-1, tIndexPO2); ...
    d.M1875_160221.(type).Oxygen_ET_20mV_2sec.PO2All(end-1, tIndexPO2)];

figure;
% There is no legend in this figure, this is the same color code as 1B.

subplot(131); hold on; title('ET 0.1% 2 s');
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
h = plot(timePO2, exps.ET20mV, 'Color', [0.5 0.5 0.5]);
% Putting the mouse shown in B & C in black and in front of the other
% traces
set(h(1), 'Color', 'k'); uistack(h(1), 'top');
plot(timePO2, mean(exps.ET20mV, 1), 'Color', colors{3}, 'LineWidth', 2);
ylabel('\Delta Po_{2} (mmHg)'); xlabel('Time (s)');
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13); 

subplot(132); hold on; title('ET 1% 2 s');
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
h = plot(timePO2, exps.ET200mV, 'Color', [0.5 0.5 0.5]);
set(h(1), 'Color', 'k'); uistack(h(1), 'top');
plot(timePO2, mean(exps.ET200mV, 1), ...
    'Color', colors{2}, 'LineWidth', 2);
xlabel('Time (s)');
ylim(yL); xlim(xL);
set(gca, 'FontSize', 13);

subplot(133); hold on; title('ET 6% 2 s');
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
h = plot(timePO2, exps.ET1p5V, 'Color', [0.5 0.5 0.5]);
set(h(1), 'Color', 'k'); uistack(h(1), 'top');
plot(timePO2, mean(exps.ET1p5V, 1), ...
    'Color', colors{1}, 'LineWidth', 2);
xlabel('Time (s)');
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);
hold off;

suptitle('Figure 1E');

%% Panel F

type = 'Delta';
xL = [0 29]; yL = [-15 50];
tIndexPO2 = d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, :) >= xL(1) ...
    & d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, :) <= xL(2);
timePO2 = d.M1397_041120.(type).Oxygen_ET_1p5V_2sec.PO2All(1, tIndexPO2);
deltaT = 0.5; % Binning time window in second
exps = struct();

% Gathering the single responses for each mouse at each concentration

% n = 14 mice, for 17 experiments and multiple acquisitions per mouse
exps.ET1p5V = {d.M1393_021020.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2); ...
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
    d.M1954_010421.(type).Oxygen_ET_1p5V_2sec.PO2All(2:end-2, tIndexPO2)]};

% n = 4 mice, for 4 experiments and multiple acquisitions per mouse
exps.ET200mV = {d.M1393_021020.(type).Oxygen_ET_200mV_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1397_131020.(type).Oxygen_ET_200mV_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1509_261020.(type).Oxygen_ET_200mV_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1511_211020.(type).Oxygen_ET_200mV_2sec.PO2All(2:end-2, tIndexPO2)};

% n = 6 mice, for 6 experiments and multiple acquisitions per mouse
exps.ET20mV = {d.M1393_021020.(type).Oxygen_ET_20mV_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1397_131020.(type).Oxygen_ET_20mV_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1509_261020.(type).Oxygen_ET_20mV_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1511_211020.(type).Oxygen_ET_20mV_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1514_051120.(type).Oxygen_ET_20mV_2sec.PO2All(2:end-2, tIndexPO2); ...
    d.M1875_160221.(type).Oxygen_ET_20mV_2sec.PO2All(2:end-2, tIndexPO2)};

integ = struct();
conds = fieldnames(exps);

% Computing the average AUC for each bin, as the average over all mice of
% the average AUC across acquisitions
for i=1:length(conds)
    integ.(conds{i}) = [];  
    
    for j=1:length(exps.(conds{i}))
        mat = exps.(conds{i}){j, 1};
        tmp = [];
        for k=1:size(mat, 1)
            for t=1:floor((timePO2(end)-timePO2(1))/deltaT)
                tmp(k, t) = sum(mat(k, timePO2 >= timePO2(1)+deltaT*(t-1) ...
                    & timePO2 < timePO2(1)+deltaT*t));
            end
        end
        
        integ.(conds{i})(j, :) = mean(tmp, 1);
    end
end

timepoints = timePO2(1)+deltaT/2:deltaT:timePO2(end)-deltaT/2;

figure;
% There is no legend in this figure, this is the same color code as 1B.

subplot(131); hold on; title('ET 0.1% 2 s');
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
h = bar(timepoints, mean(integ.ET20mV, 1), 'FaceColor', 'flat');
for i=1:size(integ.ET20mV, 2)
    h.CData(i, :) = colors{3};
end

xlim(xL); ylim(yL);
ylabel('AUC \Delta PO_{2} (mmHg.s)'); xlabel('Time (s)');
set(gca, 'FontSize', 13);

subplot(132); hold on; title('ET 1% 2 s');
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
h = bar(timepoints, mean(integ.ET200mV, 1), 'FaceColor', 'flat');
for i=1:size(integ.ET200mV, 2)
    h.CData(i, :) = colors{2};
end

xlabel('Time (s)');
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

subplot(133); hold on; title('ET 6% 2 s');
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
h = bar(timepoints, mean(integ.ET1p5V, 1), 'FaceColor', 'flat');
for i=1:size(integ.ET1p5V, 2)
    h.CData(i, :) = colors{1};
end

xlabel('Time (s)'); 
xlim(xL); ylim(yL);
set(gca, 'FontSize', 13);

suptitle('FIGURE 1F');

end

