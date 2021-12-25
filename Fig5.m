function Fig5()
% FIG5 Generates the panels of the Figure 5 of Aydin et al.
%
% function Fig5() = []
%
%   Author: Ali-Kemal Aydin, PhD student & Camille Verdier, Master's intern
%   Mail: ali-kemal.aydin@inserm.fr
%   Affiliations: 
%       * INSERM, CNRS, Institut de la Vision, Sorbonne Université, Paris, France
%   License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
%       See LICENSE.txt or <a href="matlab:web('https://creativecommons.org/licenses/by/4.0/')">here</a>
%       for a human-readable version.
%
%   DESCRIPTION: Generates the panels from Figure 5 in Aydin et al.
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
load('awakeData.mat');

%% Panel A
% Individual glom : 1951 - Glom 3
xL = [1 25];

figure; subplot(221); hold on;
d = data.M1951_200521.Delta.Glom3.Ca; yL = [-20 100];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');
plot(d(1, :), d(2:end-2, :), 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
plot(d(1, :), d(end-1, :), 'LineWidth', 2, 'Color', [0.20 0.70 0.32]);
xlabel('Time (s)'); ylabel('\Delta Ca^{2+}');
set(gca, 'FontSize', 13); xlim(xL); ylim(yL);

subplot(223); hold on;
d = data.M1951_200521.Interp.Glom3.RBC; yL = [0 1.5];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');
plot(d(1, :), d(2:end-2, :), 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
plot(d(1, :), d(end-1, :), 'LineWidth', 2, 'Color', 'k');
xlabel('Time (s)'); ylabel('RBC velocity (mm.s^{-1})');
set(gca, 'FontSize', 13); xlim(xL); ylim(yL);

subplot(222); hold on;
d = fluxfiltre.m1951.capi2.d20210520.B6_10p; yL = [0 100];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');
allFlow = [d.d13(:, 2) d.d14(:, 2) d.d18(:, 2)]; 
plot(d.d13(:, 1), allFlow, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
plot(d.d13(:, 1), mean(allFlow, 2, 'omitnan'), 'LineWidth', 2, 'Color', 'k');
xlabel('Time (s)'); ylabel('RBC flow (RBC.s^{-1})'); 
set(gca, 'FontSize', 13); xlim(xL); ylim(yL);

subplot(224); hold on;
d = po2filtre.m1951.capi2.d20210520.B6_10p; yL = [0 60];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');
allPO2 = [d.d13(:, 2) d.d14(:, 2) d.d18(:, 2)];
plot(d.d13(:, 1), allPO2, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
plot(d.d13(:, 1), mean(allPO2, 2, 'omitnan'), 'LineWidth', 2, 'Color', 'k');
xlabel('Time (s)'); ylabel('Mean Po_{2} (mmHg)');
set(gca, 'FontSize', 13); xlim(xL); ylim(yL);

suptitle('FIGURE 6A');

%% Panel B
% Mean of all, paired or not

% Data preparation
means = struct();
means.vessels = [];
means.tissues = [];
msg.vessels = [];
msg.tissues = [];
means.flows = [];

zPO2Vessel = zscoreData(po2filtre);
zPO2Tissue = zscoreData(tissue_po2filtre);
zFlow = zscoreData(fluxfiltre);

% Vessels Not MSG
manip1951_G1 = zPO2Vessel.m1951.capi1.d20210506AM.B6_10p; % 4 FLOW
means.vessels(:, 1) = mean([manip1951_G1.d11(:, 2) manip1951_G1.d12(:, 2) manip1951_G1.d13(:, 2) manip1951_G1.d16(:, 2)], 2, 'omitnan');
manip1951_G2 = zPO2Vessel.m1951.capi2.d20210506PM.B6_20p; % 1 FLOW
means.vessels(:, 2) = mean([manip1951_G2.d14(:, 2) manip1951_G2.d15(:, 2) manip1951_G2.d16(:, 2) manip1951_G2.d18(:, 2)], 2, 'omitnan');
manip1951_G3 = zPO2Vessel.m1951.capi2.d20210520.B6_10p; % 3 FLOW
means.vessels(:, 3) = mean([manip1951_G3.d13(:, 2) manip1951_G3.d14(:, 2) manip1951_G3.d18(:, 2)], 2, 'omitnan');
manip1951_G4 = zPO2Vessel.m1951.capi4.d20210709.B6_20p; % 0 FLOW
means.vessels(:, 4) = mean([manip1951_G4.d13(:, 2) manip1951_G4.d15(:, 2) manip1951_G4.d16(:, 2)  manip1951_G4.d17(:, 2)  manip1951_G4.d18(:, 2)  manip1951_G4.d19(:, 2)  manip1951_G4.d20(:, 2)  manip1951_G4.d21(:, 2)], 2, 'omitnan');
means.vessels(:, 5) = data.M1951_220721_PM.ZScore.Glom5_Vessel.PO2All(end-1, :); % 0 FLOW
means.vessels(:, 6) = data.M1951_220721_PM.ZScore.Glom7_Vessel.PO2All(end-1, :); % 3 FLOW

manip1952_G1 = zPO2Vessel.m1952.capi1.d20210505AM.B6_20p; % 3 FLOW
means.vessels(:, 7) = mean([manip1952_G1.d06(:, 2) manip1952_G1.d08(:, 2) manip1952_G1.d09(:, 2) manip1952_G1.d10(:, 2) manip1952_G1.d11(:, 2) manip1952_G1.d12(:, 2)], 2, 'omitnan');
manip1952_G2 = zPO2Vessel.m1952.capi2.d20210505PM.B6_20p; % 1 FLOW
means.vessels(:, 8) = mean([manip1952_G2.d10(:, 2) manip1952_G2.d11(:, 2) manip1952_G2.d12(:, 2) manip1952_G2.d13(:, 2)], 2, 'omitnan');
manip1952_G3 = zPO2Vessel.m1952.capi1.d20210506AM.B6_20p; % 0 FLOW
means.vessels(:, 9) = mean([manip1952_G3.d09(:, 2) manip1952_G3.d10(:, 2) manip1952_G3.d11(:, 2) manip1952_G3.d12(:, 2)], 2, 'omitnan');
manip1952_G4 = zPO2Vessel.m1952.capi4.d20210702AM.B6_20p; % 6 FLOW
means.vessels(:, 10) = mean([manip1952_G4.d09(:, 2) manip1952_G4.d10(:, 2) manip1952_G4.d11(:, 2) manip1952_G4.d12(:, 2) manip1952_G4.d14(:, 2) manip1952_G4.d15(:, 2)], 2, 'omitnan');
manip1952_G5 = zPO2Vessel.m1952.capi5.d20210702PM.B6_20p; % 8 FLOW
means.vessels(:, 11) = mean([manip1952_G5.d07(:, 2) manip1952_G5.d08(:, 2) manip1952_G5.d09(:, 2) manip1952_G5.d10(:, 2) manip1952_G5.d11(:, 2) manip1952_G5.d12(:, 2) manip1952_G5.d13(:, 2) manip1952_G5.d14(:, 2)], 2, 'omitnan');
means.vessels(:, 12) = data.M1952_230721_PM.ZScore.Glom6_Vessel.PO2All(end-1, :); % 0 FLOW
means.vessels(:, 13) = data.M1952_230721_PM.ZScore.Glom7_Vessel.PO2All(end-1, :); % 0 FLOW

manip2115_G1 = zPO2Vessel.m2115.capi1.d20210713PM.B6_20p; % 7 FLOW
means.vessels(:, 14) = mean([manip2115_G1.d17(:, 2) manip2115_G1.d18(:, 2) manip2115_G1.d19(:, 2) manip2115_G1.d21(:, 2) manip2115_G1.d23(:, 2) manip2115_G1.d24(:, 2) manip2115_G1.d25(:, 2) manip2115_G1.d26(:, 2)], 2, 'omitnan');

% Vessels MSG
msg.vessels(:, 1) = data.M1951_220721_AM.ZScore.MSG_Vessel.PO2All(end-1, :); % 0 FLOW
msg.vessels(:, 2) = data.M1952_230721_AM.ZScore.MSG_Vessel.PO2All(end-1, :); % 0 FLOW

% Tissues MSG
msg.tissues(:, 1) = data.M1951_220721_AM.ZScore.MSG_Tissue.PO2All(end-1, :); 
msg.tissues(:, 2) = data.M1952_230721_AM.ZScore.MSG_Tissue.PO2All(end-1, :);

% Tissues Not MSG
tissue1951_G2 = zPO2Tissue.m1951.capi2.d20210506PM.B6_20p;
means.tissues(:, 1) = mean([tissue1951_G2.d19(:, 2) tissue1951_G2.d20(:, 2) tissue1951_G2.d23(:, 2) tissue1951_G2.d24(:, 2)], 2, 'omitnan');
tissue1951_G3 = zPO2Tissue.m1951.capi2.d20210520.B6_10p;
means.tissues(:, 2) = mean([tissue1951_G3.d24(:, 2) tissue1951_G3.d25(:, 2) tissue1951_G3.d28(:, 2) tissue1951_G3.d29(:, 2) tissue1951_G3.d30(:, 2) tissue1951_G3.d31(:, 2)], 2, 'omitnan');
means.tissues(:, 3) = data.M1951_220721_PM.ZScore.Glom6_Tissue.PO2All(end-1, :);
means.tissues(:, 4) = data.M1951_220721_PM.ZScore.Glom7_Tissue.PO2All(end-1, :);

tissue1952_G1 = zPO2Tissue.m1952.capi3.d20210511PM.B6_20p;
means.tissues(:, 5) = mean([tissue1952_G1.d03(:, 2) tissue1952_G1.d04(:, 2) tissue1952_G1.d05(:, 2) tissue1952_G1.d06(:, 2) tissue1952_G1.d08(:, 2)], 2, 'omitnan');
means.tissues(:, 6) = data.M1952_230721_PM.ZScore.Glom6_Tissue.PO2All(end-1, :);
means.tissues(:, 7) = data.M1952_230721_PM.ZScore.Glom7_Tissue.PO2All(end-1, :);

% Calcium Vessel
vesselCa1951_G2 = [];
f = fieldnames(cafiltre.m1951.capi2.d20210506PM.B6_20p);
for i =1:length(f)
    tmp = cafiltre.m1951.capi2.d20210506PM.B6_20p.(f{i});
    bsl = mean(tmp(tmp(:, 1) < 10, 2));
    tmp(:, 2) = (tmp(:, 2) - bsl)/max(tmp(:, 2) - bsl);
    vesselCa1951_G2(:, i) = interp1(tmp(:, 1), tmp(:, 2), 1:0.05:29, 'pchip');
end
means.ca_vessels(:, 1) = mean(vesselCa1951_G2, 2, 'omitnan');

vesselCa1951_G3 = [];
f = fieldnames(cafiltre.m1951.capi2.d20210520.B6_10p);
for i =1:length(f)
    tmp = cafiltre.m1951.capi2.d20210520.B6_10p.(f{i});
    bsl = mean(tmp(tmp(:, 1) < 10, 2));
    tmp(:, 2) = (tmp(:, 2) - bsl)/max(tmp(:, 2) - bsl);
    vesselCa1951_G3(:, i) = interp1(tmp(:, 1), tmp(:, 2), 1:0.05:29, 'pchip');
end
means.ca_vessels(:, 2) = mean(vesselCa1951_G3, 2, 'omitnan');

vesselCa1952_G6 = [];
f = fieldnames(cafiltre.m1952.capi2.d20210723PM.B6_20p);
for i =1:length(f)
    tmp = cafiltre.m1952.capi2.d20210723PM.B6_20p.(f{i});
    bsl = mean(tmp(tmp(:, 1) < 10, 2));
    tmp(:, 2) = (tmp(:, 2) - bsl)/max(tmp(:, 2) - bsl);
    vesselCa1952_G6(:, i) = interp1(tmp(:, 1), tmp(:, 2), 1:0.05:29, 'pchip');
end
means.ca_vessels(:, 2) = mean(vesselCa1952_G6, 2, 'omitnan');


% Calcium tissues
tissueCa1951_G2 = [];
f = fieldnames(tissuecafiltre.m1951.capi2.d20210506PM.B6_20p);
for i =1:length(f)
    tmp = tissuecafiltre.m1951.capi2.d20210506PM.B6_20p.(f{i});
    bsl = mean(tmp(tmp(:, 1) < 10, 2));
    tmp(:, 2) = (tmp(:, 2) - bsl)/max(tmp(:, 2) - bsl);
    tissueCa1951_G2(:, i) = interp1(tmp(:, 1), tmp(:, 2), 1:0.05:29, 'pchip');
end
means.ca_tissues(:, 1) = mean(tissueCa1951_G2, 2, 'omitnan');

tissueCa1951_G3 =[];
f = fieldnames(tissuecafiltre.m1951.capi2.d20210506PM.B6_20p);
for i =1:length(f)
    tmp =  tissuecafiltre.m1951.capi2.d20210506PM.B6_20p.(f{i});
    bsl = mean(tmp(tmp(:, 1) < 10, 2));
    tmp(:, 2) = (tmp(:, 2) - bsl)/max(tmp(:, 2) - bsl);
    tissueCa1951_G3(:, i) = interp1(tmp(:, 1), tmp(:, 2), 1:0.05:29, 'pchip');
end
means.ca_tissues(:, 2) = mean(tissueCa1951_G3, 2, 'omitnan');

tissueCa1952_G1 = [];
f = fieldnames(tissuecafiltre.m1952.capi3.d20210511PM.B6_20p);
for i =2:length(f) % First acq is weird
    tmp = tissuecafiltre.m1952.capi3.d20210511PM.B6_20p.(f{i});
    bsl = mean(tmp(tmp(:, 1) < 10, 2));
    tmp(:, 2) = (tmp(:, 2) - bsl)/max(tmp(:, 2) - bsl);
    tissueCa1952_G1(:, i) = interp1(tmp(:, 1), tmp(:, 2), 1:0.05:29, 'pchip');
end
means.ca_tissues(:, 3) = mean(tissueCa1952_G1, 2, 'omitnan');


tissueCa1952_MSG = [];
f = fieldnames(tissuecafiltre.m1952.capi1.d20210721AM.B6_20p);
for i =1:length(f)
    tmp = tissuecafiltre.m1952.capi1.d20210721AM.B6_20p.(f{i});
    bsl = mean(tmp(tmp(:, 1) < 10, 2));
    tmp(:, 2) = (tmp(:, 2) - bsl)/max(tmp(:, 2) - bsl);
    tissueCa1952_MSG(:, i) = interp1(tmp(:, 1), tmp(:, 2), 1:0.05:29, 'pchip');
end
means.ca_tissues(:, 4) = mean(tissueCa1952_MSG, 2);


tissueCa1952_G6 = [];
f = fieldnames(tissuecafiltre.m1952.capi2.d20210721PM.B6_20p);
for i =1:length(f)
    tmp = tissuecafiltre.m1952.capi2.d20210721PM.B6_20p.(f{i});
    bsl = mean(tmp(tmp(:, 1) < 10, 2));
    tmp(:, 2) = (tmp(:, 2) - bsl)/max(tmp(:, 2) - bsl);
    tissueCa1952_G6(:, i) = interp1(tmp(:, 1), tmp(:, 2), 1:0.05:29, 'pchip');
end
means.ca_tissues(:, 5) = mean(tissueCa1952_G6, 2);

% Flows
flow1951_G1 = zFlow.m1951.capi1.d20210506AM.B6_10p; % 4 FLOW
means.flows(:, 1) = mean([flow1951_G1.d11(:, 2) flow1951_G1.d12(:, 2) flow1951_G1.d13(:, 2) flow1951_G1.d16(:, 2)], 2, 'omitnan');
flow1951_G3 = zFlow.m1951.capi2.d20210520.B6_10p; % 3 FLOW
means.flows(:, 2) = mean([flow1951_G3.d13(:, 2) flow1951_G3.d14(:, 2) flow1951_G3.d18(:, 2)], 2, 'omitnan');
means.flows(:, 3) = data.M1951_220721_PM_Flow.ZScore.Glom7_Vessel.Flow(end-1, :); % 3 FLOW

flow1952_G4 = zFlow.m1952.capi4.d20210702AM.B6_20p; % 6 FLOW
means.flows(:, 4) = mean([flow1952_G4.d09(:, 2) flow1952_G4.d10(:, 2) flow1952_G4.d11(:, 2) flow1952_G4.d12(:, 2) flow1952_G4.d14(:, 2) flow1952_G4.d15(:, 2)], 2, 'omitnan');
flow1952_G5 = zFlow.m1952.capi5.d20210702PM.B6_20p; % 8 FLOW
means.flows(:, 5) = mean([flow1952_G5.d07(:, 2) flow1952_G5.d08(:, 2) flow1952_G5.d09(:, 2) flow1952_G5.d10(:, 2) flow1952_G5.d11(:, 2) flow1952_G5.d12(:, 2) flow1952_G5.d13(:, 2) flow1952_G5.d14(:, 2)], 2, 'omitnan');

flow2115_G1 = zFlow.m2115.capi1.d20210713PM.B6_20p; % 7 FLOW
means.flows(:, 6) = mean([flow2115_G1.d18(:, 2) flow2115_G1.d19(:, 2) flow2115_G1.d21(:, 2) flow2115_G1.d23(:, 2) flow2115_G1.d24(:, 2) flow2115_G1.d25(:, 2) flow2115_G1.d26(:, 2)], 2, 'omitnan');

t30 = 1:0.1:29; 

% Delays
Delays = struct();
Delays.vessel(1) = getThreshold_fitSig(t30, means.vessels(:, 1), [5 10], [5 14], 0.2, false);
Delays.vessel(2) = getThreshold_fitSig(t30, means.vessels(:, 2), [5 10], [8 14.5], 0.2, false);
Delays.vessel(3) = getThreshold_fitSig(t30, means.vessels(:, 3), [5 10], [8 14.5], 0.2, false);
Delays.vessel(4) = getThreshold_fitSig(t30, means.vessels(:, 4), [5 10], [8 13.5], 0.2, false);
Delays.vessel(5) = getThreshold_fitSig(t30, means.vessels(:, 5), [5 10], [8 14.5], 0.2, false);
Delays.vessel(6) = getThreshold_fitSig(t30, means.vessels(:, 6), [5 10], [8 14.5], 0.2, false);
Delays.vessel(7) = getThreshold_fitSig(t30, means.vessels(:, 8), [5 10], [8 14], 0.2, false);
Delays.vessel(8) = getThreshold_fitSig(t30, means.vessels(:, 9), [5 10], [8 14.5], 0.2, false);
Delays.vessel(9) = getThreshold_fitSig(t30, means.vessels(:, 10), [5 10], [8 14], 0.2, false);
Delays.vessel(10) = getThreshold_fitSig(t30, means.vessels(:, 13), [5 10], [8 13.8], 0.2, false);
Delays.vessel(11) = getThreshold_fitSig(t30, means.vessels(:, 14), [5 10], [8 13.8], 0.2, false);
Delays.vessel(12) = getThreshold_fitSig(t30, msg.vessels(:, 1), [5 10], [8 15], 0.2, false);
Delays.vessel = Delays.vessel -10;

Delays.tissue(1) = getThreshold_fitSig(t30, means.tissues(:, 1), [5 10], [5 14], 0.2, false); % D2
Delays.tissue(2) = getThreshold_fitSig(t30, means.tissues(:, 2), [5 10], [5 14], 0.2, false); % D3
Delays.tissue(3) = getThreshold_fitSig(t30, means.tissues(:, 3), [5 10], [5 15], 0.2, false); % //
Delays.tissue(4) = getThreshold_fitSig(t30, means.tissues(:, 6), [5 10], [5 14], 0.2, false); % //
Delays.tissue(5) = getThreshold_fitSig(t30, means.tissues(:, 7), [5 10], [5 13.5], 0.2, false); % D10
Delays.tissue(6) = getThreshold_fitSig(t30, msg.tissues(:, 1), [5 10], [5 14], 0.2, false); % D12
Delays.tissue = Delays.tissue - 10;

Delays.flow(1) = getThreshold_fitSig(t30, means.flows(:, 1), [5 10], [5 14], 0.2, false); % V1
Delays.flow(2) = getThreshold_fitSig(t30, means.flows(:, 2), [5 10], [5 14], 0.2, false); % V3
Delays.flow(3) = getThreshold_fitSig(t30, means.flows(:, 3), [5 10], [5 14], 0.2, false); % V6
Delays.flow(4) = getThreshold_fitSig(t30, means.flows(:, 4), [5 10], [5 14], 0.2, false); % V10
Delays.flow(5) = getThreshold_fitSig(t30, means.flows(:, 5), [5 10], [5 14], 0.2, false); % V11
Delays.flow(6) = getThreshold_fitSig(t30, means.flows(:, 6), [5 10], [5 14], 0.2, false); % V14
Delays.flow = Delays.flow - 10;

% Figure
figure; 
subplot(311); hold on; set(gca,'defaultAxesColorOrder',[[1 1 1]; [1 1 1]]);
yL = [-2 4];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');
plot(t30, mean([msg.vessels(:, 2) means.vessels(:, [2 3 12])], 2, 'omitnan'), 'LineWidth', 3, 'DisplayName', 'Vessel');
plot(t30, mean([msg.tissues(:, 2) means.tissues(:, [1 2 5])], 2, 'omitnan'), 'LineWidth', 3, 'DisplayName', 'Tissue');
ylabel('Po_{2} (SD)'); ylim(yL);
xlabel('Time (s)'); xlim(xL); legend();
set(gca, 'FontSize', 13); 


subplot(312); hold on;
yL = [-0.5 1.2];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');
m = mean([means.ca_vessels], 2, 'omitnan');
plot(1:0.05:29, m/(max(m)), 'LineWidth', 2,'DisplayName', ...
    'Ca^{2+} Vessel', 'LineWidth', 2, 'Color', [0.20 0.70 0.32]);
m = mean([means.ca_tissues], 2, 'omitnan');
plot(1:0.05:29, m/(max(m)), 'LineWidth', 2, 'DisplayName', ...
    'Ca^{2+} Tissue', 'Color', [0.06 0.38 0.15]);
ylim(yL); xlim([5 20]);
xlabel('Time (s)'); ylabel('Normalized Ca^{2+}'); legend();
set(gca, 'FontSize', 13); 

subplot(313); hold on; yL = [-1 4];
patch('XData', [10 12 12 10], 'YData', [yL(1) yL(1) yL(2) yL(2)], ...
    'FaceColor', colorStim, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'DisplayName', 'Stimulus');
plot(t30, mean([msg.vessels means.vessels], 2, 'omitnan'), ...
    'LineWidth', 3, 'DisplayName', 'Vessel');
plot(t30, mean([msg.tissues means.tissues], 2, 'omitnan'), ...
    'LineWidth', 3, 'DisplayName', 'Tissue');
xlabel('Time (s)'); ylabel('Po_{2} (SD)'); xlim(xL); ylim(yL);
set(gca, 'FontSize', 13); 

suptitle('FIGURE 6B');

figure; hold on;
b = bar([mean(Delays.vessel) mean(Delays.tissue)], 'FaceColor', 'flat');
b.CData(1, :) = colors{1}; b.CData(2, :) = colors{2};
scatter(repmat(1, 1, length(Delays.vessel)), Delays.vessel, 40, ...
    'k', 'filled');
scatter(repmat(2, 1, length(Delays.tissue)), Delays.tissue, 40, ...
    'k', 'filled');
line([repmat(1, 1, length(Delays.vessel([2 3 10 12]))); ...
    repmat(2, 1, length(Delays.vessel([2 3 10 12])))], ...
    [Delays.vessel([2 3 10 12]); Delays.tissue([1 2 5 6])], ...
    'LineWidth', 1, 'Color', 'k');
xticklabels({'Vessel'; 'Tissue'}); xticks([1:2]); xtickangle(45);
xlim([0.5 2.5]); ylim([0 3.5]); 
ylabel('Delay Ca^{2+} - Po_{2} (s)');
set(gca, 'FontSize', 13);

[h, p] = ttest(Delays.vessel([2 3 10 12]), Delays.tissue([1 2 5 6]), 'tail', 'left');
disp('Paired t-test - One-Sided for PO2 onset in vessel < in tissue');
disp(['Significant difference, p = ' ...
    num2str(p)])
[h, p] = ttest2(Delays.vessel, Delays.tissue, 'tail', 'left');
disp('Two-sampled t-test - One-Sided for PO2 onset in vessel < in tissue');
disp(['Significant difference, p = ' ...
    num2str(p)])
end

