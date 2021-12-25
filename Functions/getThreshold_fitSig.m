function [threshold, output] = getThreshold_fitSig(time, y, baselineTime, fitTime, propThresh, toPlot)
% GETTHRESHOLD_FITSIG
%
% function [threshold, output] = getThreshold_fitSig(time, y, baselineTime, fitTime, propThresh, toPlot)
%
%   Author : Ali-Kemal Aydin, PhD student
%   Date : April 17th, 2020
%   Mail: ali-kemal.aydin@inserm.fr
%   Affiliation : INSERM U1128, Paris & U968, Institut de la Vision, Paris
%   License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
%       See LICENSE.txt or <a href="matlab:web('https://creativecommons.org/licenses/by/4.0/')">here</a> 
%       for a human-readable version.
%
%   DESCRIPTION : Out the propThresh threshold on the sigmoidal fit 
%   on (time,y).
%   Will fit a sigmoid function on the (time, y) data and output the
%   threshold as the time when the function gets above propThresh of 
%   the maximum, and an output structure containing details about the fit 
%   and the error made on the threshold. toPlot is a boolean to set if the 
%   fitting plots are to be plotted out.
%
%__________________________________________________________________________
% PARAMETERS:
%	time ([]): Time vector.
%
% 	y ([]): Data vector.
%
%	baselineTime ([]): time extremities of the baseline [Beginning
% 	End].
%
%	fitTime ([]): time extremities of the period to fit [Beginning
%   End].
%
%   propThresh (float): proportion of the threshold to reach to
%   determine activation (< 1).
%
%   toPlot (bool) : if true, plot the output.
%
%__________________________________________________________________________
% RETURN:
%   threshold (float): time point of the activation, as determined by
%   the fit and propThresh parameter.
%
%   output (struct): various properties of the fit (error, fit, gof,
%   ...).
%
%__________________________________________________________________________

%% Parameters
dt = mean(diff(time));
bslIdx = [find(time>baselineTime(1), 1) find(time>baselineTime(2), 1)];
fitIdx = [find(time>fitTime(1), 1) find(time>fitTime(2), 1)];

%% Computing starting values for the optimization
% Param 'a' is the baseline of the sigmoid, starting is the mean of the
% baseline
bslValue = mean(y(bslIdx(1):bslIdx(2)));
% Param 'b' is the maximum of the sigmoid, starting is the maximum of the
% trace in the fitting interval
maxValue = max(y(fitIdx(1):fitIdx(2)));
% Param 'c' is the time at half max, starting is the half between end of
% the baseline and maximum of the trace
maxTime = time(find(y' == maxValue, 1));
halfTime = maxTime - (maxTime - baselineTime(2))/2;
% Param 'd' is the rate of the sigmoid, starting is the slope between max
% and end of the baseline.
rate = (maxValue - bslValue)/((maxTime - baselineTime(2)));

%% Fitting


ft = fittype( 'a+b/(1+exp((c-x)/d))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf -Inf];
opts.MaxIter = 400000;
opts.StartPoint = [bslValue maxValue halfTime rate];
opts.Upper = [Inf Inf Inf Inf];
opts.Exclude = excludedata(time, y', 'domain', fitTime);
    
[fitresult, gof, outFit] = fit(time', y, ft, opts);

gof.nmrse = gof.rmse/mean(y(fitIdx(1):fitIdx(2)));

if toPlot
    % Plot fit with data.
    figure('Name', 'Sigmoïd fit' );
    h = plot(fitresult, time, y);
    legend( h, 'y vs. time', 'Sigmoïd fit', 'Location', 'NorthEast' );
    % Label axes
    xlabel time
    ylabel y
    grid on
    
    figure('Name', 'Sigmoïd fit' );
    h = plot(fitresult, time(fitIdx(1):fitIdx(2)), y(fitIdx(1):fitIdx(2)));
    legend(h, 'y vs. time', 'Sigmoïd fit', 'Location', 'NorthEast' );
    % Label axes
    xlabel time
    ylabel y
    grid on
end

%% Computing Threshold and error
% Threshold is computed as the solution to the equation yThresh = f(x) when
% solving for x. f(x) is the fitted sigmoid equation, yThresh is the 
% thresholded value.
% After some unpredicted bugs, I defined yThresh as 
% limit(sig, x->0+) + propThresh * amplitude_of_the_sigmoid
% to avoid problems of negativity of the sigmoid. 
% Error on the guessed threshold time is computed using the confidence
% interval of the fitted parameters, which I used as approximations for
% STD.
params = coeffvalues(fitresult);
paramsConfInt = confint(fitresult);
deltaParams = (paramsConfInt(2, :) - paramsConfInt(1, :))/2;

a = params(1); b = params(2); c = params(3); d = params(4);

sig = @(x)params(1)+params(2)/(1+exp((params(3)-x)/params(4)));
yThresh = sig(1E-20) + ...
    propThresh * (sig(1E20) - sig(1E-20));
threshold  = c - d * log(b/(yThresh-a)-1);

% Partial derivatives for the error computation
parDev_a = (b*d)/((b/(a - yThresh) + 1)*(a - yThresh)^2); 
parDev_b = -d/((b/(a - yThresh) + 1)*(a - yThresh));
parDev_c = 1;
parDev_d = -log(- b/(a - yThresh) - 1);

err = sqrt(...
    (parDev_a * deltaParams(1))^2 + ...
    (parDev_b * deltaParams(2))^2 + ...
    (parDev_c * deltaParams(3))^2 + ...
    (parDev_d * deltaParams(4))^2);

%% Outputs
threshold = double(threshold);
output.errorOnThresh = double(err);
output.fitresult = fitresult;
output.gof = gof;
output.outFit = outFit;

end

