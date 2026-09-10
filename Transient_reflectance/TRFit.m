function fitData = TRFit()
%% TRFit.m fits the transient reflectivity data and produces values
%% and uncertainties of [tau_f (ps), A_f+A_s (arb. unit), f_OSC (GHz)]

dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Data');
dataFile = fullfile(dataDir, 'TRData.mat');
data = load(dataFile, 'temperatures', 'time', 'signals');
nTemperatures = numel(data.temperatures);
fitData.temperatures = data.temperatures(:);
fitData.traces = cell(nTemperatures, 1);
fitData.values = nan(nTemperatures, 3);
fitData.uncertainties = nan(nTemperatures, 3);
fitData.perScan = cell(nTemperatures, 1);
for i = 1:nTemperatures
    [fitData.traces{i}, fitData.values(i, :), ...
        fitData.uncertainties(i, :), fitData.perScan{i}] = ...
        fitTemperature(data.time(:), data.signals(:, :, i) / -1);
end
end

function [trace, values, uncertainties, perScan] = fitTemperature(time, signals)
nScans = size(signals, 2)-1;
for column = 1:size(signals, 2)
    background = isfinite(time) & isfinite(signals(:, column)) & time >= 10 & time <= 70;
    signals(:, column) = signals(:, column)-mean(signals(background, column));
end
scans = signals(:, 2:end);
valid = isfinite(time) & time > 10 & all(isfinite(scans), 2);
fitTime = time(valid);
fitScans = scans(valid, :);
du = median(diff(fitTime));
timeSpan = fitTime(end)-fitTime(1);
periodStart = min(max(0.078*timeSpan, 4*du), 0.5*timeSpan);
nSearch = max(3, round(0.25*numel(fitTime)));
scanT0Start = nan(nScans, 1);
for scanIndex = 1:nScans
    [~, jumpIndex] = max(abs(diff(fitScans(1:nSearch, scanIndex))));
    scanT0Start(scanIndex) = fitTime(jumpIndex+1);
end

indices = reshape(1:4*nScans, nScans, 4);
lowerRates = [max(0.5, du); 200; max(du, 300); max([4*du, 70, 1e-4])];
upperRates = [30; 2500; 2500; min(0.5*timeSpan, 130)];
sharedInitial = min(max([1000; 900; periodStart; 1.05], ...
    [lowerRates(2:4); max(0.5*du, 0.5)]), [upperRates(2:4); 2]);
[tauFStart, t0Initial] = estimateEarlyStarts(fitTime, fitScans, scanT0Start, sharedInitial);
lower = [repelem(lowerRates, nScans); t0Initial-1e-4; sharedInitial(4)-1e-4];
upper = [repelem(upperRates, nScans); t0Initial+1e-4; sharedInitial(4)+1e-4];
initial = [tauFStart; repelem(sharedInitial(1:3), nScans); t0Initial; sharedInitial(4)];
initial = min(max(initial, lower), upper);
options = optimoptions('lsqnonlin', 'Display', 'off', ...
    'MaxIterations', 700, 'MaxFunctionEvaluations', 16000, ...
    'FunctionTolerance', 1e-9, 'StepTolerance', 1e-9);
objective = @(p) modelResidual(p, fitTime, fitScans);
bestNorm = Inf;
for startIndex = 1:4
    if startIndex == 1
        start = initial;
    else
        start = randomStart(lower, upper, nScans);
    end
    [parameters, resnorm, residual, exitflag, ~, ~, jacobian] = ...
        lsqnonlin(objective, start, lower, upper, options);
    if exitflag > 0 && isfinite(resnorm) && resnorm < bestNorm
        nonlinear = parameters;
        bestNorm = resnorm;
        bestResidual = residual;
        bestJacobian = jacobian;
    end
end
interval = nlparci(nonlinear, bestResidual, 'jacobian', bestJacobian);
ci95 = 0.5*(interval(:, 2)-interval(:, 1));
ci95(~isfinite(ci95)) = NaN;
[~, components, initialResponse, initialCI95] = ...
    modelResidual(nonlinear, fitTime, fitScans);
perScan = [nonlinear(indices(:, 1)), initialResponse, nonlinear(indices(:, 4))];
halfWidths = [ci95(indices(:, 1)), initialCI95, ci95(indices(:, 4))];
values = nan(1, 3);
uncertainties = nan(1, 3);
for quantity = 1:3
    [values(quantity), uncertainties(quantity)] = ...
        randomEffectsMean(perScan(:, quantity), halfWidths(:, quantity));
end
period = values(3);
values(3) = 1000/period;
uncertainties(3) = 1000*uncertainties(3)/period^2;
perScan(:, 3) = 1000./perScan(:, 3);
scanComponents = nan(numel(time), 3, nScans);
for scanIndex = 1:nScans
    scanComponents(:, :, scanIndex) = interp1(fitTime, ...
        components(:, :, scanIndex), time, 'linear', NaN);
end
scanTimeZero = repmat(nonlinear(end-1) + ...
    nonlinear(end)*erfinv(2*0.01-1), 1, nScans);
trace = [time-mean(scanTimeZero, 'omitnan'), signals(:, 1), ...
    sum(mean(scanComponents, 3, 'omitnan'), 2)];
end

function [tauStart, t0Initial] = estimateEarlyStarts(time, scans, t0Start, sharedInitial)
nScans = size(scans, 2);
tauStart = repmat(7, nScans, 1);
t0Refined = t0Start;
du = median(diff(time));
tauFloor = max(0.25*du, 0.1);
tauGrid = unique([linspace(tauFloor, 20, 50), linspace(20, 30, 31)]);
tauGrid = tauGrid(tauGrid > 0.5 & tauGrid <= 30);
for scanIndex = 1:nScans
    relativeTime = time-t0Start(scanIndex);
    mask = relativeTime >= -5 & relativeTime <= 9;
    t0Grid = linspace(t0Start(scanIndex)-0.75, t0Start(scanIndex)+0.75, 13);
    bestSSE = Inf;
    for t0Candidate = t0Grid
        for tauCandidate = tauGrid
            design = carrierDesignMatrix(time(mask), ...
                [tauCandidate; sharedInitial(1:3); t0Candidate; sharedInitial(4)]);
            linear = lsqminnorm(design, scans(mask, scanIndex));
            residual = design*linear-scans(mask, scanIndex);
            sse = sum(residual.^2);
            if isfinite(sse) && sse < bestSSE
                bestSSE = sse;
                tauStart(scanIndex) = tauCandidate;
                t0Refined(scanIndex) = t0Candidate;
            end
        end
    end
end
t0Initial = median(t0Refined, 'omitnan');
end

function [residual, components, initialResponse, initialCI95] = ...
        modelResidual(nonlinear, time, scans)
nScans = size(scans, 2);
residualMatrix = nan(size(scans));
if nargout > 1
    components = nan(numel(time), 3, nScans);
    initialResponse = nan(nScans, 1);
    initialCI95 = nan(nScans, 1);
end
for scanIndex = 1:nScans
    p = nonlinear([scanIndex+(0:3)*nScans, 4*nScans+(1:2)]);
    formula = carrierDesignMatrix(time, p);
    linear = lsqminnorm(formula, scans(:, scanIndex));
    residualMatrix(:, scanIndex) = formula*linear-scans(:, scanIndex);
    if nargout > 1
        components(:, 1, scanIndex) = formula(:, 1)*linear(1);
        components(:, 2, scanIndex) = formula(:, 2)*linear(2);
        components(:, 3, scanIndex) = formula(:, 3)*linear(3) + formula(:, 4)*linear(4);
        designInverse = pinv(formula'*formula);
        degreesOfFreedom = max(numel(time)-4, 1);
        residualMSE = sum(residualMatrix(:, scanIndex).^2)/degreesOfFreedom;
        linearCovariance = residualMSE*designInverse;
        linearCovariance = (linearCovariance+linearCovariance')/2;
        initialResponse(scanIndex) = linear(1) + linear(2);
        gradient = [1, 1, 0, 0];
        initialCI95(scanIndex) = tinv(0.975, degreesOfFreedom)*sqrt(max( ...
            gradient*linearCovariance*gradient', 0));
    end
end
residual = residualMatrix(:);
end

function design = carrierDesignMatrix(time, p)
tauF = max(p(1), eps);
shiftedTime = time-p(5);
turnOn = 0.5*(1+erf(shiftedTime/p(6)));
fastExponent = max(min(-shiftedTime/tauF, 700), -700);
fastComponent = turnOn.*exp(fastExponent);
slowComponent = turnOn.*exp(-shiftedTime/p(2));
coherentEnvelope = turnOn.*exp(-shiftedTime/p(3));
phase = 2*pi*shiftedTime/p(4);
design = [fastComponent, slowComponent, coherentEnvelope.*cos(phase), ...
    coherentEnvelope.*sin(phase)];
end

function start = randomStart(lower, upper, nScans)
start = lower + rand(size(lower)).*(upper-lower);
for parameterIndex = 1:4*nScans
    start(parameterIndex) = exp(log(lower(parameterIndex)) + rand* ...
        log(upper(parameterIndex)/lower(parameterIndex)));
end
end

function [meanValue, combinedUncertainty] = randomEffectsMean(values, fitCI95HalfWidths)
valid = isfinite(values) & isfinite(fitCI95HalfWidths) & fitCI95HalfWidths >= 0;
values = values(valid);
fitCI95HalfWidths = fitCI95HalfWidths(valid);
n = numel(values);
meanValue = mean(values);
fitStandardErrors = fitCI95HalfWidths/1.95996398454005;
withinVariance = mean(fitStandardErrors.^2);
observedVariance = var(values, 0);
betweenVariance = max(0, observedVariance-withinVariance);
typeA = sqrt(betweenVariance/n);
typeB = sqrt(sum(fitCI95HalfWidths.^2))/n;
combinedUncertainty = hypot(typeA, typeB);
end
