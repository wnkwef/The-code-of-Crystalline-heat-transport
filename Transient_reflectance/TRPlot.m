clear; close all; clc

% Fit TR data
fitData = TRFit;
temperatures = fitData.temperatures;
nTemp = numel(temperatures);

mainFigure = figure('Tag', 'TRMainFigure');
set(mainFigure, 'Name', 'Transient reflectance', 'NumberTitle', 'off', ...
    'Color', 'w', 'Position', [200, 100, 1000, 720]);
positions = [0.075, 0.54, 0.42, 0.41; 0.55, 0.54, 0.27, 0.41; ...
    0.075, 0.10, 0.245, 0.36; 0.395, 0.10, 0.245, 0.36; ...
    0.715, 0.10, 0.245, 0.36];
ax = gobjects(5, 1);
for panel = 1:5
    ax(panel) = axes(mainFigure, 'Position', positions(panel, :), ...
        'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 0.8, ...
        'Box', 'on', 'TickDir', 'in', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.14);
    hold(ax(panel), 'on');
end

%% Full trace and zoom in plots
zoomWindow = [-3, 15];
fraction = linspace(0, 1, nTemp)';
traceColors = (1-fraction).*[0.08, 0.30, 0.88] + fraction.*[0.90, 0.12, 0.08];
pointSizes = [10, 24];
curvePointCounts = [20000, 1500];
for index = 1:nTemp
    trace = fitData.traces{index};
    offset = (index-1)*0.8;
    valid = isfinite(trace(:, 1)) & isfinite(trace(:, 3));
    [fitTime, uniqueIndex] = unique(trace(valid, 1), 'sorted');
    fitSignal = trace(valid, 3);
    fitSignal = fitSignal(uniqueIndex);
    windows = [trace(1, 1), trace(end, 1); zoomWindow];
    for panel = 1:2
        mask = trace(:, 1) >= windows(panel, 1) & trace(:, 1) <= windows(panel, 2);
        scatter(ax(panel), trace(mask, 1), trace(mask, 2)+offset, pointSizes(panel), ...
            traceColors(index, :), 'filled', 'MarkerEdgeColor', 'none', ...
            'DisplayName', sprintf('%g K', temperatures(index)));
        denseTime = linspace(max(windows(panel, 1), fitTime(1)), ...
            min(windows(panel, 2), fitTime(end)), curvePointCounts(panel));
        plot(ax(panel), denseTime, interp1(fitTime, fitSignal, denseTime, 'pchip')+offset, ...
            'Color', traceColors(index, :), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end
xlim(ax(1), fitData.traces{1}([1, end], 1)');
xlim(ax(2), zoomWindow);
xticks(ax(2), [0, 5, 10, 15]);
for panel = 1:2
    ylim(ax(panel), [-0.5, nTemp*0.8+3]);
    xlabel(ax(panel), 't (ps)');
end
ylabel(ax(1), '\Delta R/R_0');
xline(ax(2), 0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(ax(2), 1.3, 0.25, '$t=0$', 'Interpreter', 'latex', 'FontSize', 11);
temperatureLegend = legend(ax(2), 'Location', 'eastoutside', 'Box', 'off', ...
    'FontName', 'Arial', 'FontSize', 11);
ax(2).Position = positions(2, :);
temperatureLegend.Position = [0.84, 0.54, 0.14, 0.41];

%% Fast lifetime, initial response, and oscillation frequency
colors = [0.10, 0.20, 1.00; 0.85, 0.10, 0.90; 1.00, 0.15, 0.15];
labels = {'\tau_f (ps)', 'A_f + A_s (arb. unit)', 'f_{OSC} (GHz)'};
scanCounts = cellfun(@(values) size(values, 1), fitData.perScan);
scanTemperatures = repelem(temperatures, scanCounts);
scanValues = vertcat(fitData.perScan{:});
T_fit = linspace(min(temperatures), max(temperatures), 600);
for quantity = 1:3
    panel = quantity+2;
    color = colors(quantity, :);
    individual = scatter(ax(panel), scanTemperatures, scanValues(:, quantity), ...
        13, color, 's', 'filled', 'MarkerFaceAlpha', 0.30, 'MarkerEdgeColor', 'none');
    plot(ax(panel), T_fit, interp1(temperatures, ...
        fitData.values(:, quantity), T_fit, 'pchip'), ...
        'Color', color, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    average = errorbar(ax(panel), temperatures, fitData.values(:, quantity), ...
        fitData.uncertainties(:, quantity), 's', 'LineStyle', 'none', ...
        'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 4, ...
        'LineWidth', 1.0, 'CapSize', 4);
    xlabel(ax(panel), 'T (K)');
    ylabel(ax(panel), labels{quantity});
    xlim(ax(panel), [min(temperatures)-7, max(temperatures)+7]);
    xticks(ax(panel), 80:40:240);
    limits = [scanValues(:, quantity); ...
        fitData.values(:, quantity)-fitData.uncertainties(:, quantity); ...
        fitData.values(:, quantity)+fitData.uncertainties(:, quantity)];
    limits = [min(limits, [], 'omitnan'), max(limits, [], 'omitnan')];
    ylim(ax(panel), limits);
    legend(ax(panel), [individual, average], {'Exp.', 'Ave.'}, ...
        'Location', 'best', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 10);
end
set(ax, 'NextPlot', 'replace');
drawnow;
