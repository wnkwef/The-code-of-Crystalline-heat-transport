clear; close all; clc

% Load photoluminescence (PL, here Em) and photoluminescence excitation
% (PLE, here Ex) Data.
Data = EmExData();

f = figure(Position=[200, 150, 800, 500]);
tl = tiledlayout(f, 1, 2, TileSpacing="compact", Padding="tight");
%% Emission spectra
ax1 = nexttile(tl, 1); hold(ax1, 'on');

for idx = 1:numel(Data.NumT)
    x = Data.emX{idx};
    y = Data.emY{idx};
    x_fit = Data.emFitX;
    y_fit = Data.emFitY{idx};
    offset = 1.1*(idx-1);
    c = Data.grayColors(idx,:);

    % Downsampling factor
    Fac = 3;
    plot(ax1, downsample(x, Fac), downsample(y + offset, Fac), '-', ...
         'LineWidth',1.5, 'Color',c);
    plot(ax1, x_fit, y_fit + offset, '--', 'LineWidth',1.5, 'Color',c);
    fill(ax1, [x_fit; flip(x_fit)], ...
         [y_fit+offset; offset*ones(size(y_fit))], c, ...
         'FaceAlpha',0.3, 'EdgeColor','none');

    plot(ax1, Data.peakX(idx), Data.peakY(idx), 'p', ...
         'MarkerSize',10, 'LineWidth',1.5, 'Color',[0,0.5,1]);
    txt = sprintf('%.0f K', Data.NumT(idx));
    text(ax1, 2.7, offset+0.5, txt, Color=[0,0.5,1], FontSize=10, ...
         HorizontalAlignment='left', VerticalAlignment='middle');
end

text(ax1, 0.05, 0.98, 'a', 'Units','normalized', ...
     'FontSize',16, 'FontWeight','bold');
xlabel(ax1, 'E (eV)', 'FontSize',12);

yq = linspace(min(Data.peakY), max(Data.peakY), 100);
xq = pchip(Data.peakY, Data.peakX, yq);
plot(ax1, xq, yq, 'r--', LineWidth=2);
xLower = pchip(Data.peakY, Data.muLower, yq);
plot(ax1, xLower, yq, 'b:', LineWidth=2);
xUpper = pchip(Data.peakY, Data.muUpper, yq);
plot(ax1, xUpper, yq, 'b:', LineWidth=2);
xlim(ax1, [1.4 3]);
ylim(ax1, [0, 1.1*numel(Data.NumT)]);
hold(ax1, 'off');

%% Excitation spectra
ax2 = nexttile(tl, 2); hold(ax2, 'on');
offset = zeros(numel(Data.NumT), 1);

for idx = 1:numel(Data.NumT)
    x_ex = Data.exX{idx};
    y_ex = Data.exY{idx};
    offset(idx) = 1.1*(idx-1);
    c = Data.grayColors(idx,:);

    plot(ax2, Data.ExIncpt(idx), offset(idx), 'v', ...
         'MarkerSize',8, 'LineWidth',1.5, 'Color',[0,0.5,1]);
    Fac = 2;
    plot(ax2, downsample(x_ex, Fac), downsample(y_ex + offset(idx), Fac), ...
         '-', 'LineWidth',1.5, 'Color',c);

    txt = sprintf('%.0f K', Data.NumT(idx));
    loc = offset(idx) + 0.5;
    switch idx
        case 17
            loc = loc + 0.2;
        case 18
            loc = loc + 0.3;
        case 19
            loc = loc + 0.45;
        case 20
            loc = loc - 0.15;
    end
    text(ax2, 2.7, loc, txt, Color=[0,0.5,1], FontSize=10, ...
         HorizontalAlignment='left', VerticalAlignment='middle');
end

text(ax2, 0.05, 0.98, 'b', 'Units','normalized', ...
     'FontSize',16, 'FontWeight','bold');
xlabel(ax2, 'E (eV)', 'FontSize',12);

yq = linspace(0, offset(end), 100);
xq = pchip(offset, Data.ExIncpt, yq);
plot(ax2, xq, yq, 'r--', LineWidth=2);
xLower = pchip(offset, Data.ExIncpt - Data.ExErrLow, yq);
plot(ax2, xLower, yq, 'b:', LineWidth=2);
xUpper = pchip(offset, Data.ExIncpt + Data.ExErrHigh, yq);
plot(ax2, xUpper, yq, 'b:', LineWidth=2);
xlim(ax2, [2.5 4.0]);
ylim(ax2, [0, 1.1*numel(Data.NumT)]);
hold(ax2, 'off');
