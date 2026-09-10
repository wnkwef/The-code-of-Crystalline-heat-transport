scriptDir = fileparts(mfilename('fullpath')); addpath(scriptDir);
dataDir = fullfile(fileparts(scriptDir),'Data');

% Load Raman Data
data = load(fullfile(dataDir,'Raman_Parallel.mat'),'data_parallel');
% Load an sorted array of Temperatures.
temps = load(fullfile(dataDir,'Temperature_numeric.mat'),'NumT','Temperatures');
NumT = temps.NumT(:);
% Load the positions of specified Raman peaks 
peaks = load(fullfile(dataDir,'pks.mat'),'Locations_parallel');
% Invoke LorentzianFitWL.m function to fit data
[centers, widths] = LorentzianFitWL(data.data_parallel, ...
    NumT, temps.Temperatures, peaks.Locations_parallel);

%% Plot details
frequencyRanges = [86,103; 116,134; 140,159; 234,246; 268,288];
markers = ["o","s","^","d","p"];
colors = [0,114,178; 230,85,13; 0,158,115; 204,121,167; 92,70,156]/255;
% Number of frequency ranges
nRanges = size(frequencyRanges,1);
lifetimeData = nan(nRanges,numel(NumT)); peakCounts = zeros(size(lifetimeData));
for index = 1:numel(NumT)
    selected = centers{index} > frequencyRanges(:,1)' & centers{index} < frequencyRanges(:,2)';
    lifetimes = 1e12./(2*pi*299792458*100*widths{index});
    peakCounts(:,index) = sum(selected,1)';
    lifetimeData(:,index) = (lifetimes'*selected)'./peakCounts(:,index);
end
% The average lifetime
meanLifetime = mean(lifetimeData,1);
% Wigner limit τ
tau = 3*14/(2*pi*8.6);
f = figure('Position',[600,160,660,365],'Color','w');
ax = axes(f); hold(ax,'on'); box(ax,'on');
handles = gobjects(1,nRanges);
for band = 1:nRanges
    % Filled markers indicate the same peak COUNT as the first input temperature (8K).
    % Hollow markers indicate less peak COUNT than the first input temperature (8K).
    filled = peakCounts(band,:) == peakCounts(band,1);
    plot(ax,NumT,lifetimeData(band,:),markers(band),'Color',colors(band,:),'MarkerFaceColor','none');
    % Plot Filled markers
    handles(band) = plot(ax,NumT(filled),lifetimeData(band,filled),markers(band), ...
        'Color',colors(band,:),'MarkerFaceColor',colors(band,:));
end
hMean = plot(ax,NumT,meanLifetime,'-o','Color',[.15,.15,.15], ...
    'LineWidth',1.5,'MarkerSize',5,'MarkerFaceColor','w');
hWigner = yline(ax,tau,'--','Color',[.861,.071,.088],'LineWidth',2);
xlabel(ax,'Temperature (K)'); ylabel(ax,'Lifetime (ps)');
yLimits = [min(.15,.95*min([lifetimeData(:);tau])), max(1.55,1.05*max([lifetimeData(:);tau]))];
set(ax,'FontName','Arial','FontSize',12,'LineWidth',1,'XTick',0:25:300, ...
    'YTick',0.2:0.2:yLimits(2),'TickLength',[.005,.005],'XLim',[0,300],'YLim',yLimits);
ax.Toolbar.Visible = 'off';
labels = compose('%g–%g cm^{-1}',frequencyRanges(:,1),frequencyRanges(:,2));
lgd = legend(ax,[handles,hMean,hWigner],[labels; "Average"; "Wigner Limit"], ...
    'Location','northeast','FontSize',10,'FontName','Arial');
lgd.ItemTokenSize = [10,10];
