function Data = EmExData()
%% EmExData prepares data fitting and Jacobian conversion

%% Load Data
ScriptDir = fileparts(mfilename('fullpath'));
DataDir = fullfile(fileparts(ScriptDir), 'Data');

% demo
load(fullfile(DataDir,"SSData.mat"))

%% Load T-dependent Emission data
% D = readtable(fullfile(DataDir, 'TempEm.txt'), ...
%     DataLine=3, VariableNamingRule='preserve');
% 
% NumT = regexp(D.Properties.VariableNames, '\d+', 'match', 'once');
% NumT = str2double(NumT);
% NumT = NumT(~isnan(NumT));
% nTemperatures = numel(NumT);
% 
% % Jacobian conversion of Em data
% D.Energy = 1240 ./ D.Wavelength;
% for idx = 1:nTemperatures
%     D{:,idx+1} = D{:,idx+1} .* 1240 ./ D.Energy.^2;
% end
% 
% % Load T-dependent Excitation data
% SheetFile = fullfile(DataDir, 'GeneralSheet.xlsx');
% nameRow = readcell(SheetFile, 'Sheet',1, 'Range','1:1');
% typeRow = readcell(SheetFile, 'Sheet',1, 'Range','2:2');
% raw = readmatrix(SheetFile, 'Sheet',1, 'NumHeaderLines',2);
% exCols = find(strcmpi(typeRow, 'Excitation Scan'));
% exNames = nameRow(exCols);
% exJD = arrayfun(@(c) rmmissing(raw(:,[c-1 c])), exCols, uni=0);
% exJD = cellfun(@(xy) [1240./xy(:,1), 1240.*xy(:,2)./(1240./xy(:,1)).^2], exJD, uni=0);
% 
% % Remove several indices not used for further analysis due to repetition or
% % extremely low quality.
% removeIndices = [26,24,21,22,6,4,1];
% for idx = removeIndices
%     exJD(idx) = [];
%     exNames(idx) = [];
% end
% 
% % Jacobian conversion of Ex data
% exJD = cellfun(@(xy) [xy(:,1), (xy(:,1).*xy(:,2)).^2], exJD, uni=0);
% exTemperatures = cellfun(@(name) str2double(regexp(name, '(\d+)K', 'tokens', 'once')), exNames);
% [~, idx_sort] = sort(exTemperatures);
% exJD = exJD(idx_sort);
% 
% % Pre-defined bandgaps by Tauc plot
% ex_data = readtable(fullfile(DataDir, 'ExPeaks_IntcpResults.xlsx'), ...
%     VariableNamingRule='preserve');
% ex_data.Temp = cellfun(@(x) str2double(regexp(x, ...
%     '(\d+)(?=K)', 'match', 'once')), ex_data.Sample);

Data.NumT = NumT;
grayLevels = linspace(0.75, 0.25, nTemperatures).';
Data.grayColors = [grayLevels grayLevels grayLevels];
Data.ExIncpt = arrayfun(@(T) ex_data.xIntercepts(ex_data.Temp == T), NumT);
Data.ExErrHigh = arrayfun(@(T) ex_data.ErrorHigh(ex_data.Temp == T), NumT);
Data.ExErrLow = arrayfun(@(T) ex_data.ErrorLow(ex_data.Temp == T), NumT);

[Data.peakX, Data.peakY, Data.muLower, Data.muUpper, ...
    Data.peakWidth, Data.emRSquared] = deal(zeros(1, nTemperatures));
[Data.emX, Data.emY, Data.emFitY, Data.exX, Data.exY] = ...
    deal(cell(1, nTemperatures));
Data.emFitX = (1.4:1e-2:3)';

% Gaussian fitting of Em data to determine the emission peaks
ft = fittype(@(A1, mu1, sigma1, x) ...
    A1 * exp(-((x - mu1).^2) / (2 * sigma1^2)), ...
    'independent', 'x', 'coefficients', {'A1','mu1','sigma1'});

for idx = 1:nTemperatures
    x = D.Energy;
    y = D{:, idx+1};

    edges = [110 170 200 inf];
    xHiTab = [3.00 3.00 2.30 2.30];
    SPsTab = [0.5 1.95 0.5;
              0.5 1.65 0.5;
              0.5 1.66 0.5;
              0.5 1.65 1.0];
    idx_row = find(NumT(idx) <= edges, 1);
    sel = x > 1.60 & x < xHiTab(idx_row);
    SPs = SPsTab(idx_row,:);
    x = x(sel);
    y = y(sel);
    y = y ./ max(y);

    edges = [1 3 8 9 10 19];
    values = [2.3 2.03 1.87 1.72 1.65 1.70];
    me = values(sum(idx >= edges));
    st = 0.15;
    weights = exp(-((x - me).^2) / (2 * st^2));
    weights = weights / max(weights);

    opts = fitoptions('Method','NonlinearLeastSquares', ...
        'StartPoint',SPs, ...
        'MaxIter',1e15, ...
        'TolFun',1e-15, ...
        'TolX',1e-15, ...
        'Lower',[0, 0, 0], ...
        'Upper',[2, 2.5, 2], ...
        'Weights',weights);
    [fitresult, gof] = fit(x, y, ft, opts);
    ci = confint(fitresult);

    Data.emX{idx} = x;
    Data.emY{idx} = y;
    Data.emFitY{idx} = feval(fitresult, Data.emFitX);
    Data.peakX(idx) = fitresult.mu1;
    Data.peakY(idx) = fitresult.A1 + 1.1*(idx-1);
    Data.muLower(idx) = ci(1,2);
    Data.muUpper(idx) = ci(2,2);
    Data.peakWidth(idx) = fitresult.sigma1;
    Data.emRSquared(idx) = gof.rsquare;
    fprintf('T = %dK  R^2 = %.3f\n', NumT(idx), gof.rsquare);

    x_ex = exJD{idx}(:,1);
    y_ex = exJD{idx}(:,2);
    idx_range = x_ex >= 2.5 & x_ex <= 4.0;
    x_ex = x_ex(idx_range);
    y_ex = y_ex(idx_range);
    Data.exX{idx} = x_ex;
    Data.exY{idx} = y_ex ./ max(y_ex);
end
end
