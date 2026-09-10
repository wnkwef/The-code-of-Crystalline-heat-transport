clear; close all; clc

scriptDir = fileparts(mfilename("fullpath"));
dataDir = fullfile(fileparts(scriptDir),"Data");
modeColors = [15,107,109;198,93,58;110,90,138]/255;
% Calibrate Degrees of different angular dependent input
% calibrationTemperatureK = [8,50,125,290];
% calibrationOffsetDeg = [0,-15,-15,-15];

% demo
calibrationTemperatureK = 8;
calibrationOffsetDeg = 0;

agBgSettings = struct("ratioBounds",[-20,20],"ratioGridPoints",2001, ...
    "optimOptions",optimset(Display="off",TolX=1e-10,MaxIter=1000,MaxFunEvals=3000));

parallelInput = loadRequiredInputs(dataDir,"Raman_Parallel.mat","data_parallel");
crossInput = loadRequiredInputs(dataDir,"Raman_Cross.mat","data_cross");
temperatureInput = loadRequiredInputs(dataDir,"Temperature_numeric.mat",["NumT","Temperatures"]);
data_parallel = parallelInput.data_parallel;
data_cross = crossInput.data_cross;
NumT = temperatureInput.NumT;
Temperatures = temperatureInput.Temperatures;
[peakData,globalResult,thetaByTemperature,appliedAngleOffsetsDeg] = analyze( ...
    scriptDir,data_parallel,data_cross,NumT,Temperatures, ...
    calibrationTemperatureK,calibrationOffsetDeg,agBgSettings);
targetTemperatures = peakData.targetTemperatures(:).';
angles = peakData.angles(:);
intensityParallel = peakData.intensityParallel;
intensityCross = peakData.intensityCross;
peakPositions = peakData.referencePeaks(:).';
peakLabels = string(peakPositions)+" cm^{-1}";
nTemperatures = numel(targetTemperatures);
thetaFine = linspace(0,2*pi,1000).';
fitTable = buildFitTable(targetTemperatures,appliedAngleOffsetsDeg,globalResult,peakPositions);

decompositionFigure = figure(Tag="RamanPolarimetry");
set(decompositionFigure,Position=[250,200,1200,250],Color="w",Name="Raman polarimetry")
decompositionLayout = tiledlayout(decompositionFigure,nTemperatures,6, ...
    TileSpacing="compact",Padding="compact");
dataSets = {intensityParallel,intensityCross};
suffixes = ["Parallel","Cross"];
polarizations = ["pp","cp"];
legendHandles = gobjects(1,numel(peakPositions));
for t = 1:nTemperatures
    for g = 1:2
        suffix = suffixes(g);
        rawData = reshape(dataSets{g}(t,:,:),numel(angles),[]);
        panel = buildPanel(thetaFine,rawData,globalResult,t,g);
        weights = globalResult.("weightAg"+suffix)(t,:);
        curves = {panel.total,panel.ag,panel.bg};
        panelTitles = [sprintf("%g K | %s", ...
            targetTemperatures(t),polarizations(g)),"A_g","B_g"];
        for p = 1:3
            ax = polaraxes(decompositionLayout);
            ax.Layout.Tile = (t-1)*6+(g-1)*3+p;
            hold(ax,"on")
            for m = 1:numel(peakPositions)
                lineHandle = polarplot(ax,thetaFine,curves{p}(:,m), ...
                    Color=modeColors(m,:),LineWidth=1.6,Tag="RamanFit",UserData=[t,g,p,m]);
                if t==1 && g==1 && p==1, legendHandles(m) = lineHandle; end
                if p==1
                    polarplot(ax,thetaByTemperature(:,t),panel.data(:,m),"o", ...
                        MarkerSize=3.5,Color=modeColors(m,:),MarkerFaceColor=modeColors(m,:), ...
                        Tag="RamanData",UserData=[t,g,m]);
                else
                    fractions = [weights(m),1-weights(m)];
                    labels = ["f_A","f_B"];
                    text(ax,1.20,0.03+0.09*(numel(peakPositions)-m+1), ...
                        sprintf("%s=%.2f",labels(p-1),fractions(p-1)), ...
                        Units="normalized",Color=modeColors(m,:),FontName="Arial", ...
                        FontSize=8,HorizontalAlignment="right",Clipping="off");
                end
            end
            radius = 1.1;
            if p==1, radius = 1.1*max([panel.data(:);panel.total(:)],[],"omitnan"); end
            if ~isfinite(radius) || radius<=0, radius = 1; end
            set(ax,RLim=[0,radius],ThetaZeroLocation="right",ThetaDir="counterclockwise", ...
                ThetaTick=0:90:270,RTick=[],ThetaGrid="on",RGrid="off", ...
                GridLineStyle="--",GridAlpha=0.5,LineWidth=1,FontName="Arial",FontSize=7);
            ax.Toolbar.Visible = "off";
            title(ax,panelTitles(p),FontSize=8,FontWeight="normal")
        end
    end
end
peakLegend = legend(legendHandles,peakLabels,Orientation="horizontal",Box="off",FontSize=9);
peakLegend.Layout.Tile = "south";

function [peakData,result,theta,offsets] = analyze(folder,parallel,cross,NumT,Temperatures,calT,calOffset,settings)
    previousFolder = pwd;
    previousPath = path;
    cleanup = onCleanup(@() restoreEnvironment(previousFolder,previousPath));
    cd(folder)
    addpath(folder,"-begin")
    peakData = MainPeaksEgTensorInput(parallel,cross,NumT,Temperatures);
    [~,index] = ismember(peakData.targetTemperatures,calT);
    offsets = calOffset(index);
    theta = deg2rad(peakData.angles(:)+offsets(:).');
    result = FitAgBgModel(theta,peakData.intensityParallel,peakData.intensityCross,settings);
end

function restoreEnvironment(folder,searchPath)
    cd(folder)
    path(searchPath)
end

function panel = buildPanel(theta,data,result,t,g)
    suffixes = ["Parallel","Cross"];
    r = result.ratio;
    c = cos(theta);
    s = sin(theta);
    if g==1
        ag = (c.^2+s.^2.*r).^2;
        bg = (2*c.*s).^2;
    else
        ag = (c.*s.*(r-1)).^2;
        bg = (c.^2-s.^2).^2;
    end
    ag = ag.*result.("coefAg"+suffixes(g))(t,:);
    bg = bg.*result.("coefBg"+suffixes(g))(t,:);
    total = ag+bg;
    scale = max(total,[],"all","omitnan");
    if ~isfinite(scale) || scale<=0, scale = 1; end
    panel = struct("data",data/scale,"ag",ag/scale,"bg",bg/scale,"total",total/scale);
end

function fitTable = buildFitTable(temperatures,offsets,result,peaks)
    n = numel(temperatures);
    fitTable = table(temperatures(:),offsets(:), ...
        VariableNames=["Temperature_K","AppliedAngleOffset_deg"]);
    names = ["b_over_a","CoefAg_pp","CoefBg_pp","WeightAg_pp","WeightBg_pp", ...
        "CoefAg_cp","CoefBg_cp","WeightAg_cp","WeightBg_cp"];
    for m = 1:numel(peaks)
        values = [repmat(result.ratio(m),n,1),result.coefAgParallel(:,m),result.coefBgParallel(:,m), ...
            result.weightAgParallel(:,m),1-result.weightAgParallel(:,m), ...
            result.coefAgCross(:,m),result.coefBgCross(:,m),result.weightAgCross(:,m), ...
            1-result.weightAgCross(:,m)];
        fitTable(:,cellstr(names+"_"+peaks(m))) = array2table(values);
    end
end

function contents = loadRequiredInputs(dataDir,fileName,requiredVariables)
    names = cellstr(requiredVariables);
    contents = load(fullfile(dataDir,fileName),names{:});
end
