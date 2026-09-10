function peakData = MainPeaksEgTensorInput(data_parallel,data_cross,NumT,Temperatures,options)
% Fit Raman spectra in memory; intensities are component heights at x = v.
% Wavenumbers are in cm^-1, temperatures in K, and g is a damping parameter.
arguments
    data_parallel
    data_cross
    NumT
    Temperatures
    options.MaxIter = 5000
    options.MaxFunEvals = 10000
end

% targets = [8,50,125,290];

% demo
targets = 8;
angles = (0:15:360).';
peaks = [271,280,287];
shiftRates = [0.008,0.004,0.015; 0.010,0.003,0.015];
settings = struct('fitRange',[260,300],'referenceCenterStep',4, ...
    'polarCenterStep',0.25,'widthTolerance',0.05,'widthLowerBound',0.2, ...
    'referenceWidthUpperBound',25);
NumT = double(NumT(:));
if isstring(Temperatures)
    Temperatures = cellstr(Temperatures);
end
Temperatures = reshape(Temperatures,1,[]);
[~,indices] = ismember(targets,NumT);
temperatures = Temperatures(indices);
degrees = compose("D%d",angles);
data = {data_parallel,data_cross};
nT = numel(targets); nA = numel(angles); nP = numel(peaks);
locations = zeros(nT,nP,2); widths = locations;
intensities = zeros(nT,nA,nP,2);
spectra = cell(nT,2);
peakIndex = 1:nP;
terms = compose('c%d*g%d.^3./(x.^2.*g%d.^2 + (x.^2 - v%d.^2).^2)', ...
    peakIndex(:),peakIndex(:),peakIndex(:),peakIndex(:));
coefficientNames = cellstr([reshape(["c"+peakIndex;"v"+peakIndex;"g"+peakIndex],1,[]),"b0","b1"]);
formula = sprintf('%s + b0 + b1*(x - %.15g)',strjoin(terms,' + '),mean(settings.fitRange));
scaleNames = cellstr("s"+(1:numel(coefficientNames)));
for index = 1:numel(coefficientNames)
    name = coefficientNames{index};
    formula = regexprep(formula,['\<',name,'\>'],sprintf('(%s*%s)',name,scaleNames{index}));
end
model = fittype(formula,'independent','x','coefficients',coefficientNames,'problem',scaleNames);
solver = fitoptions('Method','NonlinearLeastSquares','MaxIter',options.MaxIter, ...
    'MaxFunEvals',options.MaxFunEvals,'TolFun',1e-10,'TolX',1e-10, ...
    'Display','off','Robust','Bisquare');

for t = 1:nT
    for g = 1:2
        raw = data{g}.(temperatures{t});
        x = raw.wavenumber(:);
        y = zeros(numel(x),nA);
        for a = 1:nA
            y(:,a) = raw.(degrees(a));
        end
        exponent = 6.62607015e-34*299792458*100.*x./(1.380649e-23*targets(t));
        y = y.*(-expm1(-exponent)./x);
        spectra{t,g} = struct('x',x,'y',y);
        scale = max(y(:,1));
        parameters = fitPeaks( ...
            x,y(:,1)./scale,peaks-targets(t)*shiftRates(g,:),[],settings,model,solver);
        locations(t,:,g) = parameters(2,:);
        widths(t,:,g) = parameters(3,:);
    end
end
for t = 1:nT
    for g = 1:2
        spectrum = spectra{t,g};
        for a = 1:nA
            parameters = fitPeaks(spectrum.x,spectrum.y(:,a), ...
                locations(t,:,g),widths(t,:,g),settings,model,solver);
            intensities(t,a,:,g) = parameters(1,:).*parameters(3,:)./parameters(2,:).^2;
        end
    end
end
peakData = struct('targetTemperatures',targets,'temperatures',{temperatures}, ...
    'angles',angles,'referencePeaks',peaks,'locationsParallel',locations(:,:,1), ...
    'widthsParallel',widths(:,:,1),'locationsCross',locations(:,:,2), ...
    'widthsCross',widths(:,:,2),'intensityParallel',intensities(:,:,:,1), ...
    'intensityCross',intensities(:,:,:,2),'fitSettings',settings);
peakData.fitSettings.solver = struct('maxIter',options.MaxIter,'maxFunEvals',options.MaxFunEvals, ...
    'tolFun',1e-10,'tolX',1e-10,'robust',"Bisquare", ...
    'scaling',"angular max(abs(y)), then diagonal parameter scales; both undone in outputs");
end

function parameters = fitPeaks(x,y,locations,referenceWidths,settings,model,solver)
fitRange = settings.fitRange;
mask = x>fitRange(1) & x<fitRange(2);
x = x(mask); y = y(mask);
nPeaks = numel(locations);
nEdge = min(10,floor(numel(x)/4));
leftLevel = mean(y(1:nEdge));
backgroundSlope = (mean(y(end-nEdge+1:end))-leftLevel)/(x(end)-x(1));
xCenter = mean(fitRange);
backgroundLevel = leftLevel+backgroundSlope*(xCenter-x(1));
yRange = max(y)-min(y);
backgroundAtPeaks = backgroundLevel+backgroundSlope*(locations(:)-xCenter);
heights = max(interp1(x,y,locations(:),'nearest','extrap')-backgroundAtPeaks,0.02*yRange);
if isempty(referenceWidths)
    initialWidths = repmat(5,1,nPeaks);
    centerStep = settings.referenceCenterStep;
    lowerWidths = repmat(settings.widthLowerBound,1,nPeaks);
    upperWidths = repmat(settings.referenceWidthUpperBound,1,nPeaks);
else
    initialWidths = referenceWidths;
    centerStep = settings.polarCenterStep;
    lowerWidths = max(initialWidths*(1-settings.widthTolerance),settings.widthLowerBound);
    upperWidths = min(initialWidths*(1+settings.widthTolerance),max(30,1.2*max(initialWidths)));
end
initial = [heights(:).'.*locations.^2./initialWidths; locations; initialWidths];
lowerCenters = max(locations-centerStep,fitRange(1));
upperCenters = min(locations+centerStep,fitRange(2));
midpoints = (locations(1:end-1)+locations(2:end))/2;
upperCenters(1:end-1) = min(upperCenters(1:end-1),midpoints);
lowerCenters(2:end) = max(lowerCenters(2:end),midpoints);
lower = [zeros(1,nPeaks); lowerCenters; lowerWidths];
upper = [max(initial(1,:)*20,eps); upperCenters; upperWidths];
margin = max(yRange,eps);
slopeMargin = 2*margin/diff(fitRange);
solver.StartPoint = [initial(:).',backgroundLevel,backgroundSlope];
solver.Lower = [lower(:).',min(y)-margin,-slopeMargin];
solver.Upper = [upper(:).',max(y)+margin,slopeMargin];
fitScale = max(abs(y));
if isempty(referenceWidths) || fitScale==0
    fitScale = 1;
end
linearCoefficients = [1:3:3*nPeaks,3*nPeaks+(1:2)];
solver.StartPoint(linearCoefficients) = solver.StartPoint(linearCoefficients)/fitScale;
solver.Lower(linearCoefficients) = solver.Lower(linearCoefficients)/fitScale;
solver.Upper(linearCoefficients) = solver.Upper(linearCoefficients)/fitScale;
parameterScale = max(abs(solver.StartPoint),0.1*(solver.Upper-solver.Lower));
parameterScale(parameterScale<=0) = 1;
solver.StartPoint = solver.StartPoint./parameterScale;
solver.Lower = solver.Lower./parameterScale;
solver.Upper = solver.Upper./parameterScale;
fitted = fit(x,y/fitScale,model,solver,'problem',num2cell(parameterScale));
coefficients = coeffvalues(fitted).*parameterScale;
coefficients(linearCoefficients) = coefficients(linearCoefficients)*fitScale;
parameters = reshape(coefficients(1:3*nPeaks),3,[]);
end

