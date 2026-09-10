function [centers, widths] = LorentzianFitWL(data, NumT, Temperatures, locations)
%% This function fits Raman peaks with Lorentzian functions.

fitRange = [0,300]; centerStep = 5; widthRange = [1e-6,30]; initialWidth = 5;
hcOverK = 6.62607015e-34 * 299792458 * 100 / 1.380649e-23; % K cm
[centers, widths] = deal(cell(numel(NumT),1));
for index = 1:numel(NumT)
    spectrum = data.(Temperatures{index});
    x = spectrum.wavenumber(:); y = spectrum.D0(:); p = locations{index}(:)';
    selected = x > fitRange(1) & x < fitRange(2);
    x = x(selected); y = y(selected);
    y = y.*(-expm1(-hcOverK*x/NumT(index)))./x;
    peak = (1:numel(p))'; heights = interp1(x,y,p,'nearest');
    start = [heights.*p.^2/initialWidth; p; initialWidth*ones(size(p))];
    lower = [0.05*start(1,:); max(p-centerStep,fitRange(1)); widthRange(1)*ones(size(p))];
    upper = [10*start(1,:); min(p+centerStep,fitRange(2)); widthRange(2)*ones(size(p))];
    terms = compose("c%d*g%d.^3./(x.^2.*g%d.^2+(x.^2-v%d.^2).^2)",peak,peak,peak,peak);
    names = [compose("c%d",peak), compose("v%d",peak), compose("g%d",peak)]';
    % Joint model of Lorentzian functions.
    model = fittype(strjoin(terms,'+'),'independent','x','coefficients',cellstr(names(:)));
    options = fitoptions('Method','NonlinearLeastSquares','StartPoint',start(:), ...
        'Lower',lower(:),'Upper',upper(:),'MaxIter',2000,'MaxFunEvals',20000, ...
        'TolFun',1e-10,'TolX',1e-10,'Display','off');
    fitted = fit(x,y,model,options);
    params = reshape(coeffvalues(fitted),3,[]);
    centers{index} = params(2,:)'; widths{index} = params(3,:)';
end
end
