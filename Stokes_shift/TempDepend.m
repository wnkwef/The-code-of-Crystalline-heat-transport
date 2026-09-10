% =========================================================================
% WARNING: THIS SCRIPT CANNOT RUN DIRECTLY WITH DEMO DATA!
% The demo dataset contains insufficient data points, which will 
% trigger a model fitting error. Please replace it with the full dataset.
% =========================================================================
clear; close all; clc

% Load photoluminescence (PL, here Em) and photoluminescence excitation
% (PLE, here Ex) Data.
Data = EmExData();
NumT = Data.NumT;
peakX = Data.peakX;
muLower = Data.muLower;
muUpper = Data.muUpper;
ExIncpt = Data.ExIncpt;
ExErrLow = Data.ExErrLow;
ExErrHigh = Data.ExErrHigh;

f = figure(Position=[200, 100, 500, 600]);
ax = axes(f); hold on; grid minor;

% Plot bandgaps
yyaxis right
ylabel('Amplitude (eV)')
errorbar(NumT, ExIncpt, ExErrLow, ExErrHigh, '^--', ...
         'LineWidth', 1.5, 'CapSize',15, ...
         'DisplayName','Excitation Intcp', Color=[0.85 0.33 0.10]);
ylim([2.3,3.7])
yticks(2.3:0.1:3.7)

% Plot emission data
yyaxis left
ylabel('Amplitude (eV)')
errLow  = peakX - muLower;
errHigh = muUpper - peakX;
errorbar(NumT, peakX, errLow, errHigh, 'o--',...
         'LineWidth', 1.5, 'CapSize',12, ...
         'DisplayName','Emission Peaks w/ Error', MarkerFaceColor='b', ...
         Color=[0.00 0.45 0.70]);

% Plot and fit the Stokes shift
SS = ExIncpt - peakX;
SSHigh = sqrt(ExErrLow.^2 + errHigh.^2);
SSLow = sqrt(ExErrHigh.^2 + errLow.^2);
errorbar(NumT, SS, SSLow, SSHigh, '^--', ...
         'LineWidth', 1.5, 'CapSize',15, ...
         'DisplayName','Stokes Shift w/ Error', Color=[0.80 0.40 0.60]);
Tthr = 35; Thigh = 170;
Tthr_Arre = 155; Thigh_Arre = 290;
sigmaSS = max(abs(SSLow), abs(SSHigh));
fitresult_PD = fitSS_Point_dipole(NumT, SS, Tthr, Thigh);
fitresult_Arre = fitArre(NumT, SS, Tthr_Arre, Thigh_Arre, sigmaSS);

[Tcross, rootResidual, rootExitFlag] = fzero( ...
    @(T) feval(fitresult_PD,T) - feval(fitresult_Arre,T), [Tthr,Thigh]);
T_PD = linspace(Tthr, Tcross, 1000);
T_Arre = linspace(Tcross, Thigh_Arre, 1000);
plot(T_PD, feval(fitresult_PD,T_PD), '-', 'LineWidth',2, ...
    DisplayName='Point dipole model', Color=[0.2 0.15 0.8]);
plot(T_Arre, feval(fitresult_Arre,T_Arre), '-', 'LineWidth',2, ...
    DisplayName='Exponential fitting', Color=[0.36, 0.04, 0.16]);

xlabel('Temperature (K)');
ylabel('Amplitude (eV)');

ylim([1.0,2.4])
yticks(1.0:0.1:2.4)
legend(Location="best", FontSize=10)
hold off

%% Functions
function fitresult = fitSS_Point_dipole(NumT, SS, Tthr, Thigh)
%% Dielectric solvation model for point dipole

    Ebgd = 1.35942;         % Background energy (eV)
    epsilon_inf = 4;        % High-frequency dieletric constant
    p = 10;                 % Dipole
    a = 100;                % Radius of the cell
    v = a^3;                % Volume of the cell
    A = 20;                 % Pre-factor
    Ea = 30;                % Activation energy (meV)
    k = 8.6173e-5;          % Boltzmann constant (eV/K)
    mu = -35;               % Temperature shift

    % Startpoints
    SPs = [Ebgd, p, v, epsilon_inf, A, Ea, mu];
    ft = fittype(@(Ebgd, p, v, epsilon_inf, A, Ea, mu, T) ...
        Ebgd + 27.21 * (8*pi*p.^2./(3*v)) .* ( ...
                ((epsilon_inf + (A./(exp(Ea/1e3./(k*(T+mu))) - 1))) - 1) ./ ...
                  (2*(epsilon_inf + (A./(exp(Ea/1e3./(k*(T+mu))) - 1))) + 1) ...
              - (epsilon_inf - 1) ./ (2*epsilon_inf + 1) ...
            ), ...
    independent='T');

    mask = (NumT >= Tthr) & (NumT <= Thigh);
    NumTfit = NumT(mask);
    SSfit = SS(mask);

    me = 120; st = 30;
    weights = exp(-((NumTfit - me).^2) / (2 * st^2));
    weights = weights / max(weights);

    opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                      'StartPoint', SPs, ...
                      'MaxIter', 1e18, ...
                      'TolFun', 1e-18, ...
                      'TolX', 1e-18, ...
                      'Lower', [1.35942, 0, 0, 0, 0, 0, -35], ...
                      'Upper', [1.35942, Inf, 300^3, Inf, 1e2, 100, -35], ...
                      'Weights', weights);
    fitresult = fit(NumTfit(:), SSfit(:), ft, opts);
end

function fitresult = fitArre(NumT, SS, Tthr, Thigh, sigmaSS)
%% Arrhenius-type fitting
% SS(T) = Ebgd - A * exp[-Ea/(kB*T)]
% Ea in meV, and kB in eV/K.

    mask = (NumT >= Tthr) & (NumT <= Thigh);
    NumTfit = NumT(mask);
    SSfit = SS(mask);
    sigmaSSfit = sigmaSS(mask);

    kB = 8.617333262e-5;    % Boltzmann constant (eV/K)
    Ebgd = 1.719;           % Background energy (eV)
    A = 6;                  % Pre-factor (eV)
    Ea = 80;                % Activation energy (meV)

    ft = fittype(@(Ebgd, A, Ea, T) ...
        Ebgd - A.*exp(-(Ea/1e3)./(kB.*T)), ...
        independent='T');

    weights = 1 ./ (sigmaSSfit.^2);
    weights = weights / max(weights);

    % Startpoints
    SPs = [Ebgd, A, Ea];
    opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                  'StartPoint', SPs, ...
                  'MaxIter', 1e15, ...
                  'TolFun', 1e-18, ...
                  'TolX', 1e-18, ...
                  'Lower', [1.6, 0, 0], ...
                  'Upper', [1.8, 10, 100], ...
                  'Weights', weights);
    fitresult = fit(NumTfit(:), SSfit(:), ft, opts);
end
