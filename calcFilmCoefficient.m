%% TITLE: Calculating Convective Film Coefficient In Terms of Engine Length
% Author: Evan Olsen
% Date: 01/25/2024
% Description: Uses pre-calculated interior contour and CEA transport
% properties to calculate heat transfer coefficient, according to the
% equation described in source: https://jffhmt.avestia.com/2015/PDF/006.pdf

clc; clear; close all;

%% INPUT DEFINITION
% R_x.mat is the pre-calculated interior contour of the engine, calculate
% this yourself or use the Parametric-Engine-Generator script on the GitHub
% at https://github.com/RIT-Launch-Initiative/Parametric-Liquid-Engine
load("R_x.mat")

% Radius of the combustion chamber
R_ch = 0.0254; % [m]

% Lateral surface area of the interior geometry revolved about the z-axis.
% You can find this easiest by using the measure tool in Fusion360 when you
% have imported the generated geometry.
A_lat = 0.0225786645; % [m^2]

% Calculates the local cross-sectional area of the engine based on the
% local radius.
A_x = pi.*R_x(2,:).^2; % [m^2]

% Mass flowrate
m_dot = 0.315; % [kg/s]

% Viscosity/Dynamic Viscosity as provided by NASA CEA as a special output.
% Ensure you convert cP to Pa*s before running.
mu = 0.00005832; % [Pa*s]

%% SEGMENTING GEOMETRY
% This code is meant to segment the geometry into the sections defined by
% NASA CEA so that the proper transport values can be applied in each
% region.

contour_disp_from_comb = R_x(2,:)-R_x(2,1);
i=1;
while contour_disp_from_comb(i)>=-0.0001
    comb_end = i;
    i = i+1;
end

[~,throat] = min(R_x(2,:));

nozzle_end = length(R_x(2,:));
%% INITAILIZING TRANSPORT VALUES
% EXAMPLE NASA CEA TRANSPORT VALUES:
  % WITH EQUILIBRIUM REACTIONS
  % 
  %   Cp(kJ/kg-K):     1.9857    1.9856    2.0068    6.4084
  %         Cond.:     2.1899    2.1895    1.9871    3.8802
  %       Prandtl:    0.52888   0.52889   0.53668   0.66488
  % 
  % WITH FROZEN REACTIONS
  % 
  %   Cp(kJ/kg-K):     1.9501    1.9501    1.9083    1.7877
  %         Cond.:     2.1593    2.1589    1.9216    1.3441
  %       Prandtl:    0.52678   0.52679   0.52771   0.53543

% C_p [J/(kg*K)] --- Specific Heat Constant Pressure
C_p_CC_eq = 1985.7; C_p_CCE_eq = 1985.6; C_p_TH_eq = 2006.8; C_p_NE_eq = 6408.4;
C_p_CC_fr = 1950.1; C_p_CCE_fr = 1950.1; C_p_TH_fr = 1908.3; C_p_NE_fr = 1787.7;

% Pr [-] --- Prantyl Number
Pr_CC_eq = 0.52888; Pr_CCE_eq = 0.52889; Pr_TH_eq = 0.53668; Pr_NE_eq = 0.66488;
Pr_CC_fr = 0.52678; Pr_CCE_fr = 0.52679; Pr_TH_fr = 0.52771; Pr_NE_fr = 0.53543;

%%  Eq. Before Throat, Fr. After Throat
% Found to be the most accurate method to produce heat transfer film
% coefficient: Uses equillibrium transport values prior to throat, and
% frozen transport values following the throat.
% 
% Linearly interpolates between each of the values to give the best
% interpretation of the values between each point.
C_p = [linspace(C_p_CC_eq,C_p_CCE_eq,comb_end),...
         linspace(C_p_CCE_eq,C_p_TH_fr,throat-comb_end)...
         linspace(C_p_TH_fr,C_p_NE_fr,nozzle_end-throat)];
Pr =  [linspace(Pr_CC_eq,Pr_CCE_eq,comb_end),...
         linspace(Pr_CCE_eq,Pr_TH_fr,throat-comb_end)...
         linspace(Pr_TH_fr,Pr_NE_fr,nozzle_end-throat)];

%% CALCULATION
% Solving equation as is described in https://jffhmt.avestia.com/2015/PDF/006.pdf

% Correction factor based on combustion chamber radius.
Z_c= pi*R_ch^2/A_lat; % [-]
h_film = Z_c.*(m_dot./(2.*A_x)).*C_p.*mu.^0.3.*Pr.^(2/3); % [W/(m^2*K)]
%% ANSYS FORMAT OUTPUT
% Inverses the x-values due to ANSYS coordinates being inversed. Change
% this value as represents your coordinate system in ANSYS.
% Creates vertical columns for easy import into ANSYS Transient Thermal
% Analysis.
h_film_x = [-R_x(1,:);h_film]'; % [m];[W/(m^2*K]
%% GRAPHICAL ANALYSIS
% Plots the heat transfer film coefficient in terms of x for visual
% representation and for comparison with the imported values in ANSYS.

filmCoefficientPlot = figure('Name','Film Coefficient');
subplot(2,1,1)
plot(R_x(1,:),h_film,'-r');
grid on
grid minor
title('Convective Film Coefficient vs. Length')
xlabel('L_e [m]');
ylabel('h_v, Film Coefficient [W/(m^2*K)]')
xlim([R_x(1,1),R_x(1,end)])

subplot(2,1,2)
plot(R_x(1,:),R_x(2,:),'-k')
grid on
grid minor
title('Interior Engine Contour')
xlabel('L_e [m]'); ylabel('R [m]')
xlim([R_x(1,1),R_x(1,end)])
axis equal;