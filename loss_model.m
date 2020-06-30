% Sean Keenan
% Heriot Watt 4th-Year
% Fabry-Perot loss measurement for integrated waveguides

clc
close all
clear variables

%% Variables for Model

% lambda range
lambda_start = 1550;
lambda_end = 1551;
resolution = 0.001;
lambda_axis = lambda_start:resolution:lambda_end;
lambda = lambda_axis*1E-9;

% cavity length (mm) / (m)
% SiO2
str1 = ('SiO_{2}');
L(1) = 5.32; 
L_2(1) = 5.32E-3; 
% Si3N4
str2 = ('Si_{3}N_{4}');
L(2) =  7; 
L_2(2) = 7E-3; 
% loss co-efficient SiO2 (dB/cm)
loss_db(1) = 0.12;
alpha_known(1) = 0.12/4.343;
% loss co-efficient Si3N4 (dB/cm) 
loss_db(2) = 2;
alpha_known(2) = 2/4.343;
% input power
Pin = 1;
% refractive index TEEM WG (Si) 
% use https://www.siio.eu/eims.html to solve for neff
ni(1) = 1;
nt(1) = 1.53;
neff(1) = 1.5178 ;
% refractive index other WG (SiNi)
nt(2) = 1.9977;
ni(2) = 1;
% 800nm WG
neff(2) = 1.594;
% 1050nm WG
% neff = 1.667;
% 1440nm WG
% neff = 1.720;
% 1520nm WG
% neff = 1.725;
% 1620nm WG
% neff = 1.7325;
% 1800nm WG
% neff = 1.742;

%% Generate model

for n = 1:length(L)
    
    % reflectivity co-efficient
    r = (nt(n) - ni(n))./(nt(n) + ni(n));
    % delta phase shift calculation
    delta = (4*pi()*L_2(n)*neff(n))./lambda;
    % power output calculation
    Pnumerator = Pin * ((1 - r^2) * exp(-alpha_known(n)*L_2(n)))^2;
    Pdenominator = 1 + (r^4 * exp(-alpha_known(n)*4*L_2(n))) - (2 * r^2 * exp(-alpha_known(n)*2*L_2(n)) .* cos(delta));
    Pout(:,:,n) = Pnumerator./Pdenominator;
end

%% Plot Model
% plots data
hold on
figure(1)

% FP experimental data plot
subplot(2,1,1)
plot(lambda_axis, Pout(:,:,1))
title(strcat(['Theoretical Model for',' ',str1, 'Waveguide with \alpha = ', num2str(loss_db(1)),'dB/cm']),'FontSize',15)
ylabel ('Intensity (Arb. Units)','FontSize',15)
xlabel ('Lambda (nm)','FontSize',15)
axis([min(lambda_axis) max(lambda_axis) min(Pout(:,:,1))-(max(Pout(1))*0.02) max(Pout(:,:,1))+(max(Pout(1))*0.02)])

% FP theoretical data plot
subplot(2,1,2)
plot(lambda_axis, Pout(:,:,2))
title(strcat(['Theoretical Model for',' ',str2, 'Waveguide with \alpha \approx', num2str(loss_db(2)),'dB/cm']),'FontSize',15)
ylabel ('Intensity (Arb. Units)','FontSize',15)
xlabel ('Lambda (nm)','FontSize',15)
axis([min(lambda_axis) max(lambda_axis) min(Pout(:,:,2))-(max(Pout(:,:,2))*0.02) max(Pout(:,:,2))+(max(Pout(:,:,2))*0.02)])
