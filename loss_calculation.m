% Sean Keenan
% Heriot Watt 4th-Year
% Fabry-Perot loss measurement for integrated waveguides

clc
close all
clear variables

%% Load data and Generate model

% load data to analyse
[Data,Dir] = uigetfile('.txt');

% warn if user cancel 
if isequal(Data,0)
    warning('User cancelled operation')
    return;
else
end

FPdata = load(strcat(Dir,Data));
intensity = FPdata(:,1).';
lambda_axis = FPdata(:,2).';
lambda = lambda_axis*1E-9;
error = FPdata(:,3).';

% locates local maxima and minima, seperated by user defined value
TF1 = islocalmax(intensity,'MinSeparation',0.075E-9,'SamplePoints',lambda);
TF2 = islocalmin(intensity,'MinSeparation',0.075E-9,'SamplePoints',lambda);

% cavity details
% Si3N4
% str = ('Si_{3}N_{4}');
% L =  0.63;  
% L_2 = 6.3E-3;
% refractive index other WG (SiNi)
% nt = 1.9977;
% ni = 1;
% 800nm WG
% neff = 1.594;
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

% SiO2
str = ('SiO_{2}');
L = 0.532;
L_2 = 5.32E-3;
% refractive index TEEM WG (Si) 
% use https://www.siio.eu/eims.html to solve for neff
ni = 1;
nt = 1.53;
neff = 1.5178;

% loss co-efficient (dB/cm)
know_loss_db = 0.12;
alpha_known = 0.12/4.343;
% input power
Pin = 1;

% reflectivity co-efficient
r = (nt - ni)/(nt + ni);
% delta phase shift calculation
delta = (4*pi()*L_2*neff)./lambda;
% power output calculation
Pnumerator = Pin * (1 - (r^2 * exp(-alpha_known*L_2*2)))^2;
Pdenominator = 1 + (r^4 * exp(-alpha_known*4*L_2)) - (2 * r^2 * exp(-alpha_known*2*L_2) .* cos(delta));
Pout = Pnumerator ./ Pdenominator;

% choose values to use
max_index = [1,2,3,4];
min_index = [2,3,4,5];

%% Analyse Data and Calculate alpha 

% find max and min points
Imax = intensity(TF1);
Imin = intensity(TF2);
T_Imax = max(intensity(TF1));
T_Imin = min(intensity(TF2));
% determine what/where first point is
if length(Imax) > length(Imin)
    limit = length(Imin);
else
    limit = length(Imax);
end

% pre-assign vector
loss = zeros(size(max_index));
% loss = zeros(1:length(limit));
% calculate alpha from pairs of max & min points
% for n = 1:limit
for n = 1:length(max_index)
    
%     alpha_num = sqrt(Imax(n)./Imin(n)) + 1;
%     alpha_den = sqrt(Imax(n)./Imin(n)) - 1;
    alpha_num = sqrt(Imax(max_index(n))./Imin(min_index(n))) + 1;
    alpha_den = sqrt(Imax(max_index(n))./Imin(min_index(n))) - 1;
   
    loss(n) = 4.343*(1/(2*L))*log(r^2*(alpha_num/alpha_den));
end

% calculate average alpha and st.dev
loss_avg = mean(loss);
sdev = std(loss);
%calculate alpha from extremes
alpha_num = sqrt(T_Imax./T_Imin) + 1;
alpha_den = sqrt(T_Imax./T_Imin) - 1;
loss_ext = 4.343*(1/(2*L))*log(r^2*alpha_num/alpha_den);

%% Plot Model and Real data

% plots data
hold on
figure(1)

% FP experimental data plot
subplot(2,1,1)
plot(lambda_axis,intensity,lambda_axis(TF1),intensity(TF1),'r*',lambda_axis(TF2),intensity(TF2),'r*')
title(strcat(['5.32mm TEEM',' ',str, 'WG with Calculated Loss = ', num2str(loss_avg,'%.2f'),'\pm',num2str(sdev,'%.2f'),'dB/cm']),'FontSize',14)
ylabel ('Voltage (V)','FontSize',14)
xlabel ('Lambda (nm)','FontSize',14)
axis([min(lambda_axis) max(lambda_axis) min(intensity)-max(intensity)*0.01 max(intensity)+max(intensity)*0.01])

% FP theoretical data plot
subplot(2,1,2)
plot(lambda_axis, Pout(:,:,1))
title(strcat(['Theoretical Model for',' ',str, 'Waveguide with Loss = ', num2str(know_loss_db),'dB/cm']),'FontSize',15)
ylabel ('Intensity (Arb. Units)','FontSize',15)
xlabel ('Lambda (nm)','FontSize',15)
axis([min(lambda_axis) max(lambda_axis) min(Pout(:,:,1))-(max(Pout(1))*0.02) max(Pout(:,:,1))+(max(Pout(1))*0.02)])

% string output for calculations
str1 = 'Actual loss = %.3f dB/cm \n';
str2 = 'Calculated loss = %.3f dB/cm \n';
str3 = 'Sigma = %.3f \n';
str4 = strcat(str1, str2, str3);

fprintf(str4, know_loss_db, loss_avg, sdev);

