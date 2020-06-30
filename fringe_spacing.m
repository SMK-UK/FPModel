% Sean Keenan
% Heriot Watt 4th-Year
% Fabry-Perot loss measurement for integrated waveguides
% Illustration of change in fringe contrast for low-loss values

clc
clear figure
clear variables

% cavity length (cm) / (m)
L = 0.532;
L_2 = 5.32E-3;
% wavelength start (m)
l.start = 1550.48E-9;
% wavelength end (m)
l.end = 1550.62E-9;
% wavelentgth step size (m)
step = 0.002E-9;
% wavelength array (m)
lambda = l.start:step:l.end;
% loss co-efficient (dB/cm)
loss = [2 1 0.12];
% loss co-efficient (1/m)
alpha = loss*100/4.343;
% input power
Pin = 1;
% refractive index
ni = 1;
nt = 1.53;
neff = 1.5178;
% reflectivity co-efficient
r = (nt - ni )/(nt + ni);
% delta phase shift calculation
delta = (4*pi()*L_2*neff)./lambda;

dot = ['x' 'o' '*'];
xaxis = lambda*1E9;

% power output calculation
for j = 1:length(loss)
    
    Pnumerator = Pin * (1 - r^2 * exp(-alpha(j)*L_2))^2;
    Pdenominator = 1 + (r^4 * exp(-alpha(j)*4*L_2)) - (2*r^2* exp(-alpha(j)*2*L_2) .* cos(delta));
    Pout(:,:,j) = Pnumerator./Pdenominator;
    
    % calculate fringe difference
    Imax(j) = max(Pout(:,:,j));
    Imin(j) = min(Pout(:,:,j));
    
    if j >= 2
        errormax(j) = Imax(j) - Imax(j-1);
        errormin(j) = Imin(j-1) - Imin(j);
    end
    
    % plot for waveguide
    hold on
    subplot(6,1,1:2)
    plot(xaxis, Pout(:,:,j),'-','MarkerSize',3)
    hold off
    hold on
    subplot(9,1,5:9)
    plot(xaxis(49:63), Pout(:,49:63,j), dot(j),'MarkerSize',4)
    hold off
end 

% calculate fringe contrast
errormax(1) = Imax(3) - Imax(1);
errormin(1) = Imin(1) - Imin(3);
contrast = Imax./Imin;

s1 = {' '};

% legend for plot
for k = 1:length(loss)
    
    string(k) = strcat('\alpha = ',s1, num2str(loss(k)),'dB/cm');
    
end

% format for axis
ymin = min(min(Pout));
ymax = max(max(Pout));
ymin = ymin - ymin./20;
ymax = ymax + ymax./20;

% format plot
subplot(6,1,1:2)
title ('FP Spectrum over range of \lambda','FontSize',15)
legend(string,'location','best','FontSize',14)
xlabel ('Input Light Wavelength (nm)','FontSize',15)
ylabel ('Intensity (A.U)','FontSize',15)
yticks(0.85:0.05:1)
xticks(min(xaxis):0.02:max(axis))
axis([min(xaxis) max(xaxis) ymin ymax])

subplot(9,1,5:9)
title ('Close-up of Spectrum Spacing','FontSize',15)
legend(string,'location','south','FontSize',14)
xlabel ('Input Light Wavelength (nm)','FontSize',15)
ylabel ('Intensity (A.U)','FontSize',15)
yticks(0.98:0.005:1)
xticks(1550.56:0.01:1550.61)