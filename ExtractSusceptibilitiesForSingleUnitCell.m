clc 
close all; 
clear all;

%% Define some constants
I = sqrt(-1);
c = 3e8;

%% Define parameters 
% (This is the only part of the script where you should make any changes)

% Enter the dimensions of the unit cell for which you want to extract
% susceptabilities for: 
ringDiagonalLength = 2.3; % [mm]
ringWidth = 0.2;         % [mm]

% Enter the angles of incidence which you want to use for extracting the
% susceptibilities
theta1 = 0;    
theta2 = 30; 

% Enter the angle of incidence for which you want to reconstruct the unit
% cell for 
theta = 45; 

%% Extract the data from the csv files

% Search for the corresponding csv files 
fileDirectory1 = sprintf('S11-data/S11-data (%d deg incidence).csv', theta1);
fileDirectory2 = sprintf('S11-data/S11-data (%d deg incidence).csv', theta2);
data1 = readmatrix(fileDirectory1);
data2 = readmatrix(fileDirectory2);

ringDiagonalCol = 1; 
ringWidthCol = 2;
frequencyCol = 3;
realCol = 4;
imagCol = 5; 

% Extract frequency range and k0
[startRow, dummyCol] = find(data1(:, ringDiagonalCol) == ringDiagonalLength & data1(:, ringWidthCol) == ringWidth);
startRow = startRow(1);
numPoints = 401;
endRow = startRow + numPoints - 1;

frequency = data1(startRow:endRow, frequencyCol).*1e9;
lambda0 = c./frequency;
k0 = 2*pi./lambda0;

% Extract S11 data

R1 = data1(startRow:endRow, realCol) + I*data1(startRow:endRow, imagCol); 
R2 = data2(startRow:endRow, realCol) + I*data2(startRow:endRow, imagCol);

%% Perform extraction of chi terms
X_me_xy = 2*I./k0; 
X_ee_yy = (4*I./k0) .* ( cosd(theta2).*(R1+1).*(R2-1).*sind(theta1)^2 - sind(theta2)^2.*(R1-1).*(R2+1).*cosd(theta1) ) ...
    ./ ( (sind(theta1)^2 - sind(theta2)^2).*(R2+1).*(R1+1) );
X_mm_zz = -(4*I./k0) .* ( (R1-1).*(R2+1).*cosd(theta1) - cosd(theta2).*(R1+1).*(R2-1) ) ...
    ./ ( (sind(theta1)^2 - sind(theta2)^2).*(R2+1).*(R1+1) );

%% Reconstruct R for another angle 
R = ( 4.*cosd(theta) + I.*k0.*sind(theta)^2.*X_mm_zz - I.*k0.*X_ee_yy ) ... 
    ./ ( 4.*cosd(theta) - I.*k0.*sind(theta)^2.*X_mm_zz + I.*k0.*X_ee_yy );
R_dB = 20*log10(abs(R)); 
R_phase = angle(R) .* 180/pi; 

figure(1)
subplot(1,2,1);
plot(frequency/1e9, R_dB);
title('S_{11} Magnitude')
xlabel('Frequency [GHz]')
ylabel('S_{11} [dB]')

subplot(1,2,2);
plot(frequency/1e9, R_phase);
title('S_{11} Phase')
xlabel('Frequency [GHz]')
ylabel('S_{11} [Degrees]')