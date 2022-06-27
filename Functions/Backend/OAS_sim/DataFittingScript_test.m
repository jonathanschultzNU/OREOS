%% This script will take you through the steps to calculate the oscillation
% associated spectrum of a system. This script assumes you have your data
% as a 2D matrix (wavelength x time), with accompanying wavelength and
% time vectors.

%   Author: Matthew S. Kirschner
%   Email: kirschner.21 (at) gmail.com
%   Last revision date: August 6, 2019
%
%   Copyright: Matthew S. Kirschner, 2019

clear all
close all

%% Here you can set your data to correspond with this script

addpath('C:\Users\schul\Documents\Graduate\Graduate Work\Fall 2017\Wasielewski Group\Data\2018.05.01')
dataload = load('20180501#02#001.SCAN');
data = dataload; %(2:end,2:end); %early(18).sub; % put in the matrix for the 2D spectrum
data(:,1) = []; data(end,:) = [];
data = data';
wavelength = 432.65 + (1:1:1024).*0.56503;       %dataload(2:end,1); %early(18).wavelengths; % put in the wavelength vector
time = -250:10:2000; %dataload(1,2:end); %early(18).time; % put in the time vector

%% Here we will set a number of the fit terms

first_wavelength = 525; % the first wavelength in the OAS
last_wavelength = 750; % the final wavelength you want in your OAS
first_time = 80; % the first time point for your Fourier Transform
last_time = 2000; % the final time point for your fit
time_resolution = 10; % the time resolution of your measurements, your data
% will be resampled/interpolated at this rate
criterion = -1; %Information criterion for itcmp, -1 uses AIC, -2 MDL 
% see itcmp for more info
oscillation_ratio = 2; % Threshold for defining which components are 
% oscillatory.  If the exponential damping rate is less than oscThresh 
% times the frequency, it is considered an oscillation. See 
% itcmpFilterOscillations for more details


%% Calculating the OAS


% We are performing a series of discrete fourier transformations, so 
% we need to specify the period range we are interested in. We perform 
% these measurements with respect to period, but it is trivial to invert
% these values if you have a frequency range in mind.

initial_period  = 10; % the shortest period you are interested in
final_period = 1000; % the longest period you are interested in
resolution = 1; % the resolution you want

% now we convert your input values to a matrix
periods = initial_period : resolution : final_period;

tic
% The actual OAS calculation
[Amplitude,Phase,w,ExponentialFits] = OscillationAssociatedSpectrum(data, ...
    wavelength, time, first_wavelength, last_wavelength, first_time,...
    last_time, time_resolution, periods, criterion, oscillation_ratio);
toc

%% The data should now be visualized to make sure the period range chosen 
% was appropriate and determine the wavelength range to use for
% lamba > lambda max


%initializing the plot
figure(1)

% Generates a subplot for FT vs period. The period range is suitable if
% there is a well defined peak rather than a sloped line or partial peak.
% If the period range needs to be changed, the OAS should be recalculated.
subplot(1,2,1)
plot(periods, mean(Amplitude), 'Linewidth', 2, 'Color', 'k')

% making it look nicer
xlabel('Period')
ylabel('FT Magnitude')
title('FT vs Period')
set(gca,'fontsize', 18);
xlim([periods(1) periods(end)])
ylim([0, max(mean(Amplitude)) * 1.1])

% Generates a subplot for FT vs wavelength. This will allow the
% determination of a range for lambda > lambdamax. This should include the
% maximum 10 ~ 20 values on the red side of the spectrum.
subplot(1,2,2)
plot(w,mean(Amplitude,2),'Linewidth',2,'Color','k')

% making it look nicer
xlabel('Wavelength')
ylabel('FT Magnitude')
title('FT vs Wavelength')
set(gca,'fontsize', 18);
xlim([w(1), w(end)])
ylim([0, max(mean(Amplitude,2))*1.1 ])


%% Now the reference phases are calculated and then used to generate the OAS

lambda_max = 530; % put in blue edge of your lambda > lambdamax region
lambda_max_max = 600; % put in red edge of your lambda > lambdamax region

% actual calculation of theta_0
[theta] = ReferencePhase(Phase,w,lambda_max,lambda_max_max);

% converting phase into a phase map
PhaseMap = cos(Phase - ones(length(w), 1) * theta);

% finally producing the 2D Oscillation Associated Spectrum
OAS = Amplitude .* PhaseMap;

%% Finally, a plot of the OAS can be generated

c = 2.9979 * 10^-5;      % speed of light [cm/fs]

% plotting the data with a contour plot
figure(2)
contourf(w, 1./(c.*periods), abs(OAS'),'edgecolor','none')
cmap = cmap2d(10);
colormap(cmap)
% h.LineWidth = 0;

% making it look nicer
ylabel('Period')
xlabel('Wavelength')
ylim([60 400])
set(gca, 'fontsize', 18);