% Exercise 1: RF Excitation of the Bulk Magnetization
% HL2011: Magnetic Resonance Imaging
% Andrés Martínez Mora

addpath('seemri') % Include seemri folder to work with

% Useful constants
gammabar = 42.58e6;
gamma = 2*pi*gammabar;

% Creating a magnetization vector
iv = ImagingVolume(0, 0, 1e9, 1e9); % Vector at (0,0) with T1=T2=1e+9
figure
plot(iv)

% Creating an RF pulse
rf = RectPulse(5.9e-6, 42.58e6, pi/2, 1e-3); % RF pulse with B1=5.9 microT,
% a frequency of 42.58 MHz, a phase of pi/2 radians and a time width of 1
% ms.
figure
subplot(1,2,1), plot(rf), title('Envelope function of the RF pulse')
subplot(1,2,2), powspec(rf), title('Envelope power spectrum of the RF pulse')

% Apply the RF pulse to the magnetization vector
figure
seemri(iv, 1.0, rf); % Set B0=1T
iv.toEquilibrium();

% Question 1: influence of the amplitude of B1 on the flip angle
B1=[2,4,6,8,10]*(10^(-6)); % Different B1 amplitudes to test

for i=B1
    % RF pulse with varying B1 amplitude, frequency of 42.58 MHz, 
    % phase of pi/2 radians and time width of 1 ms
    rf = RectPulse(i, 42.58e6, pi/2, 1e-3);
    % Set B0=1T. Leave the same initial magnetization vector as in Section 2
    figure
    seemri(iv, 1.0, rf); 
    iv.toEquilibrium();
    
end

% Question 2: influence of the phase of B1 on the excitation process
phase=[0,pi/2,pi,3*pi/2]; % Phases to be tested, in radians

for j=phase
    % RF pulse inducing a flip angle of 90 degrees (5.8713 microT, 42.58
    % MHz frequency, time-width of 1 ms)
    rf=RectPulse(5.8713e-06, 42.58e6, j, 1e-3);
    % Set B0=1T. Leave the same initial magnetization vector as in Section 2
    figure
    seemri(iv, 1.0, rf); 
    iv.toEquilibrium();
    
end

% Pulse test from Question 4
rf_test=RectPulse(0.0029,42.58*10^6,3*pi/2,10^(-6));
% Test RF pulse: amplitude of 2.9 mT, frequency of 42.58 MHz, phase of 270
% degrees and duration of 1 microsecond

% Keep the same original magnetization from Section 2
seemri(iv, 1.0, rf_test);