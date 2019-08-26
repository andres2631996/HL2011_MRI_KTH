% Exercise 2: The signal, Slice selection
% Andrés Martínez Mora
% HL2011: Magnetic Resonance Imaging

% Add seemri folder to path:

addpath('seemri');

% Initialize useful constants:
gammabar = 42.58e6;
gamma = 2*pi*gammabar;

% Initial parameters:
B0 = 2; % External magnetic field
tp = 1e-6; % Duration of a rectangular RF pulse
iv = ImagingVolume(0, 0, 0.5, 0.02); % Magnetization vector centered in (0,0) and T1= 0.5s and T2=0.02s

% RF pulses with different flip angles and flip directions:
% For a 90 degrees flip angle: B1=5.8713e-3 T
% For a 45 degrees flip angle: B1=2.9357e-3 T
% For a B1 pulse in x': phase=0 radians
% For a B1 pulse in -y': phase= pi/2 radians
% The RF pulse must have the same frequency than the Larmor frequency to
% allow for resonance and excitation
rf1=RectPulse(5.8713e-3,gammabar*B0,0,tp); % RF pulse: 90x'
rf2=RectPulse(5.8713e-3,gammabar*B0,pi/2,tp); % RF pulse: 90-y'
rf3=RectPulse(2.9357e-3,gammabar*B0,0,tp); % RF pulse: 45x'


% Send the pulses to the original magnetization vector and record the
% respective signals

[S1, ts] = seemri(iv, B0, rf1, [], [], ADC(0.2, 0.2/100)); % Send 90x'
iv.toEquilibrium();
[S2, ts] = seemri(iv, B0, rf2, [], [], ADC(0.2, 0.2/100)); % Send 90-y'
iv.toEquilibrium();
[S3, ts] = seemri(iv, B0, rf3, [], [], ADC(0.2, 0.2/100)); % Send 45x'
iv.toEquilibrium();

% Plot the resulting signals:
figure
subplot(2,1,1)
plot(ts, abs(S1), ts, abs(S2), '--', ts, abs(S3), '-.')
xlabel('Time (s)')
ylabel('Signal magnitude')
subplot(2,1,2)
plot(ts, angle(S1), ts, angle(S2), '--', ts, angle(S3), '-.'), ylim([0,pi])
xlabel('Time (s)')
ylabel('Signal phase (radians)');

%% Spin Echo
% Imaging volume with 9x9 magnetization vectors with T1=0.5s, T2=0.02s with
% a inhomogeneity of 1 microT

iv = ImagingVolume(-4:1:4, -4:1:4, 0.5, 0.02, 1, 'dB0Gamma', 1e-6);

B0 = 2; % External magnetic field
tp = 1e-6; % Pulse duration
phi_rf = [0 pi/2]; % Pulse phase

% Compute the necessary amplitude for having flip angles of 90 and 180
% degrees with the duration given. Compute as well the frequency of the RF
% pulse required for exciting the imaging volume

B1_90 = 5.8713e-3;
B1_180 =0.01174;
f_rf =gammabar*B0;

% Combine both 90 and 180 pulses into one variable
rf = RectPulse([B1_90 B1_180], f_rf, phi_rf, tp, [0 5e-3]);

% Apply both pulses to the created Imaging Volume
[S, ts] = seemri(iv, B0, rf, [], [], ADC(0.02, 0.02/100));

%% Gradient echo

B0 = 2; % External magnetic field
G = 1.2e-6; % Gradient amplitude
tau = 5e-3; % Duration of the inverted gradient
iv = ImagingVolume(-4:1:4, -4:1:4, 0.5, 0.02, 1, ...
'dB0Gamma', 0.1e-6); % Imaging volume of 9x9 with T1=0.5s, T2=0.02s and inhomogeneities of 1 microT

% Gradient to be applied: it is -G at time tp (after the RF pulse), then it
% is changed to 2G at tp+tau (being inverted during a time tau). At last it
% is turned off at tp+3tau, being non-inverted for a time of 2tau

Gx = Gradient([tp tp+tau tp+3*tau], [-G 2*G 0]);

% Apply the GE sequence to the Imaging Volume
[S, ts] = seemri(iv, B0, rf1, Gx, [], ADC(0.02, 0.02/100));
%% Slice selection: RF pulse frequency effect
B0 = 3; % External magnetic field
iv1 = ImagingVolume(-1:0.25:5, -1:0.25:1, 0.8, 0.07, 1, 'PlotScale', 2); 
% Magnetization vectors in an array of 10x10, T1=0.8s and T2=0.07s

tp = 2e-3; % Duration
B1 = 2.9e-6; % B1 amplitude
rf1 = SincPulse(B1, gammabar*B0, 0, tp); % Sinc RF pulse of amplitude B1, central frequency w0, phase 0 and duration tp

g1 = Gradient([0 tp 1.5*tp], [70e-6 -70e-6 0]); % Slice selection gradient 
[S, ts] = seemri(iv1, B0, rf1, g1, [], ADC(1.5*tp, tp/100)); % Apply RF pulse and slice-selection gradient to the given imaging volume


rf2 = SincPulse(B1, gammabar*B0+4e3, 0, tp); % RF pulse with an extra 4 kHz of central frequency: excite different points of the volume than RF1
rf3 = SincPulse(B1, gammabar*B0+8e3, 0, tp); % RF pulse with an extra 8 kHz of central frequency: excite different points of the volume than RF1

iv2 = ImagingVolume(-1:0.25:5, -1:0.25:1, 0.8, 0.07, 1, 'PlotScale', 2); % Imaging Volume where RF2 is applied
iv3 = ImagingVolume(-1:0.25:5, -1:0.25:1, 0.8, 0.07, 1, 'PlotScale', 2); % Imaging Volume where RF3 is applied
[S, ts] = seemri(iv2, B0, rf2, g1, [], ADC(1.5*tp, tp/100)); % Apply RF2 and the same gradient to Imaging Volume 2
[S, ts] = seemri(iv3, B0, rf3, g1, [], ADC(1.5*tp, tp/100)); % Apply RF3 and the same gradient to Imaging Volume 3

% Plot of the resulting Imaging Volumes
subplot(3,1,1);
plot(iv1)
subplot(3,1,2);
plot(iv2)
subplot(3,1,3);
plot(iv3)


%% Slice selection: Gradient effect
g2 = Gradient([0 tp 1.5*tp], [110e-6 -110e-6 0]); % Slice selection gradient much more intense than the one in the previous part
g3 = Gradient([0 tp 1.5*tp], [150e-6 -150e-6 0]); % Slice selection gradient much more intense than the one in the previous part and the previous one

% Set imaging volumes to equilibrium
iv1.toEquilibrium();
iv2.toEquilibrium();
iv3.toEquilibrium();
[S, ts] = seemri(iv1, B0, rf3, g1, [], ADC(1.5*tp, tp/100)); % Excite the Imaging Volume with the RF pulse with the extra 8 kHz and the 70e-6 gradient
[S, ts] = seemri(iv2, B0, rf3, g2, [], ADC(1.5*tp, tp/100)); % Excite the Imaging Volume with the RF pulse with the extra 8 kHz and the 110e-6 gradient
[S, ts] = seemri(iv3, B0, rf3, g3, [], ADC(1.5*tp, tp/100)); % Excite the Imaging Volume with the RF pulse with the extra 8 kHz and the 150e-6 gradient

% Plot the different imaging volumes
subplot(3,1,1);
plot(iv1)
subplot(3,1,2);
plot(iv2)
subplot(3,1,3);
plot(iv3)

%% Slice selection sequence design

tp = 2e-3; % Duration of the RF pulse
Delta_f = 12/tp; % Bandwidth of the RF pulse: six zero-crossings of the sinc on each side

Gss = 70.456*10^(-6); % Slice selection gradient to be applied for a slice thickness of 2mm
g = Gradient([0 tp 1.5*tp], [Gss -Gss 0]); % Variation in time of the chosen slice selection gradient

f_rf = 127.743*10^6; % Excitation frequency in order to have the slice centered in x0=1mm
rf = SincPulse(B1, f_rf, 0, tp); % Final sinc RF pulse

iv1.toEquilibrium();
[S, ts] = seemri(iv1, B0, rf, g, [], ADC(1.5*tp, tp/100)); % Apply the computed gradient on the Imaging Volume

% Apply now a sequence with both gradients in -x' and -y': oblique gradient
iv1.toEquilibrium();
[S, ts] = seemri(iv1, B0, rf, g, g, ADC(1.5*tp, tp/100));
