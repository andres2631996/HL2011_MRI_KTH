% Exercise 4: Parallel Imaging
% HL 2011 Magnetic Resonance Imaging
% Andrés Martínez Mora

% Add seemri to path
addpath('seemri');

% Define some useful constants
gammabar = 42.58e6;
gamma = 2*pi*gammabar;

%% Brain Imaging

% Resolution of 4mm
B0 = 1.5; % External magnetic field
tp = 1e-3; % Pulse duration
B1 = (pi/2)/(gamma*tp); % RF pulse amplitude
f_rf = gammabar*B0; % RF pulse frequency
rf = RectPulse(B1, f_rf, 0, tp); % RF pulse to be applied: 90x'


TE = 10e-3; % Echo time
TR = 2; % Repetition time

tau = 4e-3; % Phase encoding gradient duration

FOV=220; % FOV size
pixel_size=4; % Pixel size


% k-space sampling parameters for the Gradient Echo Sequence
kmaxx = 1/(2*pixel_size);
kmaxy = 1/(2*pixel_size);
dkx = 1/FOV;
dky = 1/FOV;
dkx = kmaxx/ceil(kmaxx/dkx);
dky = kmaxy/ceil(kmaxy/dky);
kxa = -kmaxx:dkx:kmaxx-dkx; % Array
kya = -kmaxy:dky:kmaxy-dky; % Array
[kxg, kyg] = meshgrid(kxa,kya); % Grid


% Phase encoding gradient amplitudes
Gpexs =kxa./(gammabar*tau); 

% Applied phase-encoding gradient
gx = Gradient([tp tp+tau], {Gpexs 0});

% Frequency-encoding gradient amplitudes and applied frequency-encoding
% gradient
Gfey1 =-kmaxy/(gammabar*tau);
Gfey2 = kmaxy/(gammabar*tau);
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);

% Sampling time
dt = 1.4286e-04; % Time step
adc = ADC(TE-tau, TE+tau, dt); % Acquisition times

% Application of the Gradient Echo sequence

[S, ts] = brain_4mm_pixel(B0, rf, gx, gy, adc, TR, length(Gpexs));

% Image reconstruction and display
figure
mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true),title('Reconstructed image with a pixel size of 4mm')



% Resolution of 2mm
B0 = 1.5; % External magnetic field
tp = 1e-3; % Pulse duration
B1 = (pi/2)/(gamma*tp); % RF pulse amplitude
f_rf = gammabar*B0; % RF pulse frequency
rf = RectPulse(B1, f_rf, 0, tp); % RF pulse to be applied: 90x'


TE = 10e-3; % Echo time
TR = 2; % Repetition time

tau = 4e-3; % Phase encoding gradient duration

FOV=220; % FOV size
pixel_size=2; % Pixel size


% k-space sampling parameters for the Gradient Echo Sequence
kmaxx = 1/(2*pixel_size);
kmaxy = 1/(2*pixel_size);
dkx = 1/FOV;
dky = 1/FOV;
dkx = kmaxx/ceil(kmaxx/dkx);
dky = kmaxy/ceil(kmaxy/dky);
kxa = -kmaxx:dkx:kmaxx-dkx; % Array
kya = -kmaxy:dky:kmaxy-dky; % Array
[kxg, kyg] = meshgrid(kxa,kya); % Grid


% Phase encoding gradient amplitudes
Gpexs =kxa./(gammabar*tau); 

% Applied phase-encoding gradient
gx = Gradient([tp tp+tau], {Gpexs 0});

% Frequency-encoding gradient amplitudes and applied frequency-encoding
% gradient
Gfey1 =-kmaxy/(gammabar*tau);
Gfey2 = kmaxy/(gammabar*tau);
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);

% Sampling time
dt = 7.2072e-5; % Time step
adc = ADC(TE-tau, TE+tau, dt); % Acquisition times

% Application of the Gradient Echo sequence

[S, ts] = brain_2mm_pixel(B0, rf, gx, gy, adc, TR, length(Gpexs));

% Image reconstruction and display
figure
mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true),title('Reconstructed image with a pixel size of 2mm')

% Resolution of 0.5mm
B0 = 1.5; % External magnetic field
tp = 1e-3; % Pulse duration
B1 = (pi/2)/(gamma*tp); % RF pulse amplitude
f_rf = gammabar*B0; % RF pulse frequency
rf = RectPulse(B1, f_rf, 0, tp); % RF pulse to be applied: 90x'


TE = 10e-3; % Echo time
TR = 2; % Repetition time

tau = 4e-3; % Phase encoding gradient duration

FOV=220; % FOV size
pixel_size=0.5; % Pixel size


% k-space sampling parameters for the Gradient Echo Sequence
kmaxx = 1/(2*pixel_size);
kmaxy = 1/(2*pixel_size);
dkx = 1/FOV;
dky = 1/FOV;
dkx = kmaxx/ceil(kmaxx/dkx);
dky = kmaxy/ceil(kmaxy/dky);
kxa = -kmaxx:dkx:kmaxx-dkx; % Array
kya = -kmaxy:dky:kmaxy-dky; % Array
[kxg, kyg] = meshgrid(kxa,kya); % Grid


% Phase encoding gradient amplitudes
Gpexs =kxa./(gammabar*tau); 

% Applied phase-encoding gradient
gx = Gradient([tp tp+tau], {Gpexs 0});

% Frequency-encoding gradient amplitudes and applied frequency-encoding
% gradient
Gfey1 =-kmaxy/(gammabar*tau);
Gfey2 = kmaxy/(gammabar*tau);
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);

% Sampling time
dt = 1.8018e-05; % Time step
adc = ADC(TE-tau, TE+tau, dt); % Acquisition times

% Application of the Gradient Echo sequence

%[S, ts] = brain_halfmm_pixel(B0, rf, gx, gy, adc, TR, length(Gpexs));

% Image reconstruction and display
figure
%mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true),title('Reconstructed image with a pixel size of 0.5mm')


%% Parallel Imaging

% Add new folder to path
addpath('brain_4mm_pixel_coil');

% Resolution of 4mm
B0 = 1.5; % External magnetic field
tp = 1e-3; % Pulse duration
B1 = (pi/2)/(gamma*tp); % RF pulse amplitude
f_rf = gammabar*B0; % RF pulse frequency
rf = RectPulse(B1, f_rf, 0, tp); % RF pulse to be applied: 90x'


TE = 10e-3; % Echo time
TR = 2; % Repetition time

tau = 4e-3; % Phase encoding gradient duration

FOV=220; % FOV size
pixel_size=4; % Pixel size


% k-space sampling parameters for the Gradient Echo Sequence
kmaxx = 1/(2*pixel_size);
kmaxy = 1/(2*pixel_size);
dkx = 1/FOV;
dky = 1/FOV;
dkx = kmaxx/ceil(kmaxx/dkx);
dky = kmaxy/ceil(kmaxy/dky);
kxa = -kmaxx:dkx:kmaxx-dkx; % Array
kya = -kmaxy:dky:kmaxy-dky; % Array
[kxg, kyg] = meshgrid(kxa,kya); % Grid


% Phase encoding gradient amplitudes
Gpexs =kxa./(gammabar*tau);

% Reduce phase-encoding steps by half

Gpexs=Gpexs(1:2:end);

% Applied phase-encoding gradient
gx = Gradient([tp tp+tau], {Gpexs 0});

% Frequency-encoding gradient amplitudes and applied frequency-encoding
% gradient
Gfey1 =-kmaxy/(gammabar*tau);
Gfey2 = kmaxy/(gammabar*tau);
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);

% Sampling time
dt = 1.4286e-04; % Time step
adc = ADC(TE-tau, TE+tau, dt); % Acquisition times

% Application of the Gradient Echo sequence

[S, ts, C] = brain_4mm_pixel_coil(B0, rf, gx, gy, adc, TR, length(Gpexs));

% Image reconstruction and display
figure
[nx,ny,nc] = size(S);
for c = 1:nc
    img(:,:,c) = mrireconstruct(S(:,:,c), max(kmaxx,kmaxy), 'Plot', false);
    subplot(2,2,2*c-1),imshow(img(:,:,c),[]),title(strcat('Image obtained from coil #',num2str(c)));
    subplot(2,2,2*c),imshow(C(:,:,c),[]),title(strcat('Sensitivity map from coil #',num2str(c)));
end
% Obtain sum-of-squares reconstructed image:
sos=sum(img,3);
figure
imshow(sos,[]),title('Reconstructed sum-of-squares image');


% SENSE reconstruction
img = fftshift(img,2);
sense_img = [];
for jj=1:nx
    for ii=1:ny
        AB = [2*ii,2*ii+1]; % size(AB) = [1 2]
        S = [C(jj,ii,1),C(jj,ii,2)]; % size(S) = [2 2]
        P = [img(jj,ii,1),img(jj,ii,2)]'; % size(P) = [2 1]
        sense_img(jj,AB)= pinv(S')*P;
    end
end

figure
imshow(sense_img,[]),title('SENSE reconstructed image')