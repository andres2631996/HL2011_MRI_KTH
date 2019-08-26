% Exercise 3: K-space sampling and basic imaging
% HL2011 Magnetic Resonance Imaging
% Andrés Martínez Mora

% Add seemri to path:
addpath('seemri');

% Define some useful constants:
gammabar = 42.58e6;
gamma = 2*pi*gammabar;

%% Initialization
% Define a disc as a MATLAB function object:
u = @(x,y) sqrt(x.^2+y.^2)<=3; % Radius 3

% Digital image of the disc
figure
x = -5:0.1:5;
y = -5:0.1:5;
[xg, yg] = meshgrid(x, y);
imagesc(x, y, u(xg, yg));
axis image
colormap gray
title('Image space')

% Define the Fourier Transform of the disc:
Fu = @(kx,ky) 3*besselj(1, 2*pi*(sqrt(kx.^2+ky.^2)+1e-9*(kx==0 & ky==0))*3)...
./(sqrt(kx.^2+ky.^2)+1e-9*(kx==0 & ky==0));

% Define a sampling of the k-space:
kmaxx = 5;
kmaxy = 5;
dkx = 0.1;
dky = 0.1;
kxa = -kmaxx:dkx:kmaxx-dkx; % Array
kya = -kmaxy:dky:kmaxy-dky; % Array
[kxg, kyg] = meshgrid(kxa,kya); % Grid

% Visualization of the sampled k-space:
figure
imagesc(kxa, kya, Fu(kxg, kyg));
axis image
colormap gray
title('k-space')

% Reconstruction of the sampled k-space:
figure
mrireconstruct(Fu(kxg, kyg), max(kmaxx,kmaxy), 'Plot', true)
title(sprintf('kmaxx = %g, kmaxy = %g, dkx = %g, dky = %g', kmaxx, kmaxy, dkx, dky))


% Varying k-space extent:
figure
cont=0;
dkx=[0.05,0.1,0.2];
dky=[0.05,0.1,0.2];
for i=dkx
    cont=cont+1;
    kxa = -kmaxx:i:kmaxx-i; % Array
    kya = -kmaxy:i:kmaxy-i; % Array
    [kxg, kyg] = meshgrid(kxa,kya); % Grid
    subplot(3,2,cont),imagesc(kxa, kya, Fu(kxg, kyg)),title(strcat('k-space extent=',num2str(i)));
    cont=cont+1;
    subplot(3,2,cont),mrireconstruct(Fu(kxg, kyg), max(kmaxx,kmaxy), 'Plot', true), title(strcat('Reconstruction with k-space extent=',num2str(i)));
end

% Varying k-space max. range

figure
cont=0;
dkx = 0.1;
dky = 0.1;
kmaxx=[2,4,10];
kmaxy=[2,4,10];
for i=kmaxx
    cont=cont+1;
    kxa = -i:dkx:i-dkx; % Array
    kya = -i:dky:i-dky; % Array
    [kxg, kyg] = meshgrid(kxa,kya); % Grid
    subplot(3,2,cont),imagesc(kxa, kya, Fu(kxg, kyg)),title(strcat('k max=',num2str(i)));
    cont=cont+1;
    subplot(3,2,cont),mrireconstruct(Fu(kxg, kyg), i, 'Plot', true), title(strcat('Reconstruction with k max=',num2str(i)));
end
%% Basic Gradient Echo Imaging

B0 = 1.5; % External magnetic field
tp = 1e-3; % Pulse duration
B1 = (pi/2)/(gamma*tp); % RF pulse amplitude
f_rf = gammabar*B0; % RF pulse frequency
rf = RectPulse(B1, f_rf, 0, tp); % RF pulse to be applied: 90x'

iv = disc(3,1); % Imaging Volume where to apply the sequence: radius 3 circle

TE = 10e-3; % Echo time
TR = 2; % Repetition time

tau = 4e-3; % Phase encoding gradient duration

FOV=8; % FOV size
pixel_size=1; % Pixel size


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

% k-space sampling representation
figure
scatter(kxg(:),kyg(:),'r'),title('k-space sampling')
hold on
for i=kxa
    hold on
    plot(i+zeros(size(kxa)),kya,'b')
end
hold off

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
dt = 10^(-3); % Time step
adc = ADC(TE-tau, TE+tau, dt); % Acquisition times

% Application of the Gradient Echo sequence

[S, ts] = seemri(iv, B0, rf, gx, gy, adc, TR, length(Gpexs), 'PlotKSpace', true);


% Image reconstruction:
figure
subplot(1,2,1),imagesc(reshape(iv.Mz0, 7, 7)),title('Nominal image')
subplot(1,2,2),mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true),title('Reconstructed image');

