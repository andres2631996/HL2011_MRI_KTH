% Exercise 5: Image Contrast & Chemical Shift
% HL 2011 Magnetic Resonance Imaging
% Andrés Martínez Mora

% Add seemri to path
addpath('seemri');

% Define some useful constants
gammabar = 42.58e6;
gamma = 2*pi*gammabar;

% Take the sequence from Exercise 4 with a 4mm pixel size

B0 = 1.5; % External magnetic field
tp = 1e-3; % Pulse duration
B1 = (pi/2)/(gamma*tp); % RF pulse amplitude
f_rf = gammabar*B0; % RF pulse frequency
rf = RectPulse(B1, f_rf, 0, tp); % RF pulse to be applied: 90x'



TR = 5; % Repetition time

tau = 2e-3; % Phase encoding gradient duration

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

% Sampling time
dt = 2*tau/length(kya); % Time step

figure
cont=0;
c=[];
for TE=[5,50,100]*10^(-3)
    gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);


    adc = ADC(TE-tau, TE+tau, dt); % Acquisition times

    % Application of the Gradient Echo sequence

    [S, ts] = brain_4mm_pixel(B0, rf, gx, gy, adc, TR, length(Gpexs));

    % Image reconstruction and display
    cont=cont+1;
    subplot(1,3,cont),mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true),title(strcat('Reconstructed image with TE=',num2str(TE/10^(-3)),'ms'))
    a=mrireconstruct(S, max(kmaxx,kmaxy));
    c(cont)=abs(a(21,22)-a(41,43))/a(41,43);
end

%% Changing TR with fixed TE

% Resolution of 4mm
B0 = 1.5; % External magnetic field
tp = 1e-3; % Pulse duration
B1 = (pi/2)/(gamma*tp); % RF pulse amplitude
f_rf = gammabar*B0; % RF pulse frequency
rf = RectPulse(B1, f_rf, 0, tp); % RF pulse to be applied: 90x'


TE = 5e-3; % Echo time


tau = 2e-3; % Phase encoding gradient duration

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
dt = 2*tau/length(kya); % Time step
adc = ADC(TE-tau, TE+tau, dt); % Acquisition times
cont=0;
d=[];
figure
% Application of the Gradient Echo sequence
for TR = [2,1,0.5] % Repetition time
    [S, ts] = brain_4mm_pixel(B0, rf, gx, gy, adc, TR, length(Gpexs));
    cont=cont+1;
    % Image reconstruction and display
    subplot(1,3,cont),mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true),title(strcat('Reconstructed image with TR=',num2str(TR),'s'))
    b=mrireconstruct(S, max(kmaxx,kmaxy));
    d(cont)=abs(b(21,22)-b(41,43))/b(41,43);
end


%% Chemical shift

B0 = 1.5; % External magnetic field
tp = 1e-3; % Pulse duration
B1 = (pi/2)/(gamma*tp); % RF pulse amplitude
f_rf = gammabar*B0; % RF pulse frequency
rf = RectPulse(B1, f_rf, 0, tp); % RF pulse to be applied: 90x'

TR=2; % Repetition time
TE = 30e-3; % Echo time


FOV=20; % FOV size
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

cont=0;
d=[];

% Application of the Gradient Echo sequence
for tau=[2,4,10]*10^(-3) % Repetition time
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
    dt = 2*tau/length(kya); % Time step
    adc = ADC(TE-tau, TE+tau, dt); % Acquisition times
    
    M = waterandfat(FOV/4,pixel_size);
    cont=cont+1;
    % Image reconstruction and display
    figure
    [S,ts]=seemri(M,B0,rf,gx,gy,adc,TR,length(kya));
    figure
    mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true), title(strcat('Water and fat phantom with\tau=',num2str(tau/10^(-3)),'ms'))
end

%% Echo Planar Imaging (EPI)

tau = 0.5e-3; % ms
iv = disc(3,1); % Disc of radius 3mm and pixel size of 1mm
kmax = 1/2; % Maximum spatial frequency
dky = 1/8; % Frequency increment in phase encoding direction
dkx = kmax/ceil(kmax/dky); % Frequency increment in frequency encoding direction 
B0 = 1.5; % External magnetic field in T
tp = 1e-3; % Duration of RF pulse
B1 = pi/(2*gamma*tp); % Amplitude of RF pulse
f_rf = gammabar*B0; % Frequency of RF pulse
rf = RectPulse(B1, f_rf, 0, tp); % RF pulse

% Phase prewinder
Gpex1 = [-kmax/(gammabar*tau) 0]; % Phase-encoding gradient amplitudes during prewinder
Tpex1 = [tp tp+tau]; % Time instants for the prewinder

% Frequency prewinder
Gfey1 = [-kmax/(gammabar*tau) 0]; % Frequency-encoding gradient amplitudes during prewinder
Tfey1 = [tp tp+tau];
% First readout
Gfey2 = [kmax/(gammabar*tau) 0];
Tfey2 = [tp+tau tp+3*tau];
dt = dky/(gammabar*Gfey2(1));

Gblip = []; % Leave this empty
Tblip = []; % Leave this empty
cont=0;
for ll = 2:2*kmax/dky
    cont=cont+1;
    Gblip = cat(2,Gblip,[kmax/(gammabar*tau) 0]);
    Tblip = cat(2,Tblip, [Tfey2(end) Tfey2(end)+dt]);
    Gfey2 = cat(2,Gfey2,[((-1)^(ll-1))*kmax/(gammabar*tau) 0]);
    Tfey2 = cat(2,Tfey2, [Tblip(end) Tblip(end)+2*tau]);
end
gx = Gradient([Tpex1 Tblip], [Gpex1 Gblip]);
gy = Gradient([Tfey1 Tfey2], [Gfey1 Gfey2]);

adc = ADC(Tfey2(1), Tfey2(end)+dt, dt);
[S, ts] = seemri(iv, B0, rf, gx, gy, adc, Tfey2(end), 1,'PlotKspace',true);

% Read in the positive frequency encoding direction those unread lines
S = reshape(S,[2*kmax/dkx+1 2*kmax/dkx]);
S(:,2:2:end) = flip(S(:,2:2:end)); % Flip every second k-space line
S(end,:) = []; % Cut away the blip samples
figure
mrireconstruct(S, kmax, 'Plot', true);


%% Echo Planar Imaging (EPI) on brain image

tau = 0.5e-3; % ms

FOV = 220;
Delta_w = 4;
kmax = 0.5/Delta_w;
dk = 1/FOV;
dk = kmax/ceil(kmax/dk);
ks = -kmax:dk:kmax-dk;

B0 = 1.5; % T
tp = 1e-3; % ms

B1 = pi/(2*gamma*tp);
f_rf = gammabar*B0;
rf = RectPulse(B1, f_rf, 0, tp);

% Phase prewinder
Gpex1 = [-kmax/(gammabar*tau) 0];
Tpex1 = [tp tp+tau];

% Frequency prewinder
G=kmax/(gammabar*tau);
Gfey1 = [-G 0]; 
Tfey1 = [tp tp+tau];

% First readout
Gfey2 = [G 0];
Tfey2 = [tp+tau tp+3*tau];

dt = dk/(gammabar*G);

Gb = dk/(gammabar*dt);

Gblip = []; % Leave this empty
Tblip = []; % Leave this empty

for ll = 2:2*kmax/dk
    Gblip = cat(2,Gblip,[Gb 0]);
    Tblip = cat(2,Tblip, [Tfey2(end) Tfey2(end)+dt]);
    Gfey2 = cat(2,Gfey2,[G*(-1)^(ll-1) 0]);
    Tfey2 = cat(2,Tfey2, [Tblip(end) Tblip(end)+2*tau]);
end
    
gx = Gradient([Tpex1 Tblip], [Gpex1 Gblip]);
gy = Gradient([Tfey1 Tfey2], [Gfey1 Gfey2]);

adc = ADC(Tfey2(1), Tfey2(end)+dt, dt);

[S, ts] = brain_4mm_pixel(B0, rf, gx, gy, adc, Tfey2(end), 1);
S = reshape(S,[2*kmax/dk+1, 2*kmax/dk]);
S(:,2:2:end) = flip(S(:,2:2:end)); % Flip every second k-space line
S(end,:) = []; % Cut away the blip samples
mrireconstruct(S, kmax, 'Plot', true);

