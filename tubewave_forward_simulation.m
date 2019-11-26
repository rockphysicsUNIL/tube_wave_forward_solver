% Semi-analytical tube-wave forward simulator
% Written by Jürg Hunziker, University of Lausanne
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The code is based on the following two papers:  
% - Minato, S. and Ghose, R.: Low-frequency guided waves in a fluid-filled borehole: Simultaneous effects of generation and scattering due to
%   multiple fractures, Journal of applied Physics, 121, 104 902, 2017.
% - Minato, S., Ghose, R., Tsuji, T., Ikeda, M., and Onishi, K.: Hydraulic Properties of Closely Spaced Dipping Open Fractures Intersecting a
%   Fluid-Filled Borehole Derived From Tube Wave Generation and Scattering, Journal of Geophysical Research: Solid Earth, 122, 2017.
clear all; close all; clc; tic; 

nfrac = 3;

doCG = 0; % To determine the characteristic wavenumber kr, use conjugate gradient (1) or fminsearch (0)
doricker = 0; % Assign a Ricker wavelet (1) or load a wavelet from a file (0)
doclosepool = 0; % Close parallel pool at end of simulation (1) or not (0)
fc = 1500; % Center frequency of Ricker-wavelet [Hz]
freqvec_out = linspace(0,4095,4096); % Frequency vector for output [Hz] Start at 0 Hz!
freqvec = [linspace(1,4095,240)]; % Defines the frequencies, for which the forward problem is actually evaluated
% freqvec = [linspace(1,4095,4096)]; % Defines the frequencies, for which the forward problem is actually evaluated
R = 0.147/2; % Borehole radius [m]
theta = [68.0,62.0,77.0]; % fracture angle [°] (0° is a horizontal fracture)
rho = 2730; % Density of the surrounding rock [kg/m^3]
rhof = 1000; % Density of borehole fluid [kg/m^3]
mu = 26e9; % shear modulus of the formation [N/m^2]
mut = 6e9; % shear modulus for the tube-wave [N/m^2]
Ks = 28e9; % bulk modulus of the formation [N/m^2]
Kf = 2.251e9; % fluid bulk modulus [N/m^2]
L0 = [0.01,0.0001,0.0007]; % aperture of fracture [m]
Z = [0.37,0.96,0.53]*1e-14; % fracture compliance [m/Pa]
nu = 1e-6; % kinematic viscosity of the fluid [m^2/s]
zsrc = 1.6; % depth of the p-wave source [m]
zfrac = [23.8,23.6,25.0]; % depth of fractures [m]
zvec = linspace(17,37,201); % Depthvec along borehole
gsexp = 1.0; % geometrical spreading exponent
tshift = 0.002; % tube-wave attenuation shift factor [s] (between 0.001 and 0.02)
gsexptube = 784.1; % tube-wave attenuation exponent [-] (between 0.0 and 1000.0)
CoreNum=2; % Amount of CPU's to be used

dz = zvec(2)-zvec(1); % vertical spacing [m]
alpha_f = sqrt(Kf/rhof); % fluid velocity in the fracture [m/s]
Vs = sqrt(mu./rho); % S-wave velocity of the formation [m/s]
Vst = sqrt(mut./rho); % S-wave velocity for the tube-wave calculation [m/s]
Vp = sqrt((Ks+4/3.*mu)./rho); % P-wave velocity of the formation [m/s]
inv_alpha_eff = sqrt(1./(alpha_f.^(2))+rhof.*Z./L0); % inverse of the effective fluid velocity in the fracture (equation A3)
cT = alpha_f./sqrt(1+rhof.*alpha_f.^2./(rho.*Vst.^2)); % Speed of tube-waves (Bukka et al. (2013), below equation 8);
R1 = R./cosd(theta) + L0./2.*tand(theta); % semi-major axis of the ellipse formed by the intersection of the circular borehole with a tilted fracture plane (Tang and Cheng, 1993, page 173)
Le = pi.*(3.*(R1+R)-sqrt((3.*R1+R).*(R1+3.*R))); % Circumference of this ellipse (Approximation of Ramanujan: https://en.wikipedia.org/wiki/Ellipse) 
Rbar = Le/(2*pi); % radius of equivalent circle (Tang and Cheng, 1993, equation 20)

% Open parallel pool
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool(CoreNum);
else
    disp('Parallel pool is already open.');
end

% Matrix of distances between all positions along the borehole
distmat = zeros(length(zvec),length(zvec));
for isrc=1:length(zvec)
    for irec=1:length(zvec)
        distmat(irec,isrc) = abs(zvec(irec)-zvec(isrc))./cT;
    end
end

% Options for nonlinear optimization
opt = optimset('MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6);

% Find the depth elements at which there are fractures
fracel = zeros(size(zfrac));
for ifrac=1:length(zfrac)
    [~,fracel(ifrac)] = min(abs(zvec-zfrac(ifrac)));
end

% Initialize total pressure field on sparse frequency vector
p_sparse = zeros(length(zvec),length(freqvec));

% Initialize P-wave pressure field on sparse frequency vector
pinc_sparse = zeros(length(zvec),length(freqvec));

% Parallel for-loop over frequencies
parfor ifreq=1:length(freqvec)
    freq = freqvec(ifreq);
    omega = 2*pi*freq; % angular frequency [1/s]

    if doCG==1
        kr = zeros(size(zfrac));
        for ifrac=1:length(zfrac)
            temp = find_kr_CG(L0(ifrac),Kf,rhof,nu,omega); % radial wavenumber, equation 4 of Bakku et al. 2013
            kr(ifrac) = complex(abs(real(temp)),abs(imag(temp)));
        end
    else
        kr = zeros(size(zfrac));
        for ifrac=1:length(zfrac)
            temp = fminsearch(@(x) dispersion_relation(x,omega,L0(ifrac),alpha_f,nu),[0.0,0.0],opt); % radial wavenumber, equation 4 of Bakku et al. 2013
            if temp(1)<0
                temp(1) = abs(temp(1));
                fprintf('Warning: real part of kr is negative. Took absolute value.\n')
            end
            kr(ifrac) = complex(temp(1),temp(2));
        end
    end
    zeta = kr.*alpha_f.*inv_alpha_eff; % effective radial wavenumber (equation A3)
    H0 = besselh(0,1,zeta.*Rbar); % changed for tilted fractures
    H1 = besselh(1,1,zeta.*Rbar); % changed for tilted fractures

    % interface compliance (equation 4 in Minato and Ghose, JGR)
    eta = -(2.*zeta.*Rbar.*L0.*H1)./(R.^2*kr.^2*alpha_f.^2*rhof.*H0); % changed for tilted fractures

    % Tube wave pressure (equation 8, without sigma0 which is eliminated when calculating gammag in equation 10)
    pt = (Le./(2.*pi.*R)).*1i.*omega.*cT.*rhof.*Z.*(1./inv_alpha_eff).*H1./(kr.*alpha_f.*R.*H0); % changed for tilted fractures

    % normally-incident plane P-wave (equation 9, without sigma0 which is eliminated when calculating gammag in equation 10)
    pinc = rhof.*cT.^2./(rho.*Vs.^2).*((1-2.*Vs.^2./Vp.^2)./(1-cT.^2./Vp.^2));

    % reduce the amplitude of the incoming P-wave according to the transmission coefficient across the fracture
    pinc_z = ones(size(zvec))*pinc;
    for ifrac=1:length(fracel)
        Tp = get_TP(mu,Ks,theta(ifrac),rho,Z(ifrac),freq);
        pinc_z(fracel(ifrac):end) = pinc_z(fracel(ifrac)-1)*Tp;
    end

    % reduce the amplitude of the incoming p-wave due to geometrical spreading
    pinc_z = pinc_z./(zvec.^gsexp);

    % Amplitude ratio of incident P wave and the generated tube wave (equation 10)
    gammag = pt./pinc_z(fracel);

    % Tube-wave scattering potential (equation 20)
    phis = zeros(size(zvec));
    % Put a fracture at a certain depth
    phis(fracel) = 1i.*omega.*eta;

    % Tube-wave generation potential (equation 22)
    phig = zeros(size(zvec));
    % Put a fracture at a certain depth
    phig(fracel) = 2/(rhof*cT)*gammag;

    % Create the matrix of Greens-functions between all receiver positions (equation C1)
    G = rhof.*cT./2.*exp(-1i.*omega.*distmat);
    
    % Create matrix M (equation C4)
    phismat = repmat(phis,[length(zvec),1]);
    M = phismat.*G; % There is no factor dz, because the amplitude of the Dirac Delta function is 1 and not 1/dz. 
    
    % Create matrix K (equation C5)
    dphimat = repmat(phig-phis,[length(zvec),1]);
    K = dphimat.*G; % There is no factor dz, because the amplitude of the Dirac Delta function is 1 and not 1/dz. 
    
    % Vector of incoming plane P-waves
    pincvec = (pinc_z.').*exp(-1i.*omega.*abs(zvec.'-zsrc)./Vp);
   
    % Just the incoming P-wave
    pinc_sparse(:,ifreq) = pincvec;
 
    % Compute the total tube wavefield (C7)
    p_sparse(:,ifreq) = (eye(size(M))-M)\(eye(size(K))+K)*pincvec;
end

% Interpolate on dense frequency vector
[zgrid,freqgrid] = ndgrid(zvec,freqvec);
[zlingrid,freqlingrid] = ndgrid(zvec,freqvec_out(2:end));
p = zeros(length(zvec),length(freqvec_out));
ponly = zeros(length(zvec),length(freqvec_out));
p(:,2:end) = interpn(zgrid,freqgrid,p_sparse,zlingrid,freqlingrid,'spline');
ponly(:,2:end) = interpn(zgrid,freqgrid,pinc_sparse,zlingrid,freqlingrid,'spline');

% Assign a wavelet
if doricker==1 % Create a Ricker wavelet
    fwave = zeros(1,length(freqvec_out));
    fwave(1,:) = -(freqvec_out/fc).^2.*exp(-(freqvec_out/fc).^2);
    waveletmat = repmat(fwave,[length(zvec) 1]);
else % Load a wavelet from a file
    load wavelet10ms.mat seis time
    nfreq_wave = 2*size(p,2);
    dfreq_wave = 1/(nfreq_wave*(time(2)-time(1)));
    freqvec_wave = linspace(-nfreq_wave/2,nfreq_wave/2-1,nfreq_wave)*dfreq_wave;
    fwave = fftshift(fft(seis,2*size(p,2))*(time(2)-time(1)));
    fwave_int = interp1(freqvec_wave,fwave,freqvec_out);
    waveletmat = repmat(fwave_int,[length(zvec) 1]);
end
p = p.*waveletmat;
ponly = ponly.*waveletmat;

% Extend the frequency range by zero-padding to obtain a smaller dt. 
extfac = 8; 
dfreq = freqvec_out(2)-freqvec_out(1);
freqvec_out = linspace(0,extfac*(freqvec_out(end)+dfreq)-dfreq,extfac*length(freqvec_out));
padmat = zeros(length(zvec),length(freqvec_out)-size(p,2));
p = [p,padmat];
ponly = [ponly,padmat];

% FFT to the time domain
dfreq = freqvec_out(2)-freqvec_out(1);
temp = conj(flipdim(p,2));
temp = cat(2,zeros(length(zvec),1),temp);
datamat = cat(2,temp,p(:,2:end));
temp_ponly = conj(flipdim(ponly,2));
temp_ponly = cat(2,zeros(length(zvec),1),temp_ponly);
datamat_ponly = cat(2,temp_ponly,ponly(:,2:end));
nfreq = 2*length(freqvec_out);
xtdata = nfreq*fftshift(ifft(fftshift(datamat,2),[],2)*dfreq,2);
xtdata_ponly = nfreq*fftshift(ifft(fftshift(datamat_ponly,2),[],2)*dfreq,2);
dt = 1/(nfreq*dfreq);
tvec = linspace(-nfreq/2,nfreq/2-1,nfreq)*dt;

% As absolute amplitudes are not of interest, we normalize the data
normfac = max(abs(xtdata(:))); % normalization factor
xtdata = xtdata./normfac; % normalize the data
xtdata_ponly = xtdata_ponly./normfac; % normalize the data

% Extract a dataset that contains only the tube waves
xtdata_tonly = xtdata-xtdata_ponly;

% Apply attenuation to tube waves as a function of time. 
% Move-out correction
tmoveout = (zvec-zsrc)/Vp;
tmoveoutel = round(tmoveout/dt);
for iz=1:length(zvec)
    xtdata_tonly(iz,:)=circshift(xtdata_tonly(iz,:),[0,-tmoveoutel(iz)]);
end
% Apply attenuation
tubeattvec = exp(-gsexptube*tvec);
tubeattvec(1:length(tvec)/2) = 1;
tubeattvec=circshift(tubeattvec,[0,round(tshift/dt)]);
tubeattvec(1:length(tvec)/2) = 1;
tubeattmat = repmat(tubeattvec,[length(zvec),1]);
xtdata_tonly = xtdata_tonly.*tubeattmat;
% Undo move-out correction
for iz=1:length(zvec)
    xtdata_tonly(iz,:)=circshift(xtdata_tonly(iz,:),[0,tmoveoutel(iz)]);
end

% Plotting
figure;
imagesc(tvec,zvec,xtdata_tonly+xtdata_ponly)
colormap('gray')
colorbar
cax = caxis;
xlim([0 0.015])
xlabel('Time [s]','Fontsize',14)
ylabel('Depth [m]','Fontsize',14)
title('Complete wavefield','Fontsize',14)
set(gca,'Fontsize',14)

figure;
imagesc(tvec,zvec,xtdata_tonly)
colormap('gray')
caxis(cax)
colorbar
xlim([0 0.015])
xlabel('Time [s]','Fontsize',14)
ylabel('Depth [m]','Fontsize',14)
title('Tube-waves only','Fontsize',14)
set(gca,'Fontsize',14)

figure;
imagesc(tvec,zvec,xtdata_ponly)
colormap('gray')
caxis(cax)
colorbar
xlim([0 0.015])
xlabel('Time [s]','Fontsize',14)
ylabel('Depth [m]','Fontsize',14)
title('P-wave only','Fontsize',14)
set(gca,'Fontsize',14)

% Close current parallel pool
if doclosepool==1
    delete(gcp('nocreate'))
end

toc
