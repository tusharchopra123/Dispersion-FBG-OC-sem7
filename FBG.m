%clear;
% r -> reflection spectrum of FBG
% tau -> time delay of FBG
% w -> wavelength on which calculation has been made
%input arguements:  Lg                  length of FBG in cm
%                   neff                effective group RI at 1550nm
%                   del_neff            index change of RI
%                   Chirp_var           d(lamda)/dz in nm/cm
%                   epodize_fxn_choice  choice of epodization function
%                                       1 - uniform
%                                       2 - gaussian
%                                       3 - raised cosine
% output arguement:   FBG structure 
%                     r       reflectence coeff. of each wavelength
%                     tau     time delay of each wavelength
%                     w       wavelength array used for computation 
%                     
function [fbg_spectrum] = FBG(varargin)
%% input grating parameter
% setting up default values
if length(varargin)>5
    error('Too Many Arguements in FBG function')
end

defaults = {0.5,1.4682,0.0004,2.5,1};
defaults(1:nargin) = varargin;

[Lg,neff,del_neff,chirp_var,epodize_fxn_choice] = defaults{:};

%% calculation of chirp
%operating wavelength in nm
wd = 1550;
%wavelength precision for optical spectrum calculation in nm
del_w = 0.01;
%segmentation of FBG
N = 500;
%length array 
z = linspace(0,Lg,N);
dz = z(2) - z(1);
%chirp -> rate of phase change along z direction in FBG
chirp = (-4*pi*neff*(chirp_var*1e-7)/((wd*1e-9)^2))*(z*1e-2);

%% apodization function
switch epodize_fxn_choice
    case 1
        gz = ones(1,N); %uniform chirp
    case 2
        gz = exp(-64*((z - Lg/2)/Lg).^4); %gaussian profile
    case 3
        gz = 0.5*(1 + cos(pi*(z - Lg/2)/Lg)); %raised cosine 
    otherwise
        error('invalid choice value')
end

%% wavelength span for FBG
% 10nm wavelenght span 
%total number of points for calculation
M = round(10/del_w);
%input wavelength
w = linspace(wd-5,wd+5,M)*1e-9;
%bragg wavelengths variation along z direction
wb = linspace(wd - chirp_var*Lg/2,wd + chirp_var*Lg/2,N)*1e-9;

%% initializing reflectence, phase, delay vector
r=zeros(1,length(w));
phi=zeros(1,length(w));

%% transmisson matrix method
dz = dz*1e-2;
for i = 1:M
    %DC coupling coefficient
    sigmaz = (2*pi*neff)*(1./w(i) - 1./wb) + (2*pi*del_neff)./w(i) - chirp;
    %AC coupling coefficinet
    kz = (pi*del_neff)*(gz/w(i));
    %gamma b
    gamma = sqrt(kz.^2 - sigmaz.^2);
    %matrix elements
    f11 = cosh(gamma*dz) - 1j*(sigmaz./gamma).*sinh(gamma*dz);
    f12 = -1j*(kz./gamma).*sinh(gamma*dz);
    f21 = 1j*(kz./gamma).*sinh(gamma*dz);
    f22 = cosh(gamma*dz) + 1j*(sigmaz./gamma).*sinh(gamma*dz);
    
    %initializing transmission matrix
    F = [1,0;0,1];
    %calculation of transmission matrix for entire FBG for w(i)
    for j = 1:N
        F = F*[f11(j), f12(j);f21(j), f22(j)];
    end
    %reflective rate for wavelength w(i)
    r(i) = (abs(F(3)/F(1))).^2;
    %phase of w(i)
    phi(i) = angle(F(3)/F(1));
end

dtg = [phi(1) diff(phi)];
s = std(dtg);
% plot(diff(dtg))
dtg(abs([0 diff(dtg)])>s) = nan;

tau = -((w.^2)/(2*pi*3e-2)).*(dtg/del_w)*1e12;
%conerting back to nm
w = w*1e9;

fbg_spectrum = struct('r',r,'tau',tau,'w',w);

end %function end


%% display graph
% subplot(2,1,1)
% plot(lamda*1e9,r,'b');
% grid on
% xlabel('Wavelength/nm');
% ylabel('Reflection');
% title('Reflection Spectrum of Chirped Grating');
% 
% subplot(2,1,2)
% plot(lamda*1e9,tau,'r');
% grid on
% xlabel('Wavelength/nm');
% ylabel('Time Delay');
% title('Time Delay of Chirped Grating');