% script for pulse dispersion and attenuatoin 
% input parameters: optical laser pulse struct
%                   length of fiber in km
%                   attenuation parameter (alpha) db/km
%                   second order disperson parameter (So) in ps/nm2.km
% output = modified optic pusle struct
% optic struct field:   pulse in time domain (ps): function handle 
%                       optical spectrum: array
%                       wavelengths axis: array
function [out_pulse] = smfloss(in_pulse,L,varargin)
nopin = length(varargin);
if nopin > 2
    error('Too many input arguements in smfloss')
end

defaults = {0.14 + 0.04*rand(1),0.092};
defaults(1:nopin) = varargin;

[alpha,So] = defaults{:};

%% dispersion calculation
%FWHM width of optical spectrum of laser pulse
sigmaw = pulsewidth(10.^(in_pulse.spectrum/10),in_pulse.wavelengths);
%total dispersoin in ps/nm.km
Dt = (So*1550/4)*(1 - (1310/1550)^4);
%pulse broadening calculation in ps
pulse_inc = sigmaw*L*Dt;
%rater of change of group delay w.r.t wavelength (ps/nm)
Dt_per_w = Dt*L;

%% attenuating spetcrum
spectrumf = in_pulse.spectrum - alpha*L;
af = 10^(-alpha*L);

%% generating output pulse
t = in_pulse.t;
%calculating delay for each wavelength w.r.t 1550nm (in ps)
tg = Dt_per_w*(in_pulse.wavelengths - 1550);
%converting to index shift
Ng = round(tg/(t(2)-t(1)));
res_pulse = af*in_pulse.pulse;
%inserting delay
for i = 1:length(Ng)
    res_pulse(i,:) = circshift(res_pulse(i,:),Ng(i));
end

figure
plot(Ng)
title('index change per wavelength')

%% saving output pulse
out_pulse = struct('t',t,'pulse',res_pulse,'spectrum',spectrumf,'wavelengths',in_pulse.wavelengths);

fprintf('rms spectrum width:'); disp(sigmaw);
fprintf('Dispersion parameter:'); disp(Dt);
fprintf('dispersion inc/km:'); disp(pulse_inc/L);
fprintf('dispersion introduce (in ps):'); disp(pulse_inc);
end


