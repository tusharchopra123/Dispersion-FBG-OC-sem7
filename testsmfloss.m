%test smfloss 
close all
clear

t = linspace(0,1e3,1e3);
pulsei = @(t)gaussmf(t,[5 500]);
wavelengths = linspace(1550-10,1550+10,1e3);
spectrumi = gaussmf(wavelengths,[0.1 1550]);

%input pulse structure
in_pulse = struct('pulse',pulsei,'spectrum',spectrumi,'wavelengths',wavelengths);

%% optic fiber 
%calling smf loss function for diff lengths
out_pulse1 = smfloss(in_pulse,10);

%ploting
plot(t,pulsei(t));
hold on
plot(t,out_pulse1.pulse(t));
legend('input wave','L = 50km');
xlabel('time in ps')
ylabel('normalized intensity')

% figure
% plot(t,out_pulse1.pulse(t))
% hold on
% FWHM = pulsewidth(out_pulse1.pulse(t),t);
% sigma = FWHM/(2*sqrt(log(2)));
% [a,t0] = max(out_pulse1.pulse(t));
% t0 = t0*(t(2) - t(1));
% pulse1 = @(t) a*gaussmf(t,[sigma,t0]);
% 
% out_pulse1.pulse = pulse1;
% plot(t,out_pulse1.pulse(t))

%% FBG 
[r,tau,w] = FBG();
out_pulse2 = compensate(out_pulse1,50,r,w);

figure
subplot(211)
plot(t,out_pulse1.pulse(t));
hold on
plot(t,out_pulse2.pulse(t));
xlabel('time in ps')
ylabel('normalized intensity')
legend('FBG in','FBG out')

subplot(212)
plot(wavelengths,out_pulse1.spectrum);
hold on
plot(wavelengths,out_pulse2.spectrum);
xlabel('wavelength in nm')
ylabel('Intensity')
legend('FBG in','FBG out')





