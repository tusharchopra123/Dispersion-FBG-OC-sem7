%pulse compensation function
% input parameters:   optical laser pulse structure
%                     FBG structure:  r       reflection coeff
%                                     tau     time delay 
%                                     w       wavelength array for conputation
% ouput: modified optical pulse structure
function [out_pulse] = compensate(in_pulse,FBG)
if nargin>2
    error('Too many arguement in compesate fxn')
elseif nargin<2
    error('Not enough arguement in compensate fxn')
end


%% generating output spectrum and selecting tau to use
dt = in_pulse.t(2) - in_pulse.t(1);
w = in_pulse.wavelengths;
% delay corresponding to each wavelength
fbg_tau = FBG.tau;
fbg_tau(isnan(fbg_tau)) = 0;
tg = zeros(size(w));
% corner cases if not in FBG wavelength calculation range
edge_db_left = 10*log10(FBG.r(1));
edge_db_right = 10*log10(FBG.r(end));
tg_left = mean(fbg_tau(1:10));
tg_right = mean(fbg_tau(end-10:end));
% initializing output spectrum 
spectrumf = zeros(size(w));
r_to_use = zeros(size(w));
for i = 1:length(w)
    index = find(abs(FBG.w - w(i))<=1e-2);
    if ~isempty(index)
        if(length(index)>1)
            spectrumf(i) = in_pulse.spectrum(i) + 10*log10(FBG.r(index(1)));
            r_to_use(i) = FBG.r(index(1));
            tg(i) = fbg_tau(index(1));
        else
            spectrumf(i) = in_pulse.spectrum(i) + 10*log10(FBG.r(index));
            r_to_use(i) = FBG.r(index);
            tg(i) = fbg_tau(index);
        end
    elseif i < length(w)/2
        spectrumf(i) = in_pulse.spectrum(i) + edge_db_left;
        tg(i) = tg_left;
    else
        spectrumf(i) = in_pulse.spectrum(i) + edge_db_right;
        tg(i) = tg_right;
    end
end   

%% generating ouput pulse
% finding index for 1550nm
index_at_1550 = find(abs(w - 1550)<1e-1);
if length(index_at_1550)>1
    index_at_1550 = index_at_1550(1);
elseif isempty(index_at_1550)
    error('comepnsate: cannot find index for 1550nm in wavelength array')
end
% centralizing delay at 1550nm
tg = tg - tg(index_at_1550);
% calculating index change 
Ng = round(tg/dt);

figure
plot(w,Ng)
title('index change for each wavelength for FBG')

% initializing resulting pulse
res_pulse = (in_pulse.pulse).*r_to_use';
% modifying delay
for i = 1:length(w)
    res_pulse(i,:) = circshift(res_pulse(i,:),Ng(i));
end


%% generating output structure
out_pulse = struct('t',in_pulse.t,'pulse',res_pulse,'spectrum',spectrumf,'wavelengths',w);

end


