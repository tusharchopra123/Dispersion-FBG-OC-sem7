% script for generating structure for optical fiber paraeters
% output structur contains all the necessary performance related parametres
% for the project needed to estimate lenght of fiber
% input parameters:   Dt      Max dispersion (ps/nm.km) 
%                     So      Zero dispersoin slope (ps/nm2.km)
%                     alpha   attenuatoin/km (db/km)
%                     dalpha  Max attenuation difference (db/km) 
%                     neff    effective group index for refraction
% output structure:   Dt      Max dispersion (ps/nm.km)  
%                     So      Zero dispersoin slope (ps/nm2.km)
%                     alpha   attenuatoin/km (db/km)
%                     dalpha  Max attenuation difference (db/km) 
%                     neff    effective group index for refraction
% all parameters at 1550nm wavelength
function optic_fiber = optic_fiber_parameter(varargin)
%% seting default parameters for SMF-28
if length(varargin)>5
    error('optic_fiber_parameter: too many input arguements')
end

defaults = {18.0,0.092,0.18,0.02,1.4682};
defaults(1:nargin) = varargin;

optic_fiber = struct(   'Dt',       defaults{1},...
                        'So',       defaults{2},...
                        'alpha',    defaults{3},...
                        'dalpha',   defaults{4},...
                        'neff',     defaults{5} );
end
                




