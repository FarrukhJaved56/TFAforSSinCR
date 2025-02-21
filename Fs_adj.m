%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Farrukh Javed Chaudhary                           %
%               Centre for Advanced Studies and Engineering               %
%                         Islamabad, Pakistan                             %
%                      farrukh_javed56@yahoo.com                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%             Adjust sampling rate to desired level (Function)            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sig_out Fs duration] = Fs_adj(signal_in, Fs_init, Fs_desired);
if  Fs_desired >= Fs_init
    factor = round(Fs_desired / Fs_init); % decimation / interpolation factor
    sig_out = interp(signal_in , factor);%May incr no of samples
    Fs = Fs_init * factor; % Changed sampling rate
else
    factor = round(Fs_init / Fs_desired); % decimation / interpolation factor
    sig_out = decimate(signal_in , factor);%May decr no of samples
    Fs = Fs_init * (1/factor); % Changed sampling rate
end

duration = length(sig_out) * (1 / Fs); % each sample comes after 1 /Fs sec