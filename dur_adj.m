%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Farrukh Javed Chaudhary                           %
%               Centre for Advanced Studies and Engineering               %
%                         Islamabad, Pakistan                             %
%                      farrukh_javed56@yahoo.com                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%          Adjust duration of signals to make them equal (Function)       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sig_out = dur_adj(sig_in,Fs,offset,uc,terminal);

offset_z = zeros(offset * Fs,1); % Zero padding for offset

data = sig_in(1: (uc * Fs)); % Data values

terminal_z = zeros((terminal*Fs),1); % Terminal zeros to equate signal

sig_out = [offset_z;data;terminal_z]; % Signal for the given duration
