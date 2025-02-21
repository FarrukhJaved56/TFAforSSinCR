%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Farrukh Javed Chaudhary                           %
%               Centre for Advanced Studies and Engineering               %
%                         Islamabad, Pakistan                             %
%                      farrukh_javed56@yahoo.com                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%              Encode: Discretise and digitise (Function)                 %
%                       Not considering zeros                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sig_in_bps Fs_factor] = eadc(sig_in , enc_bits)

%%%%%%%%%%%%%% Change signal to binary
sig_in_int = uencode(sig_in, enc_bits); % Encode to integers
%   encode to 8 bit symbols. 2^enc_bits = 256 levels
sig_in_bin = dec2bin (sig_in_int); % change to bniary. 
%   bps = Fs x enc_bits matrix
sig_in_bps = bin2dec(reshape (sig_in_bin, [],1)); 
%   Change each bit into a row and change to double

Fs_factor = enc_bits;