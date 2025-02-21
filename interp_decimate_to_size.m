%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Farrukh Javed Chaudhary                           %
%               Centre for Advanced Studies and Engineering               %
%                         Islamabad, Pakistan                             %
%                      farrukh_javed56@yahoo.com                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%              interp or decimate to a given size(Function)               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sig_out = interp_decimate_to_size(sig_in size)
%%%%%%%%%%%%%% Interpolat to original size
mismatch = rem(length(sig_in),size);
if mismatch~=0 && mismatch > 1
    sig_out = [interp1([1:1:length(sig_in)],sig_in,[1:1:size])]';
elseif mismatch~=0 && mismatch < 1
    sig_out = decimate(sig_in, length(sig_in)/size);
else
    sig_out = sig_in;
end