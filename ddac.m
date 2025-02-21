%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Farrukh Javed Chaudhary                           %
%               Centre for Advanced Studies and Engineering               %
%                         Islamabad, Pakistan                             %
%                      farrukh_javed56@yahoo.com                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        Decoder (DAC) (Function)                         %
%                         Not considering zeros                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sig_out = ddac(sig_out_bps, enc_bits, length_sig_in)

%%%%%%%%%%%%%% Change binary to original signal
sig_out_bin = reshape (dec2bin(sig_out_bps), [] ,enc_bits);
sig_out_int = uint16(bin2dec(sig_out_bin));
sig_out = udecode (sig_out_int , enc_bits);

%%%%%%%%%%%%%% Interpolat to original size

mismatch = rem(length(sig_out),length_sig_in);
if mismatch~=0
    sig_out = [interp1([1:1:length(sig_out)],sig_out,[1:1:length_sig_in])]';
    [i,j] = find(isnan(sig_out)); % Index of the NaN elements
    sig_out(i.*j , 1) = 0;
end

