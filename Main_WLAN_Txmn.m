%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Farrukh Javed Chaudhary                           %
%               Centre for Advanced Studies and Engineering               %
%                         Islamabad, Pakistan                             %
%                      farrukh_javed56@yahoo.com                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Comparison of FFT and TFA for spectrum sensing in a CR            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Load input signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
load('Input3')
wl_in = signal_in';
wl_Fs_init = Fs_init;

%%%%%%%%%%%%%%%%%%%%%%%%% Adjust sampling rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate or decimate to adjust sampling rate

wl_Fs_desired = 75e3; % We are using a max bit rate of 375 kbps so preferrably
% use a sampling rate that can divide 375 kbps. 375 kbps is also the
% maximum possible bit rate for input signal

[wl_sig_in wl_Fs wl_duration] = Fs_adj(wl_in, wl_Fs_init, wl_Fs_desired);

length_wl_sig_in=length(wl_sig_in);
%%%%%%%%%%%%%%%%%%%%%%%%%% Encode and digitise %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analog to digital conversion and encoding

wl_enc_bits = 8;
[wl_in_bps wl_Fs_factor] = eadc(wl_sig_in , wl_enc_bits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Modulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 4;
% n depends no the number of cdma users.i.e. the number of codes required.
% n is also the spreading factor, Ts = n * Tc where Ts = 1 / (Fs * enc_bits * mod_bits)
% Every sample will be spreaded by n element "bipolar" code
[wl_txmn PN_seq_seed orth_code wl_Fs] = wl_txr(wl_in_bps,wl_Fs,wl_Fs_factor,n);
%wl_Fs_factor = 1;
%wl_Fs is now the new WLAN Fs (1920(txd_frame) * 12.5e3 (frame/ sec)) and factor is now 1 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Maedium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Medium Effects
% Add noise frames in end and beginning to see that rxr works
% We are assuming that data is trnamitted at multiples of frame time only
% even if there are pauses in between. This means that noise will also be in
% frames of 625*mod_bits / 1920 samples. This is to avoid long calc for each row/sample

path_loss = 20; SNR = 40; % In dbs
wl_med = medium_proj(wl_txmn, path_loss,SNR);

save('W-LAN Txmn','wl_txmn','wl_Fs','wl_Fs_factor');

%%%%%%%%%%%%%%%%%%%%%%%%%%%% De-modulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wl_out_bps wl_Fs_factor noise_frame_index]= wl_rxr(wl_med,wl_Fs,wl_Fs_factor,PN_seq_seed,orth_code,n);

%%%%%%%%%%%%%%%%%%%%%%%%%% Digital to Analog %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wl_sig_out = ddac(wl_out_bps, wl_enc_bits, length_wl_sig_in);

%%%%%%%%%%%%%%%%%%%%%%% Analysis of received signal %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% BER Computation
[wl_no_of_errors,wl_BER] = biterr(wl_in_bps(1:min((length(wl_out_bps)),(length(wl_in_bps)))),wl_out_bps(1:min((length(wl_out_bps)),(length(wl_in_bps)))))
wl_error = mean (wl_sig_out(1:min((length(wl_sig_out)),(length(wl_sig_in)))) - wl_sig_in(1:min((length(wl_sig_out)),(length(wl_sig_in))))) 

figure(1)
subplot(3,2,5);plot(wl_sig_in(1:1000))
title('Input Signal')

subplot(3,2,6);plot(wl_sig_out(1:1000))
title('Output Signal')