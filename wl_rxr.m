%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Farrukh Javed Chaudhary                           %
%               Centre for Advanced Studies and Engineering               %
%                         Islamabad, Pakistan                             %
%                      farrukh_javed56@yahoo.com                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       WLAN Receiver (Function)                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sig_out_bps Fs_factor noise_frame_index]= wl_rxr(sig_med,Fs,Fs_factor,PN_seq_seed,H,n);

load ('wlan_parameters')

%%%%%%%%%%%%%% Generate and multiply with conjugate carrier freqs
no_of_rxd_frames = ceil(length(sig_med)/txd_frame);
Fs_c = 2 * 128e6 ; % Total band 2 * 128 MHz
carrier_conj = [repmat(conj(fskmod(39,128,1e6,txd_frame, Fs_c)),1,no_of_rxd_frames)]'; 
sig_bb = sig_med .* carrier_conj(1:length(sig_med));

frames = no_of_rxd_frames;
noise_frame_index = [];
%%%%%%%%%%%%%% Demultiplex OFDM frames
sig_rxd = [reshape(sig_bb,(data_cols+trg_seq_rep),[])]';

%%%%%%%%%%%%%% FFT de-modulation and removal of CP

for z = 1:frames
    OFDM_out((64*(z-1))+(1:64),:) = (1/(sqrt(OFDM_frame_rows)*sqrt(OFDM_frame_rows/(OFDM_pilot_rows -1))))* round(fft(sig_rxd((80*(z-1))+(17:80),:),[],2));
% Remove the normalisation factor and round off to remove the error
% The CP is removed by simply taking FFT till index 64. The trg seq has to
% be FFTed otherwise it could be removed here also
end

%%%%%%%%%%%%%% OFDM frame
%%% Remove trg seq
OFDM_out(:,1:trg_seq_rep) = []; % We are not using the seq in this code

%%% Remove pilot seq and the padded zeros for IFFT
pad_zero_index = [[54:64:64*frames] [55:64:64*frames] [56:64:64*frames] [57:64:64*frames] [58:64:64*frames] [59:64:64*frames] [60:64:64*frames] [61:64:64*frames] [62:64:64*frames] [63:64:64*frames] [64:64:64*frames]];
pilot_index = [[6:64:64*frames] [20:64:64*frames] [27:64:64*frames] [34:64:64*frames] [48:64:64*frames]];
delete_index = [pad_zero_index pilot_index];
OFDM_out(delete_index,:) = [];

%%%%%%%%%%%%%% Remove zero padding from each frame
mod_frame = reshape(OFDM_out' , data_frame_size, []); % Each frame is now a col
del_rows = 0;z_rows = sum(mod_frame,2); % zero rows will have zero sum
while z_rows(end,1) == 0;
    z_rows(end,:) = []; % Delete last zero row
    del_rows = del_rows + 1;
end
mod_frame(end-del_rows+1:end,:) = []; % Delete last row of data_frame
[mod_data_frame_row,no_of_frames] = size(mod_frame);

sig_demod_bb = reshape(mod_frame, [],1); % Arrange in a single col

%%%%%%%%%%%%%% Baseband dememodulate DQPSK to binary signal
sig_cdma = dpskdemod(sig_demod_bb , mode);

%%%%%%%%%%%%%% DS - CDMA Demodulate
sig_d_cdma = repmat(H(1:n),length(sig_cdma) / n,1) .* reshape(sig_cdma,length(sig_cdma)/n , n);
sig_out_bps = sig_d_cdma(: , 1);

Fs_factor = Fs_factor / n;
%%%%%%%%%%%%%%%%%%%%%%% Analysis of received signal %%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(3,2,2);specgram (sig_demod_bb , Fs_factor*n ,1,[],0 )
% Demodulated DQPSK sig
title('Demodulated signal')

subplot(3,2,4);specgram (sig_med , txd_frame ,1,[],0 )
title('Received Signal')
