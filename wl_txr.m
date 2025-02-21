%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Farrukh Javed Chaudhary                           %
%               Centre for Advanced Studies and Engineering               %
%                         Islamabad, Pakistan                             %
%                      farrukh_javed56@yahoo.com                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      WLAN Transmitter (Function)                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modulate an input speech signal on IEEE 802.22 with specs               %
%% Air interface DS - CDMA                                                %
%% Freq range 5.5 GHz (but we are assuming ISM band 2.4 - 2.4835 GHz)     %
%% Tx Power = 25 mW                                                       %
%% It should reduce after spreading                                       %
%% Channel bandwidth 20 MHz                                               %
%% Channels 13                                                            %
%% Access mode FDMA / TDMA                                                %
%% Duplex mode Half Duplex                                                %
%% Modulation OFDM with subcarrier modulation BPSK/QPSK/16QAM/64QAM       %
%% As per paper CCQ - DQPSK modulation                                    %
%% Error correction code Convolutional                                    %
%% chip rate 6/9/12/18/24/36/48/54 Mbps                                   %
%% Number of chips/slot = 52 modulated symbols per OFDM slot (means 52    %
%% sub-carriers)                                                          %
%% Frame duration Packets of several 100 micro-sec                        %
%% slot per frame : Variable                                              %
%% slot duration : 1 OFDM symbol of 3.3 micro-sec + 0.8 micro-sec guard   %
%% Maximum cell radius Max 100 m                                          %
%% Net data rate upto 25 Mbps (at least 2 chips for each bit)             %
%% Evolutionary concepts — IEEE 802.11n                                   %  
%% Comparable systems HiperLAN/2                                          %
%                                                                         %
%%%%%%%%%%%%%% Steps                                                      %
%% Change to sig to bps                                                   %
%% Spread code with quasi-random orthogonal code                          %
%% DQPSK modulate the sig                                                 %
%% Modulate 52 symbols on 52 sub_carriers. Keep in baseband so the        %
%% sub_carriers will be distributed on sides of zero.                     %
%% Pad zeros to complete last frame and last column                       %
%% Add pilot and trg seq                                                  %
%% IFFT                                                                   %
%% Add CP                                                                 %
%% Reverse for receiver                                                   %
%                                                                         %
%%%%%%%%%%%%% Specs we are using                                          %
%% Bit rate = .375 Mbps (The maximum bit rate that can be accommodated)   %
%% Bit rate after encoding by 8 bits = 3 Mbps                             %
%% Chip rate = 12 Mbps (4 bit spreading)                                  %
%% Block rate = 12.5 k blocks/sec                                         %
%% Symbol rate = 12 M symbols / sec (9600 symbols / block)                %
%% Chips per block = 960 (12 M / 12.5k)                                   %
%% Mod bits = 1  (Mod order is 1 and not 2 as it should be)               %
%% Tx bits per block = 960                                                %
%% Tx bits per symbol = 1                                                 %
%% Data cols per frame = 20                                               %
%% Data rows per frame = 960 / 20 = 48                                    %
%% Pilot rows per frame = 5 (Total 48 + 5 = 53)                           %
%% Total rows per frame = 64                                              %
%% CP rows per frame = 16                                                 %
%% Total rows per transmission frame = 80                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sig_txr PN_seq_seed H Fs]  = wl_txr(sig_in_bps,Fs,Fs_factor,n);

load ('wlan_parameters')
%%%%%%%%%%%%%% DS - CDMA modulation
H = hadamard(n); % Each row / coulumn is orthogonal
sig_cdma = reshape((repmat(H(1:n),length(sig_in_bps),1) .* repmat(sig_in_bps,1,n)) , n*length(sig_in_bps),1);
% Select one row/col (one code) and multiply with each element to spread

Fs_factor = Fs_factor * n;
%%%%%%%%%%%%%% Find sub_carriers to be switched on
nz_per_frame = 960 * (Fs * Fs_factor / selected_symbol_rate); 
% How many slot out of 960 to be filled in each frame of 960
save('nz_per_frame','nz_per_frame');
%%%%%%%%%%%%%% Base band modulate binary sig by DQPSK
% mode is 4 hence QPSK
sig_mod_bb = dpskmod(sig_cdma,mode); 

%%%%%%%%%%%%%% OFDM frame
%%% Make data frames

data_frame_nz = reshape(sig_mod_bb,nz_per_frame,[]); % Each col is now a frame with nz elements
no_of_frames = size(data_frame_nz , 2);
data_frame_col = [data_frame_nz; complex(zeros(data_frame_size - nz_per_frame, no_of_frames))];
data_frame_col2 = reshape(data_frame_col,1,[]);
% Each col is now a frame of data_frame_size
data_frame = [reshape(data_frame_col2,20,[])]'; 
%Frames are now arranged in data_cols of size no_of_frames * data_rows
data_rows = data_frame_size / data_cols;
% The rows in the frame will now be data_frame_size / data_cols = 48

%%% Add pilots, zeros for IFFT and trg seq and perform IFFT
%%% Pilot specs
% The pilot is repeated after 5, 13, 6, 6, 13 rows in each frame ???
% The third pilot rep is complex zeros ???
% The last pilot rep is inverted (maybe to indicate end of frame)???
% The number of rows is now (data_frame_size / data_cols) + pilot_rep = 53
PN_seq_seed = randseed; % Generate a valid seed for random sequence
rand('state',PN_seq_seed); % Reset state to initial seed
PN_seq_unipolar = round(rand(no_of_frames , data_cols)); % Generate a PN sequence of 0 & 1
PN_seq = 2*PN_seq_unipolar - 1;

%%% Trg seq specs
% The trg_seq is to be a vector of size 
        % (frame_size / sym_size) + pilot_rep , trg_seq_rep = 53,4
        % we have loaded a col of size 53
trg_seq = [trg_seq;complex(zeros(OFDM_frame_rows - OFDM_pilot_rows,1))];
trg_seq = repmat(trg_seq,1, trg_seq_rep);
% Zeros have been padded in this step to make frame rows as 64. This allows
% IFFT to be performed

%%% Frame operations (Pilot addition, trq seq addition, IFFT, CP addition)
for z = 1:no_of_frames
    OFDM_frame((64*(z-1))+(1:64),:) =[trg_seq [data_frame((48*(z-1))+(1:5),:);PN_seq(z,:);data_frame((48*(z-1))+(6:18),:);PN_seq(z,:);data_frame((48*(z-1))+(19:24),:);complex(zeros(1,data_cols));data_frame((48*(z-1))+(25:30),:);PN_seq(z,:);data_frame((48*(z-1))+(31:43),:);-1*PN_seq(z,:);data_frame((48*(z-1))+(44:48),:);complex(zeros(OFDM_frame_rows - OFDM_pilot_rows,data_cols))]];
    % Making of an OFDM frame

    % ifft and normalise as per demo ?????? OFDM_pilot_rows are 53-1 (1 DC and 52 sub_carriers)
    sig_mod((64*(z-1))+(1:64),:) = sqrt(OFDM_frame_rows)*sqrt(OFDM_frame_rows/(OFDM_pilot_rows -1))*(ifft(OFDM_frame((64*(z-1))+(1:64),:),[],2));

    sig_cp((80*(z-1))+(1:80),:) = [sig_mod((64*(z-1))+(1:16),:);sig_mod((64*(z-1))+(1:64),:)];
    % Add cyclic prefix. Repeat first 16 rows
end

Fs = txd_frame * 12.5e3;
%%%%%%%%%%%%%% Multiplex OFDM frames
sig_mux = reshape(sig_cp',[],1);

%%%%%%%%%%%%%% Generate and multiply with conjugate carrier freqs
Fs_c = 2 * 128e6 ; % Total band 2 * 128 MHz
carrier = [repmat(fskmod(39,128,1e6,txd_frame, Fs_c),1,no_of_frames)]';
sig_txr = sig_mux .* carrier;

%%%%%%%%%%%%%%%%%%%%%%% Analysis of received signal %%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(3,2,1);specgram (sig_mod_bb , Fs_factor  ,1,[],0 );
title('Modulating signal')

subplot(3,2,3);specgram (sig_txr , txd_frame  ,1,[],0 )
title('Transmitted Signal')