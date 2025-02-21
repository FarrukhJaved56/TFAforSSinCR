%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Farrukh Javed Chaudhary                           %
%               Centre for Advanced Studies and Engineering               %
%                         Islamabad, Pakistan                             %
%                      farrukh_javed56@yahoo.com                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           SS for OFDM                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The simulation is for 1000 frames.
% We cannot go for spectrum holes because the zeros cannot be used for
% transmission. Otherwise PU will also consider them as data.
% We can use underlay scanario where the oppertunity will be defined as txmn
% below interference cap.
% The buffer size would be the number of frames which are observed to
% define an interference cap.
% The WLAN parameters are assumed to be known and hence we can know which
% carriers are permanently off and thus have noise only. We calc a cap and
% based on these carriers.
% use it on the data block to see which carriers are off.
% We are not considering pilot and trq seq for calc error
%%%%%%%%%%%%%%%%%%%%%%%%%% Load input signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
load ('W-LAN txmn')
load ('WLAN_parameters')
load ('nz_per_frame')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make small vectors %%%%%%%%%%%%%%%%%%%%%
txmn = wl_txmn;
fuc = 100; % Frames under consideration for one cycle
PFD_avg = zeros(1,1000);PLD_avg = zeros(1,1000);
for m = 1:10
    wl_txmn = txmn((m-1)*(txd_frame*fuc)+1:m*txd_frame*fuc);
    frames = floor(length(wl_txmn) / txd_frame);

    col = data_cols+trg_seq_rep;% cols per frame
    zr = 11; % zero rows per frame (54:64). These will be the noise only rows as they have noise and trg seq

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Demodulate from BT freq %%%%%%%%%%%%%%%%%%%%%
    Fs_c = 2 * 128e6 ; % Total band 2 * 128 MHz
    carrier_conj = [repmat(conj(fskmod(39,128,1e6,txd_frame, Fs_c)),1,frames)]'; 
    wl_txmn = wl_txmn .* carrier_conj(1:length(wl_txmn));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        chan = ricianchan(1/(1e6),100,1);
        wl_med = filter(chan, wl_txmn);

        %%%%%%%%%%%%%%%%%%%%%%%%%% Changing buffer size %%%%%%%%%%%%%%%%%%%%%%%%%%%
        k = 5;  % Buffer size. These many frames are analysed for calc cap
            l = 1*k;% k is buffer size in number of frames
            noise = reshape(wl_med,CP_frame_rows,[]);
            noise = noise(70:80,:);
            noise = reshape(noise,1,[]);
            noise_adj = noise(1:l*zr*CP_frame_rows);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            B0 = (1/(sqrt(OFDM_frame_rows)*sqrt(OFDM_frame_rows/(OFDM_pilot_rows -1))))* specgram(noise_adj,24,1,[],0);
            PSD0 = B0 .* conj(B0) / OFDM_frame_rows; % Reduce amp of FFT

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT Sig and remove CP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            B1 = (1/(sqrt(OFDM_frame_rows)*sqrt(OFDM_frame_rows/(OFDM_pilot_rows -1))))* specgram(wl_med,24,1,[],0);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Delete padded zeros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            B1_data = reshape(B1,txd_frame,[]); % We have frames col wise
            B1_data(1:16*(data_cols+trg_seq_rep),:) = [];        % Removing CP. Each frame is now 53*24
            B1_data(53*(data_cols+trg_seq_rep)+1:end,:) = [];        % Removing padded zeros 54:64 in each col. Each frame is now 53*24
            B1_data = reshape (B1_data, (data_cols+trg_seq_rep), []); % Frames are row wise now
            trg_rxd = B1_data(1:trg_seq_rep,:); % A trg matrix of 52 * 4 with frames row wise
            B1_data = [B1_data(trg_seq_rep+1:end,:)]'; % Delete trg seq and transpose
            B1_data = reshape (B1_data, OFDM_pilot_rows, []); % Frames are row wise now
            B1_data([6 20 27 34 48],:) = [];        % Removing pilots. Each frame is now 48*24
            B1_data = reshape (B1_data, [], data_cols); % Frames are row wise now

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc PSD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PSD1 = B1_data .* conj(B1_data) / OFDM_frame_rows; % Consider B1 without the trg seq
            PSD1_NB = PSD1(1:l*(OFDM_pilot_rows - pilot_rep),:); % These many cols of each frames considered for calc ST

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc ST FFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            m_s0 = mean(reshape(PSD0,1,[])); % Mean of noise
            sd_s0 = std(reshape(PSD0,1,[]));

            m_s1 = mean(reshape(PSD1_NB,1,[])); % Mean of noise
            sd_s1 = std(reshape(PSD1_NB,1,[]));

            %%%  Noise threshold (ST), SEF, PFD,PLD theoratically calc, ref paper [14]
            ST_th = ((m_s0*sd_s1) + (m_s1*sd_s0)) ./ (sd_s0+sd_s1);
            temp_PSD1 = PSD1; PFD_actual = []; PLD_actual = [];
            ST_min = min(m_s0 , m_s1) - min(sd_s0,sd_s1);
            ST_max = max(m_s0 , m_s1) + max(sd_s0,sd_s1);

            for i = 1:1000
                PSD1 = temp_PSD1;
                ST = ST_min + ((i-1)  * (ST_max - ST_min)/1000);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Compare with ST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            white_space_1 = double(PSD1 >= ST);% The values of PSD less than ST will now become zeros

            %%%%%%%%%%%%%%%%%%%%%% Calculate actual PFD and PLD %%%%%%%%%%%%%%%%%%%%%%%
            %%% Creat actual occupation of each frame
            % CP is only for synch and will be disregarded in calc of actual PFD and PLD

            actual_frame = [ones(nz_per_frame,1);zeros(data_frame_size- nz_per_frame,1)];
            actual_frame = [reshape(actual_frame, data_cols,[])]';
            cho = repmat(actual_frame,frames,1);

            %%% Actual PFD and PLD
            PFD_actual(i) = sum(sum(double((cho == 1) & (cho ~= white_space_1)))) / (length(find(cho)));
            PLD_actual(i) = sum(sum(double((cho == 0) & (cho ~= white_space_1)))) / ((size(cho,1) * size(cho,2))- length(find(cho)));
            end
    PFD_avg = PFD_actual + PFD_avg;
    PLD_avg = PLD_actual + PLD_avg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PFD = [PFD_avg./m];
PLD = [PLD_avg./m];

xlswrite('OFDM1_RC',PFD',3);
xlswrite('OFDM1_RC',PLD',4);