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
    SNR_factor = 9; % The changing noise
        SNR = (SNR_factor * 5) - 45; % SNR changes from -40 to 5

        wl_med = awgn(wl_txmn,SNR,'measured'); % Add noise as per SNR given
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Corr of sig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for w = 1 : frames * CP_frame_rows
            AC11((w-1)*(2*col-1)+1:w*(2*col-1)) = xcorr(wl_med(((w-1)*col)+1:w*col));
            AC1(w) = max(AC11((w-1)*(2*col-1)+1:w*(2*col-1)));
        end

        k = 5; 
            l = 1*k;% k is buffer size in number of frames
            noise = reshape(wl_med,CP_frame_rows,[]);
            noise = noise(70:80,:);
            noise = reshape(noise,1,[]);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Auto Correlation noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            AC0 = xcorr(noise(1:l*zr*CP_frame_rows)); % Consider L frames for calc threshold

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Auto Correlation Sig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            AC1_k = AC1(1:l*CP_frame_rows); % Consider L frames for calc threshold

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc ST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            m_s0 = mean(reshape(abs(AC0),1,[])); % Mean of noise
            sd_s0 = std(reshape(abs(AC0),1,[]));

            m_s1 = mean(reshape(abs(AC1_k),1,[])); % Mean of noise
            sd_s1 = std(reshape(abs(AC1_k),1,[]));

            %%%  Noise threshold (ST), SEF, PFD,PLD theoratically calc, ref paper [14]
            ST_th = ((m_s0*sd_s1) + (m_s1*sd_s0)) ./ (sd_s0+sd_s1);
            temp_AC1 = AC1; PFD_actual = []; PLD_actual = [];
            ST_min = min(m_s0 , m_s1) - min(sd_s0,sd_s1);
            ST_max = max(m_s0 , m_s1) + max(sd_s0,sd_s1);

            for i = 1:1000
                AC1 = temp_AC1;
                ST = ST_min + ((i-1)  * (ST_max - ST_min)/1000);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Compare with ST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            white_space = double(AC1 >= ST);% The values of PSD less than ST will now become zeros
            white_space = reshape(white_space, CP_frame_rows,[]);
            white_space(1: (CP_frame_rows-OFDM_frame_rows),:) = [];
            white_space(54:64,:) = [];
            white_space([6 20 27 34 48],:) = [];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Remove padding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cho = [ones(ceil(nz_per_frame/data_cols),1);zeros(48 - ceil(nz_per_frame/data_cols),1)];
            cho = repmat(cho,1,frames);

            %%% Actual PFD and PLD
            PFD_actual(i) = sum(sum(double((cho == 1) & (cho ~= white_space)))) / (length(find(cho)));
            PLD_actual(i) = sum(sum(double((cho == 0) & (cho ~= white_space)))) / ((size(cho,1) * size(cho,2))- length(find(cho)));
        end
    PFD_avg = PFD_actual + PFD_avg;
    PLD_avg = PLD_actual + PLD_avg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PFD = [PFD_avg./m];
PLD = [PLD_avg./m];

xlswrite('OFDM1_ST',PFD',7);
xlswrite('OFDM1_ST',PLD',8);