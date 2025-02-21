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
PFD_avg = zeros(10);PLD_avg = zeros(10);
for m = 1:10
    wl_txmn = txmn((m-1)*(txd_frame*fuc)+1:m*txd_frame*fuc);
    frames = floor(length(wl_txmn) / txd_frame);

    col = data_cols+trg_seq_rep;% cols per frame
    zr = 11; % zero rows per frame (54:64). These will be the noise only rows as they have noise and trg seq

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Demodulate from BT freq %%%%%%%%%%%%%%%%%%%%%
    Fs_c = 2 * 128e6 ; % Total band 2 * 128 MHz
    carrier_conj = [repmat(conj(fskmod(39,128,1e6,txd_frame, Fs_c)),1,frames)]'; 
    wl_txmn = wl_txmn .* carrier_conj(1:length(wl_txmn));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Demux %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wl_txmn = [reshape(wl_txmn,(data_cols+trg_seq_rep),[])]';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for SNR_factor = 1:10 % The changing noise
        SNR = (SNR_factor * 5) - 45; % SNR changes from -40 to 5

        wl_med = awgn(wl_txmn,SNR,'measured'); % Add noise as per SNR given


        %%%%%%%%%%%%%%%%%%%%%%%%%% Changing buffer size %%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:10  % Buffer size. These many frames are analysed for calc cap 
            l = 1*k;% k is buffer size in number of frames
            noise = reshape(wl_med,CP_frame_rows,[]);
            noise = noise(70:80,:);
            noise = reshape(noise,1,[]);
            noise_adj = noise(1:l*zr*CP_frame_rows);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            B0 = (1/(sqrt(OFDM_frame_rows)*sqrt(OFDM_frame_rows/(OFDM_pilot_rows -1))))* fft(noise_adj,CP_frame_rows);
            PSD0 = B0 .* conj(B0) / OFDM_frame_rows; % Reduce amp of FFT

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT Sig and remove CP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for z = 1:frames
                B1((64*(z-1))+(1:64),:) = (1/(sqrt(OFDM_frame_rows)*sqrt(OFDM_frame_rows/(OFDM_pilot_rows -1))))* fft(wl_med((80*(z-1))+(17:80),:),[],2);
                % FFT is for each row to see which rows are vaccant. Each row is multiplied
                % with one sub_carrier, so that sub_carriers is free for that frame,
                % Provided it is below interference cap we may transmit
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Delete padded zeros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            B1_data = reshape(B1,OFDM_frame_rows,[]); % We have first col each frame then sec . . .  
            B1_data(54:64,:) = [];        % Removing padded zeros 54:64 in each col. Each frame is now 53*24
            B1_data([6 20 27 34 48],:) = [];        % Removing pilots. Each frame is now 48*24
            B1_data = reshape (B1_data, frames*(OFDM_pilot_rows - pilot_rep), []); % Frames are row wise now
            trg_rxd = B1_data(:,1:trg_seq_rep); % A trg matrix of 52 * 4 with frames row wise
            B1_data = B1_data(:,trg_seq_rep+1:end); % Delete trg seq

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Equalize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        trg_txd = repmat(trg_seq,frames,trg_seq_rep);

            %%% Calc equalisation factor
    %        u = trg_rxd ./ trg_txd;
    %        B1_eq = (repmat((conj(mean(u,2)) ./ mean(u .^ 2 , 2)) , 1,data_cols)) .* B1_data ;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc PSD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PSD1 = B1_data .* conj(B1_data) / (OFDM_frame_rows); % Consider B1 without the trg seq

            PSD1_NB = PSD1(1:l*(OFDM_pilot_rows - trg_seq_rep),:); % These many cols of each frames considered for calc ST

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc ST FFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            m_s0 = mean(reshape(PSD0,1,[])); % Mean of noise
            sd_s0 = std(reshape(PSD0,1,[]));

            m_s1 = mean(reshape(PSD1_NB,1,[])); % Mean of noise
            sd_s1 = std(reshape(PSD1_NB,1,[]));

            %%%  Noise threshold (ST), SEF, PFD,PLD theoratically calc, ref paper [14]
            ST_th = ((m_s0*sd_s1) + (m_s1*sd_s0)) ./ (sd_s0+sd_s1);
            PFD_th(k,SNR_factor) = qfunc((ST_th-m_s0)/sd_s0);
            PLD_th(k,SNR_factor) = qfunc(-(ST_th-m_s1)/sd_s1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Compare with ST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            white_space_1 = double(PSD1 >= ST_th);% The values of PSD less than ST will now become zeros

            %%%%%%%%%%%%%%%%%%%%%% Calculate actual PFD and PLD %%%%%%%%%%%%%%%%%%%%%%%
            %%% Creat actual occupation of each frame
            % CP is only for synch and will be disregarded in calc of actual PFD and PLD

            actual_frame = [ones(nz_per_frame,1);zeros(data_frame_size- nz_per_frame,1)];
            actual_frame = [reshape(actual_frame, data_cols,[])]';
            cho = repmat(actual_frame,frames,1);

            %%% Actual PFD and PLD
            PFD_actual(k,SNR_factor) = sum(sum(double((cho == 1) & (cho ~= white_space_1)))) / (length(find(cho)));
            PLD_actual(k,SNR_factor) = sum(sum(double((cho == 0) & (cho ~= white_space_1)))) / ((size(cho,1) * size(cho,2))- length(find(cho)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    PFD_avg = PFD_actual + PFD_avg;
    PLD_avg = PLD_actual + PLD_avg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PFD = [PFD_avg./m;PFD_th];
PLD = [PLD_avg./m;PLD_th];

xlswrite('OFDM1',PFD,1);
xlswrite('OFDM1',PLD,2);