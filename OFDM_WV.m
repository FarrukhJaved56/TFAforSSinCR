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
    col = data_cols+trg_seq_rep;% cols per frame
    zr = 11; % zero rows per frame (54:64). These will be the noise only rows as they have noise and trg seq

    frames = floor(length(wl_txmn) / txd_frame);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Demodulate from BT freq %%%%%%%%%%%%%%%%%%%%%
    Fs_c = 2 * 128e6 ; % Total band 2 * 128 MHz
    carrier_conj = [repmat(conj(fskmod(39,128,1e6,txd_frame, Fs_c)),1,frames)]'; 
    wl_txmn = wl_txmn .* carrier_conj(1:length(wl_txmn));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Remove CP and trg seq %%%%%%%%%%%%%%%%%%%%%
    wl_txmn = reshape(wl_txmn, txd_frame,[]); % Remove CP from each frame
    wl_txmn(1:col*(CP_frame_rows-OFDM_frame_rows),:) = []; % each frame is now 64 * 24

    wl_txmn = reshape(wl_txmn, col,[]); % Frames col wise

    %%%%%%%%%%%%%%%%%%%%%%%%%% Changing buffer size %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 1:10

        NB = k ; % Buffer size. These many frames are analysed for calc cap
        % Changes from 2 to 19 (remove these many rows out of 20 in each frame)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for SNR_factor = 1:10 % The changing noise
            SNR = (SNR_factor * 5) - 45; % SNR changes from -20 to 25

            wl_med = awgn(wl_txmn,SNR,'measured'); % Add noise as per SNR given

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc PSD1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % FFT is for each row to see which rows are vaccant. Each row is multiplied
            % with one sub_carrier, so that sub_carriers is free for that frame,
            % Provided it is below interference cap we may transmit 
            wl_med = reshape(wl_med,col,[]);

            N = col;
            taumax= round(N/2)-1;
            tau2=-taumax:taumax; 
            indices2 = rem(N+tau2,N)+1;

            PSD1 = []; % Reset PSD
            for i = 1:frames*OFDM_frame_rows
                x = wl_med(:,i);
                tfr= zeros (N,col);  
                for ti=1:round(N/2)-1
                     tau=-(ti-1):(ti-1); 
                     indices= rem(N+tau,N)+1;
                     tfr(indices,ti) = x(ti+tau,1) .* conj(x(ti-tau,1));
                end
                for ti= round(N/2): (col) - round(N/2) +1
                     tfr(indices2,ti) = x(ti+tau2,1) .* conj(x(ti-tau2,1));
                end 
                for ti=(col) - round(N/2)+2 : col
                     tau=-(col-ti):(col-ti); 
                     indices= rem(N+tau,N)+1;
                     tfr(indices,ti) = x(ti+tau,1) .* conj(x(ti-tau,1));
                end 
                tfr= fft(tfr);
                    PSD1(:,i) = [sum(tfr.* conj(tfr)/ N,2)];
            end

            %%%%%%%%%%%%%%%%%%%%%% Remove trg seq and pilots and get PSD0,NB %%%%%%%%%%%%%%
            PSD1(1:trg_seq_rep,:) = [];
            PSD1 = reshape(PSD1',OFDM_frame_rows,[]);

            PSD0 = PSD1(54:64,:);
            PSD0 = reshape(PSD0,[],data_cols);
            PSD0 = PSD0(1:NB*11,:);
            
            PSD1(54:64,:) = [];
            PSD1([6 20 27 34 48],:)=[];
            PSD1 = reshape(PSD1,[],data_cols);
            
            PSD1_NB = PSD1(1:NB*(OFDM_pilot_rows-pilot_rep),:);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc ST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PSD1_NB = reshape(PSD1_NB,1,[]);
            PSD1_NB = sort((PSD1_NB)); % Sort is ascending order
           
            PSD1_NB = PSD1_NB((0.8 * NB * data_frame_size):end); % Keep the higher 20 % values

            m_s0 = mean(reshape(PSD0,1,[])); % Mean of noise
            sd_s0 = std(reshape(PSD0,1,[])); % Std dev of noise

            m_s1 = mean(reshape(PSD1_NB,1,[])); % Mean of noise
            sd_s1 = std(reshape(PSD1_NB,1,[])); % Std dev of noise

            %%%  Noise threshold (ST), SEF, PFD,PLD theoratically calc, ref paper [14]
            ST_th = ((m_s0*sd_s1) + (m_s1*sd_s0)) ./ (sd_s0+sd_s1);
            PFD_th(k,SNR_factor) = qfunc((ST_th-m_s0)/sd_s0);
            PLD_th(k,SNR_factor) = qfunc(-(ST_th-m_s1)/sd_s1);

            ST_th = (m_s0+m_s1)/2; % delete
            %%%%%%%%%%%%%%%%%%%%%% Compare with ST %%%%%%%%%%%%%%%%%%%%%%%
            white_space = double(PSD1 >= ST_th);% The values of PSD less than ST will now become zeros

            %%%%%%%%%%%%%%%%%%%%%% Calculate actual PFD and PLD %%%%%%%%%%%%%%%%%%%%%%%
            %%% Creat actual occupation of each frame
            % CP is only for synch and will be disregarded in calc of actual PFD and PLD

            actual_frame = [ones(nz_per_frame,1);zeros(data_frame_size- nz_per_frame,1)];
            actual_frame = [reshape(actual_frame, data_cols,[])]';
            cho = repmat(actual_frame,frames,1);

            %%% Actual PFD and PLD
            PFD_actual(k,SNR_factor) = sum(sum(double((cho == 1) & (cho ~= white_space)))) / (length(find(cho)));
            PLD_actual(k,SNR_factor) = sum(sum(double((cho == 0) & (cho ~= white_space)))) / ((size(cho,1) * size(cho,2))- length(find(cho)));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    PFD_avg = PFD_actual + PFD_avg;
    PLD_avg = PLD_actual + PLD_avg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PFD = [PFD_avg./m;PFD_th];
PLD = [PLD_avg./m;PLD_th];

xlswrite('OFDM1',PFD,5);
xlswrite('OFDM1',PLD,6);