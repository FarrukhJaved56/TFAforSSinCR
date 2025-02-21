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
% We are checking it with FFT because specgram is nothing but STFT
% We have kept the hop time exactly as hop time, hence there are not cross
% terms. Results are to be compared with other TFA techniques.
% The theoratical value of ST gives poor results as the overlap of PDFs is
% not as it should theoratically be. Hence I have used the PDF of max in
% each frame only
% The simulation is for 128 frames

%%%%%%%%%%%%%%%%%%%%%%%%%% Load input signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
load ('bt_txmn')
load ('bluetooth_parameters')
load ('hop_freqs')

PFD = zeros(1,1000); PLD = zeros(1,1000); 
for reps = 1:100
sh = max_bits_per_hop * mod_bits; % samples per hop
huc = 128; % Hops under consideration
bt_txmn = bt_txmn(1:(huc*sh)); % Consider only first huc frames for analysis
hops = floor(length(bt_txmn) / sh);

%%%%%%%%%%%%%%%%%%%%%%%%%% Changing buffer size %%%%%%%%%%%%%%%%%%%%%%%%%%%
factor = 5;

    NB = round(sh / (factor*5)); % Buffer size. These many samples per hop are analysed for decision
    bt_txmn_NB2 = reshape(bt_txmn, [],hops);
    bt_txmn_NB2(NB+1:end,:) = [];
    bt_txmn_NB = reshape(bt_txmn_NB2,[],1); % Only considerin NB elements in each frame

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SNR_factor = 5; % The changing noise
        SNR = (SNR_factor * 5) - 50; % SNR changes from -45 to 50

        bt_med_NB = awgn(bt_txmn_NB,SNR,'measured'); % Add noise as per SNR given

        noise = bt_med_NB - bt_txmn_NB; % Extract noise

        B0 = specgram(noise(1:NB), NB, 1, [],0); % Get FFT of noise
        PSD0 = B0 .* conj(B0) / NB; % Get PSD of noise

        m_s0 = mean(reshape(PSD0,1,[])); % Mean of noise for first frame
        sd_s0 = std(reshape(PSD0,1,[])); % Std dev of noise

        %%%%%%%%%%%%%%%%%%%%%%% Analyse signal with SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%
        B1 = specgram(bt_med_NB, NB, 1, [],0); % fs = 1 so index will be 0 to 1
        PSD1 = B1 .* conj(B1) / NB; % Reduce amp of FFT

        m_s1 = mean(reshape(PSD1(:,1),1,[])); % Mean of sig for first frame
        sd_s1 = std(reshape(PSD1(:,1),1,[])); % Std dev of sig

        %%%  Noise threshold (ST), SEF, PFD,PLD theoratically calc, ref paper [14]
        ST_th = ((m_s0*sd_s1) + (m_s1*sd_s0)) ./ (sd_s0+sd_s1);
        PFD_th(factor,SNR_factor) = qfunc((ST_th-m_s0)/sd_s0);
        PLD_th(factor,SNR_factor) = qfunc(-(ST_th-m_s1)/sd_s1);

%        ST_th = ST_th * abs(SNR);
        temp_PSD1 = PSD1; PFD_actual = []; PLD_actual = [];
        ST_min = min(m_s0 - sd_s0, m_s1 - sd_s1);
        ST_max = max(m_s0 + sd_s0, m_s1 + sd_s1);
        
        for i = 1:1000
            PSD1 = temp_PSD1;
            ST = ST_min + ((i-1) * 10*(ST_max - ST_min)/1000);
         
            gray_space = PSD1 .* (PSD1 > ST);% The values of PSD less than ST will now become zeros
            [r,c,v] = find(gray_space);
            loc = floor(rem(((r / NB) * 128)+64,128)); % channels judged as occupied
            if length(loc) ~= 0
                if loc(1) == hop1(c(1))
                    found = 1;
                else found = 0;
                end
                for k = 2: length(loc)
                    if (loc(k) == hop1(c(k))) && (((loc(k)~=loc(k-1)) && (c(k) == c(k-1))) || (c(k) ~= c(k-1)))
                        found = found +1;
                    end
                end
            PFD_actual(i) = 100 * sum(double(loc ~= hop1(c))) / length(loc); % Incorrectly judged
            PLD_actual(i) = ceil(100 * (1 - (found/huc))); % 1 - probability of detection
            else 
            PFD_actual(i) = 0; % Incorrectly judged
            PLD_actual(i) = 100; % 1 - probability of detection
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PFD = PFD + PFD_actual;
    PLD = PLD + PLD_actual;
end

PFD = PFD / reps;
PLD = PLD / reps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlswrite('BT1_ST',PFD',3);
xlswrite('BT1_ST',PLD',4);