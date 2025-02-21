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

Fs_c = fhss_carriers_modulo_2 * freqsep ; % Total band 128 MHz

%%%%%%%%%%%%%%%%%%%%%%%%%% Changing buffer size %%%%%%%%%%%%%%%%%%%%%%%%%%%
%for factor = 1:10
factor = 5;
    NB = round(sh / (factor*5)); % Buffer size. These many samples per hop are analysed for decision
    bt_txmn_NB2 = reshape(bt_txmn, [],hops);
    bt_txmn_NB2(NB+1:end,:) = [];
    bt_txmn_NB = reshape(bt_txmn_NB2,[],1); % Only considerin NB elements in each frame
    
    carrier = fskmod([0:1:127],fhss_carriers_modulo_2,freqsep,NB, Fs_c); 
    carrier = reshape(carrier,NB,[]);
    % All freqs generated

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    for SNR_factor = 1:10 % The changing noise
SNR_factor = 5;
        SNR = (SNR_factor * 5) - 50; % SNR changes from -45 to 0

        bt_med_NB = awgn(bt_txmn_NB,SNR,'measured'); % Add noise as per SNR given

        noise = bt_med_NB - bt_txmn_NB; % Extract noise

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Auto Correlation noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AC0 = xcorr(noise); % Consider L frames for calc threshold

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Auto Correlation Sig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for w = 1 : huc
            for c = 1:fhss_carriers_modulo_2
                AC11((w-1)*(2*NB-1)+1:w*(2*NB-1)) = xcorr(bt_med_NB(((w-1)*NB)+1:w*NB), carrier(:,c));
                AC1(c,w) = max(AC11((w-1)*(2*NB-1)+1:w*(2*NB-1)));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc ST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        m_s0 = mean(reshape(abs(AC0),1,[])); % Mean of noise
        sd_s0 = std(reshape(abs(AC0),1,[]));

        m_s1 = mean(reshape(abs(AC1),1,[])); % Mean of noise
        sd_s1 = std(reshape(abs(AC1),1,[]));

        %%%  Noise threshold (ST), SEF, PFD,PLD theoratically calc, ref paper [14]
        ST_th = ((m_s0*sd_s1) + (m_s1*sd_s0)) ./ (sd_s0+sd_s1);
        PFD_th(factor,SNR_factor) = qfunc((ST_th-m_s0)/sd_s0);
        PLD_th(factor,SNR_factor) = qfunc(-(ST_th-m_s1)/sd_s1);

        temp_AC1 = AC1; PFD_actual = []; PLD_actual = [];
%        ST_th = ST_th * abs(SNR);
        ST_min = min(m_s0 - sd_s0, m_s1 - sd_s1);
        ST_max = max(m_s0 + sd_s0, m_s1 + sd_s1);
        
        for i = 1:1000
            AC1 = temp_AC1;
            ST = ST_min + ((i-1) * (ST_max - ST_min)/1000);
            gray_space = abs(AC1) .* (abs(AC1) > ST);% The values of PSD less than ST will now become zeros
            [r,c,v] = find(gray_space);
            loc = r-1; % channels judged as occupied
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
xlswrite('BT1_ST',PFD',7);
xlswrite('BT1_ST',PLD',8);