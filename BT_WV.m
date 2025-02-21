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
% The simulation is for 128 frames with Wigner villey
% The buffer size is greatlt reduced though the computation time has
% enhanced
% Upper half of freqs ?????

%%%%%%%%%%%%%%%%%%%%%%%%%% Load input signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
load ('bt_txmn')
load ('bluetooth_parameters')
load ('hop_freqs')

sh = max_bits_per_hop * mod_bits; % samples per hop
huc = 128; % Hops under consideration
bt_txmn = bt_txmn(1:(huc*sh)); % Consider only first 10 frames for analysis
hops = floor(length(bt_txmn) / sh);

N = fhss_carriers_modulo_2;
%%%%%%%%%%%%%%%%%%%%%%%%%% Changing buffer size %%%%%%%%%%%%%%%%%%%%%%%%%%%
for factor = 1:10

    NB = round(sh / (factor*12)); % Buffer size. These many samples per hop are analysed for decision
    bt_txmn_NB2 = reshape(bt_txmn, [],hops);
    bt_txmn_NB2(NB+1:end,:) = [];
    bt_txmn_NB = reshape(bt_txmn_NB2,[],1); % Only considerin NB elements in each frame

    taumax= round(N/2)-1;
    tau2=-taumax:taumax; 
    indices2 = rem(N+tau2,N)+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for SNR_factor = 1:10 % The changing noise

        SNR = (SNR_factor * 10) - 50; % SNR changes from -40 to 50

        bt_med_NB = awgn(bt_txmn_NB,SNR,'measured'); % Add noise as per SNR given

        noise = bt_med_NB(1:NB) - bt_txmn_NB(1:NB);
        x = noise;
        tfr= zeros (N,NB);  
            for ti=1:round(N/2)-1
                 tau=-(ti-1):(ti-1); 
                 indices= rem(N+tau,N)+1;
                 tfr(indices,ti) = x(ti+tau,1) .* conj(x(ti-tau,1));
            end
            for ti= round(N/2): NB - round(N/2) +1
                 tfr(indices2,ti) = x(ti+tau2,1) .* conj(x(ti-tau2,1));
            end 
            for ti=NB - round(N/2)+2:NB
                 tau=-(NB-ti):(NB-ti); 
                 indices= rem(N+tau,N)+1;
                 tfr(indices,ti) = x(ti+tau,1) .* conj(x(ti-tau,1));
            end 
            tfr= fft(tfr);
            PSD0 = sum(tfr.* conj(tfr) /N,2);

        m_s0 = mean(reshape(PSD0,1,[])); % Mean of noise
        sd_s0 = std(reshape(PSD0,1,[])); % Std dev of noise

        %%%%%%%%%%%%%%%%%%%%%%% Analyse signal with SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%
        bt_med_NB = reshape(bt_med_NB,NB,[]);
        for i = 1:huc

            x = bt_med_NB(:,i);
            [m,loc] = max(fft(x).*conj(fft(x))); % Cater for upper half of freqs because they are lost in conj mult
            if (loc/NB) < 0.5
                shift = 1;
            else
                shift = 0;
            end

            tfr= zeros (N,NB);  
            for ti=1:round(N/2)-1
                 tau=-(ti-1):(ti-1); 
                 indices= rem(N+tau,N)+1;
                 tfr(indices,ti) = x(ti+tau,1) .* conj(x(ti-tau,1));
            end
            for ti= round(N/2): NB - round(N/2) +1
                 tfr(indices2,ti) = x(ti+tau2,1) .* conj(x(ti-tau2,1));
            end 
            for ti=NB - round(N/2)+2:NB
                 tau=-(NB-ti):(NB-ti); 
                 indices= rem(N+tau,N)+1;
                 tfr(indices,ti) = x(ti+tau,1) .* conj(x(ti-tau,1));
            end 
            tfr= fft(tfr);
            if shift==0
                PSD1(:,i) = [sum(tfr.* conj(tfr)/ N,2);zeros(N,1)];
            else
                PSD1(:,i) = [zeros(N,1);sum(tfr.* conj(tfr)/ N,2)];
            end
        end

        m_s1 = mean(reshape(PSD1(:,1),1,[])); % Mean of sig for one frame
        sd_s1 = std(reshape(PSD1(:,1),1,[])); % Std dev of sig

        %%%  Noise threshold (ST), SEF, PFD,PLD theoratically calc, ref paper [14]
        ST_th = ((m_s0*sd_s1) + (m_s1*sd_s0)) ./ (sd_s0+sd_s1);
        PFD_th(factor,SNR_factor) = qfunc((ST_th-m_s0)/sd_s0);
        PLD_th(factor,SNR_factor) = qfunc(-(ST_th-m_s1)/sd_s1);
        SEF_th(factor,SNR_factor) = qfunc((m_s1-m_s0)/(sd_s1+sd_s0));


        gray_space = PSD1 .* (PSD1 > ST_th);% The values of PSD less than ST will now become zeros

        [C,I] = max(gray_space); % The indices of values greater than noise

        cho = floor(I'/2)-1; % channels  occupied

        PFD_actual(factor,SNR_factor) = sum(double(cho ~= hop1)) / huc;
        PLD_actual(factor,SNR_factor) = sum(double(cho == 0)) / huc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PFD = [PFD_actual;PFD_th];
PLD = [PLD_actual;PLD_th];


xlswrite('BT1',PFD_actual,5);
xlswrite('BT1',PLD_actual,6);