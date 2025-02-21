clc; clear all; close all

%% Assumptions
nb = 10; % Buffer size and bit size in no of samples
nc = 10*100; % Total channels
nos = nc*nb; % Total no of samples
Fs = 10; % Sampling rate in samples/sec
t=1/Fs* (0:1:nos-1); % Total simulated time nos/ Fs secs
fc=2* Fs; % Carrier frequency
SNR = -15:1:15; 
y = zeros(length(SNR),nos);n=y;D=zeros(length(SNR),nc);D2=D;

%% Analysis with changing Beta
Bvar = 2; % Beta variations implying changing priors
APeS = zeros(length(SNR),Bvar+1);APeN = APeS;APeP = APeS;
APeO = APeS;APeAN = APeS;APeso = APeS;
for l = 0:Bvar
    %% Generating Test Data R and decision Metric D
    ph1 = 0.3+l*(0.2);ph0 = 1-ph1;% Priors. ph0 should be greater than ph1
    beta = ph0/ph1; % Current Beta

    %% BPSK
    d=2*randint(1,round(ph0*nc))-1; % Data sequence
    dd=repmat(d',1,nb); % replicate each bit nb times
    dw=dd'; dw=dw(:)';
    % Convert dw to a column vector (colum by column) and convert to a row vector
    w=dw.*cos(2*pi*(1/fc)*t(1:length(dw))); % carrier waveform
    s1 = std(w); m1 = mean(w);
    x=[w, zeros(1,nos-length(dw))]; % modulated waveform
    for i = 1:length(SNR)
        y(i,:) = awgn(x,SNR(i),'measured');
        n = y-repmat(x,length(SNR),1);
        s0 = std(n,0,2); m0 = mean(n,2);
        for k =1:nc
            D(:,k) = sum((abs(y(:,((k-1)*nb)+1:k*nb))).^2,2);
        end
    end

%% STs
    muD1 = mean(D(:,1:round(ph0*nc)),2);sD1 = std(D(:,1:round(ph0*nc)),0,2);
    muD0 = mean(D(:,round(ph0*nc)+1:end),2);sD0 = std(D(:,round(ph0*nc)+1:end),0,2);
 
    pe1 = (muD1+muD0)./(sD1.^2 - sD0.^2);K = sD0 ./ sD1;
    pe2 = ((sD1.^2 .* muD0.^2) + (sD0.^2.* muD1.^2) + ...
        (sD0.^2 .* sD1.^2 .* (log((beta^2)*(1./K.^2)))))./(sD1.^2 - sD0.^2);
    STpe = pe1 + sqrt((pe1.^2)+ pe2);

    %% PFD and PLD
    PLD3 = sum(double(D(:,1:round(ph0*nc)) <= repmat(STpe,1,round(ph0*nc))),2)/(ph0*nc);
    PFD3 = sum(double(D(:,(round(ph0*nc))+1:end) >= repmat(STpe,1,nc-(round(ph0*nc)))),2)/(nc-(round(ph0*nc)));
    APeP(:,l+1) = smooth(smooth(smooth((ph0*PFD3)+(ph1*PLD3))));

end
plot(SNR,APeP(:,1),'g',SNR,APeP(:,2),'--r',SNR,APeP(:,3),':k','LineWidth',5);
title('P_e(\lambda_{Pe}) for Changing Priors','fontsize',18,'fontweight','b');
legend( 'P(h_0)=0.7', 'P(h_0)=0.5', 'P(h_0)=0.3');ylim([0 1]);
xlabel('SNR','fontsize',16,'fontweight','b');
ylabel('P_e','fontsize',16,'fontweight','b');

