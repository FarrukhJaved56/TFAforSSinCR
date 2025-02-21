clc; clear all; close all

%% Assumptions
nb = 10; % Buffer size. For Soft threshold nb = 1 beause each sample is tested
nc = 20*1000; nos = nc*nb; % no of channels (multiple of 100) and samples
mu0=0; % Noise is zero mean
s1 = 6;
SNR_db = -15:1:15; SNR = (db2pow(SNR_db))'; % SNR is mu1^2 / s0^2
s0 = sqrt((s1^2)./SNR); % sigma for noise according to SNR
%s0 = ones(size(SNR));
%% Test statistic

ts1 = (s1^2./((2*s0.^2).*(s0.^2+s1^2)));
ts2 = (0.5*log(s0.^2./(s0.^2+s1^2)));


%% Analysis with changing Beta
Bvar = 2; % Beta variations implying changing priors
APeAN = zeros(length(SNR),Bvar+1);APeso = APeAN;APeso2 = APeso;
for l = 0:Bvar
    %% Generating Test Data R and decision Metric D
    ph1 = 0.3+l*(0.2);
    ph0 = 1-ph1;% Priors. ph0 should be greater than ph1
    beta = ph0/ph1; % Current Beta
   
    z = zeros(length(SNR),round(ph0*nos));y = zeros(length(SNR),round(ph1*nos));s=y;
    for i = 1:length(SNR)
        z(i,:) = normrnd(0,s0(i),1,round(ph0*nos));% Generate noise. 
        s(i,:) = normrnd(0,s1,1,round(ph1*nos));% Rxd sig without noise
        y(i,:) = s(i,:)+ (normrnd(mu0,s0(i),1,round(ph1*nos)));% Rxd sig
    end % Generate two disn
    R = [z,y]; % Rxd sig and noise over complete duration
    D = zeros(length(SNR),nc);
    for k =1:nc
%        D(:,k) = sum((abs(R(:,((k-1)*nb)+1:k*nb))).^2,2);
%        D(:,k) = max((abs(R(:,((k-1)*nb)+1:k*nb))).^2,[],2);
        D(:,k) = mean((abs(R(:,((k-1)*nb)+1:k*nb))).^2,2);
    end
    TS = (D .* repmat(ts1,1,nc)) + repmat(ts2,1,nc) + log(1/beta);
    STs = zeros(size(SNR));

    %% PFD and PLD
    
    PFD = sum(double(TS(:,1:round(ph0*nc)) >= repmat(STs,1,round(ph0*nc))),2)/(ph0*nc);
    PLD = sum(double(TS(:,(round(ph0*nc))+1:end) <= repmat(STs,1,nc-(round(ph0*nc)))),2)/(nc-(round(ph0*nc)));
    APeso(:,l+1) = smooth(smooth(smooth((ph0*PFD)+(ph1*PLD))));
    if ph0 == 0.5
        PF = smooth(smooth(smooth(PFD))); PL = smooth(smooth(smooth(PLD)));
    end
end
%subplot(2,1,1)
plot(SNR_db,APeso(:,1),'g',SNR_db,APeso(:,2),'--r',SNR_db,APeso(:,3),':k','LineWidth',5);
title('P_e(\lambda_{Soft}) for Changing Priors','fontsize',16,'fontweight','b');
legend( 'P(h_0)=0.7', 'P(h_0)=0.5', 'P(h_0)=0.3');ylim([0 1]);
xlabel('SNR','fontsize',16,'fontweight','b');
ylabel('P_e','fontsize',16,'fontweight','b');

figure(2);
semilogy(SNR_db,PF,SNR_db,PL,'g','LineWidth',5);
title('ROC','fontsize',16,'fontweight','b');
xlabel('PLD','fontsize',16,'fontweight','b');
ylabel('PFD','fontsize',16,'fontweight','b');