clc; clear all; close all;

%% Assumptions
mu1 = 50; % signal mean original 300, 2 for anp, 0 for soft
nb = 60; % Buffer size. For Soft threshold nb = 1 beause each sample is tested
nc = 1*100; nos = nc*nb; % no of channels (multiple of 100) and samples
s1 =0; mu0=0; % Noise is zero mean and sig is zero std dev
SNR_db = -15:1:15; SNR = (db2pow(SNR_db))'; % SNR is mu1^2 / s0^2
s0 = sqrt((mu1^2)./SNR); % sigma for noise according to SNR
sD0 = sqrt(2*nb * (s0.^4));% Var for squared noise samples is 4(mu0^2)(s0.^2)+2(s0.^4)

%% Analysis with changing Beta
PA = 0.1; % Acceptable Pfa
Bvar = 2; % Beta variations implying changing priors
APeAN = zeros(length(SNR),Bvar+1);
for l = 0:Bvar
    %% Generating Test Data R and decision Metric D
    ph1 = 0.3+l*(0.2);ph0 = 1-ph1;% Priors. ph0 should be greater than ph1
    beta = ph0/ph1; % Current Beta
    
    STanp = s0.^2 * chi2inv(PA/ph0,nc);%For mu1 = 20, nb = 10
%    STanp = sD0.^2 * qfuncinv(PA/ph0);%For mu1 = 20, nb = 10
%    STanp = s0.^2 * qfuncinv(PA/ph0);
    z = zeros(length(SNR),round(ph0*nos));y = zeros(length(SNR),round(ph1*nos));s=y;
    for i = 1:length(SNR)
        z(i,:) = normrnd(mu0,s0(i),1,round(ph0*nos));% Generate noise. 
        s(i,:) = normrnd(mu1,s1,1,round(ph1*nos));% Rxd sig without noise
        y(i,:) = s(i,:)+ (normrnd(mu0,s0(i),1,round(ph1*nos)));% Rxd sig
    end % Generate two disn
    R = [z,y]; % Rxd sig and noise over complete duration
    D = zeros(length(SNR),nc);D2=D;
    for k =1:nc
        D(:,k) = sum((abs(R(:,((k-1)*nb)+1:k*nb))).^2,2);
    end

    %% PFD and PLD
    PFD5 = sum(double(D(:,1:round(ph0*nc)) >= repmat(STanp,1,round(ph0*nc))),2)/(ph0*nc);
    PLD5 = sum(double(D(:,(round(ph0*nc))+1:end) <= repmat(STanp,1,nc-(round(ph0*nc)))),2)/(nc-(round(ph0*nc)));
    APeAN(:,l+1) = smooth(smooth(smooth((ph0*PFD5)+(ph1*PLD5))));
%    plot(STanp); hold on;
end

plot(SNR_db,APeAN(:,1),'g',SNR_db,APeAN(:,2),'--r',SNR_db,APeAN(:,3),':k','LineWidth',5);
title('P_e(\lambda_{ANP}) for Changing Priors','fontsize',16,'fontweight','b');
legend( 'P(h_0)=0.7', 'P(h_0)=0.5', 'P(h_0)=0.3');ylim([0 1]);
xlabel('SNR','fontsize',16,'fontweight','b');
ylabel('P_e','fontsize',16,'fontweight','b');
