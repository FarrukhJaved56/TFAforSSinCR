clc; clear all; close all

%% Assumptions
mu = 9; % mu1 = 3 and nb = 1
nb = 1; % Buffer size. For Soft threshold nb = 1 beause each sample is tested
nc = 20*100; nos = nc*nb; % no of channels (multiple of 100) and samples
s1 =0; mu0=0; % Noise is zero mean and sig is zero std dev
SNR_db = -15:1:15; SNR = (db2pow(SNR_db))'; % SNR is mu1^2 / s0^2
s0 = sqrt((mu^2)./SNR); % sigma for noise according to SNR
sD1 = sqrt(2*nb*((2*mu^2 * (s0.^2)) + (s0.^4)));% Var for squared sig samples is 
sD0 = sqrt(2*nb * (s0.^4));% Var for squared noise samples is 4(mu0^2)(s0.^2)+2(s0.^4)

%% STs
pe5 = ((sD1.^2).*(sD1.^2./((2*sD0.^2).*(sD0.^2+sD1.^2))))+(0.5*log(sD0.^2./(sD0.^2-sD1.^2)));
% mu1 = 3 and nb = 1
% mu = 9, mu1 = 3 and nb = 1

%pe52 = (s1.^2./((2*s0.^2).*(s0.^2+s1.^2)));
%pe53 = (0.5*log(s0.^2./(s0.^2-s1.^2)));
%pe54 = (sD1.^2./((2*sD0.^2).*(sD0.^2+sD1.^2)));
%pe55 = (0.5*log(sD0.^2./(sD0.^2-sD1.^2)));


%% Analysis with changing Beta
Bvar = 2; % Beta variations implying changing priors
APeAN = zeros(length(SNR),Bvar+1);APeso = APeAN;APeso2 = APeso;
for l = 0:Bvar
    %% Generating Test Data R and decision Metric D
    ph1 = 0.3+l*(0.2);ph0 = 1-ph1;% Priors. ph0 should be greater than ph1
    beta = ph0/ph1; % Current Beta
   
    % similar to 100 times inc in nb but not possible there because out of memory
    STs = pe5 + log(1/beta);

    z = zeros(length(SNR),round(ph0*nos));y = zeros(length(SNR),round(ph1*nos));s=y;
    for i = 1:length(SNR)
        z(i,:) = normrnd(mu0,s0(i),1,round(ph0*nos));% Generate noise. 
        s(i,:) = normrnd(mu,s1,1,round(ph1*nos));% Rxd sig without noise
        y(i,:) = s(i,:)+ (normrnd(mu0,s0(i),1,round(ph1*nos)));% Rxd sig
    end % Generate two disn
    R = [z,y]; % Rxd sig and noise over complete duration

    D = zeros(length(SNR),nc);
    for k =1:nc
        D(:,k) = sum((abs(R(:,((k-1)*nb)+1:k*nb))).^2,2);
    end
%    STs2 = (D .* repmat(pe52,1,nc)) + repmat(pe53,1,nc) + log(1/beta);

    %% PFD and PLD
    
    PFD6 = sum(double(D(:,1:round(ph0*nc)) >= repmat(STs,1,round(ph0*nc))),2)/(ph0*nc);
    PLD6 = sum(double(D(:,(round(ph0*nc))+1:end) <= repmat(STs,1,nc-(round(ph0*nc)))),2)/(nc-(round(ph0*nc)));
    APeso(:,l+1) = smooth(smooth(smooth((ph0*PFD6)+(ph1*PLD6))));
    
%    PFD7 = sum(double(D(:,1:round(ph0*nc)) >= STs2(:,1:round(ph0*nc))),2)/(ph0*nc);
%    PLD7 = sum(double(D(:,(round(ph0*nc))+1:end) <= STs2(:,(round(ph0*nc))+1:end)),2)/(nc-(round(ph0*nc)));
%    APeso2(:,l+1) = smooth(smooth(smooth((ph0*PFD7)+(ph1*PLD7))));

end
%subplot(2,1,1)
plot(SNR_db,APeso(:,1),'g',SNR_db,APeso(:,2),'--r',SNR_db,APeso(:,3),':k','LineWidth',5);
title('P_e(\lambda_{Soft}) for Changing Priors','fontsize',16,'fontweight','b');
legend( 'P(h_0)=0.7', 'P(h_0)=0.5', 'P(h_0)=0.3');ylim([0 1]);
xlabel('SNR','fontsize',16,'fontweight','b');
ylabel('P_e','fontsize',16,'fontweight','b');

%subplot(2,1,2)
%plot(SNR_db,APeso2(:,1),'g',SNR_db,APeso2(:,2),'--r',SNR_db,APeso2(:,3),':k','LineWidth',3);
%title('P_e(\lambda_{Soft}) for Changing Priors','fontsize',16,'fontweight','b');
%legend( 'P(h_0)=0.7', 'P(h_0)=0.5', 'P(h_0)=0.3');ylim([0 1]);
%xlabel('SNR','fontsize',16,'fontweight','b');
%ylabel('P_e','fontsize',16,'fontweight','b');
