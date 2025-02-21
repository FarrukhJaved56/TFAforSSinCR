clc; clear all; close all

%% Assumptions
mu1 = 3; % signal mean original 300, 2 for anp, 0 for soft
nb = 10; % Buffer size. For Soft threshold nb = 1 beause each sample is tested
nc = 20*100; nos = nc*nb; % no of channels (multiple of 100) and samples
s1 =0; mu0=0; % Noise is zero mean and sig is zero std dev
SNR_db = -15:1:15; SNR = (db2pow(SNR_db))'; % SNR is mu1^2 / s0^2
s0 = sqrt((mu1^2)./SNR); % sigma for noise according to SNR
%s0 = ones(size(SNR));
%% Distributions
muD0 = nb * ((s0.^2));% mean for squared noise samples is (s0.^2)+(mu0.^2)
muD1 = nb* (mu1^2 + s0.^2);% Mean of nb squared nb samples (as above)
sD0 = sqrt(2*nb * (s0.^4));% Var for squared noise samples is 4(mu0^2)(s0.^2)+2(s0.^4)
sD1 = sqrt(2*nb*((2*mu1^2 * (s0.^2)) + (s0.^4)));% Var for squared sig samples is 

%% STs
PA = 0.1; % Acceptable Pfa

%pe5 = ((muD1.^2).*(s1^2./((2*s0.^2).*(s0.^2+s1^2))))+(0.5*log(s0.^2./(s0.^2-s1^2)));
pe5 = ((sD1.^2).*(sD1.^2./((2*sD0.^2).*(sD0.^2+sD1.^2))))+(0.5*log(sD0.^2./(sD0.^2-sD1.^2)));

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
   
%    STanp = s0.^2 * qfuncinv(PA/ph0);
%    STanp = sD0.^2 * qfuncinv(PA/ph0);% for all nb!=1
    STanp = 10 * sD0.^2 * qfuncinv(PA/ph0);% for all nb!=1
    % similar to 100 times inc in nb but not possible there because out of memory
    STs = pe5 + log(1/beta);

    z = zeros(length(SNR),round(ph0*nos));y = zeros(length(SNR),round(ph1*nos));s=y;
    for i = 1:length(SNR)
        z(i,:) = normrnd(mu0,s0(i),1,round(ph0*nos));% Generate noise. 
        s(i,:) = normrnd(mu1,s1,1,round(ph1*nos));% Rxd sig without noise
        y(i,:) = s(i,:)+ (normrnd(mu0,s0(i),1,round(ph1*nos)));% Rxd sig
    end % Generate two disn
    R = [z,y]; % Rxd sig and noise over complete duration
    R2 = [z,(1/mu0)*y]; % Normalized for single sample case

    D = zeros(length(SNR),nc);D2=D;
    for k =1:nc
        D(:,k) = sum((abs(R(:,((k-1)*nb)+1:k*nb))).^2,2);
        D2(:,k) = sum((abs(R2(:,((k-1)*nb)+1:k*nb))).^2,2);
    end
%    STs2 = (D .* repmat(pe52,1,nc)) + repmat(pe53,1,nc) + log(1/beta);

    %% PFD and PLD
    PFD5 = sum(double(D(:,1:round(ph0*nc)) >= repmat(STanp,1,round(ph0*nc))),2)/(ph0*nc);
    PLD5 = sum(double(D(:,(round(ph0*nc))+1:end) <= repmat(STanp,1,nc-(round(ph0*nc)))),2)/(nc-(round(ph0*nc)));
    APeAN(:,l+1) = smooth(smooth(smooth((ph0*PFD5)+(ph1*PLD5))));

    PFD6 = sum(double(D2(:,1:round(ph0*nc)) >= repmat(STs,1,round(ph0*nc))),2)/(ph0*nc);
    PLD6 = sum(double(D2(:,(round(ph0*nc))+1:end) <= repmat(STs,1,nc-(round(ph0*nc)))),2)/(nc-(round(ph0*nc)));
    APeso(:,l+1) = smooth(smooth(smooth((ph0*PFD6)+(ph1*PLD6))));
    
%    PFD7 = sum(double(D(:,1:round(ph0*nc)) >= STs2(:,1:round(ph0*nc))),2)/(ph0*nc);
%    PLD7 = sum(double(D(:,(round(ph0*nc))+1:end) <= STs2(:,(round(ph0*nc))+1:end)),2)/(nc-(round(ph0*nc)));
%    APeso(:,l+1) = smooth(smooth(smooth((ph0*PFD7)+(ph1*PLD7))));

end
figure(5)
subplot(2,1,1);
plot(SNR_db,APeAN(:,1),'g',SNR_db,APeAN(:,2),'--r',SNR_db,APeAN(:,3),':k','LineWidth',3);
title('P_e(\lambda_{ANP}) for Changing Priors','fontsize',16,'fontweight','b');
legend( 'P(h_0)=0.7', 'P(h_0)=0.5', 'P(h_0)=0.3');ylim([0 1]);
xlabel('SNR','fontsize',16,'fontweight','b');
ylabel('P_e','fontsize',16,'fontweight','b');

%figure(6)
subplot(2,1,2);
plot(SNR_db,APeso(:,1),'g',SNR_db,APeso(:,2),'--r',SNR_db,APeso(:,3),':k','LineWidth',3);
title('P_e(\lambda_{Soft}) for Changing Priors','fontsize',16,'fontweight','b');
legend( 'P(h_0)=0.7', 'P(h_0)=0.5', 'P(h_0)=0.3');ylim([0 1]);
xlabel('SNR','fontsize',16,'fontweight','b');
ylabel('P_e','fontsize',16,'fontweight','b');
