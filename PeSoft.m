clc; clear all; close all

%% Assumptions
nc = 100; % no of channels (multiple of 100) and samples
s1 =6; % Noise is zero mean and sig is zero std dev
SNR_db = -15:1:15; SNR = (db2pow(SNR_db))'; % SNR is mu1^2 / s0^2
s0 = sqrt((s1^2)./SNR); % sigma for noise according to SNR
%s0=1;
%% Distributions
pe52 = s1.^2 ./((2*s0.^2).*(s0.^2+s1.^2));
pe53 = 0.5*log(s0.^2./(s1.^2-s0.^2));

%% Analysis with changing Beta
Bvar = 2; % Beta variations implying changing priors
APeso = zeros(length(SNR),Bvar+1);

for l = 0:Bvar
    %% Generating Test Data R and decision Metric D
    ph1 = 0.3+l*(0.2);ph0 = 1-ph1;% Priors. ph0 should be greater than ph1
    beta = ph0/ph1; % Current Beta

    z = zeros(length(SNR),nc);y = z;s=y;
    for i = 1:length(SNR)
        z(i,:) = normrnd(0,s0(i),1,nc);% Generate noise. 
        s(i,:) = normrnd(0,s1,1,nc);% Rxd sig without noise
        y(i,:) = s(i) + normrnd(0,s0(i),1,nc);% Rxd sig
    end % Generate two disn
%    R = [z,y].^2; % Rxd sig and noise over complete duration
%    STs = (R .* repmat(pe52,1,nc)) + repmat(pe53,1,nc) - log(beta);
    STsz = (abs(z).^2 .* repmat(pe52,1,nc)) + repmat(pe53,1,nc) - log(beta);
    STsy = (abs(y).^2 .* repmat(pe52,1,nc)) + repmat(pe53,1,nc) - log(beta);

    %% PFD and PLD
%    PFD6 = sum(double(R(:,1:round(ph0*nc)) >= STs(:,1:round(ph0*nc))),2)/(ph0*nc);
%    PLD6 = sum(double(R(:,(round(ph0*nc))+1:end) <= STs(:,(round(ph0*nc))+1:end)),2)/(nc-(round(ph0*nc)));
    PFD6 = sum(double(z.^2 >= log(STsz)),2)/nc;
    PLD6 = sum(double(y.^2 <= log(STsy)),2)/nc;
    APeso(:,l+1) = smooth(smooth(smooth((ph0*PFD6)+((1-ph0)*PLD6))));
    
end
plot(SNR_db,APeso(:,1),'g',SNR_db,APeso(:,2),'--r',SNR_db,APeso(:,3),':k','LineWidth',3);
title('P_e(\lambda_{Soft}) for Changing Priors','fontsize',16,'fontweight','b');
legend( 'P(h_0)=0.7', 'P(h_0)=0.5', 'P(h_0)=0.3');ylim([0 1]);
xlabel('SNR','fontsize',16,'fontweight','b');
ylabel('P_e','fontsize',16,'fontweight','b');
