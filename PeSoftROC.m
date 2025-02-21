clc; clear all; close all

%% Assumptions
nb = 1; % Buffer size. For Soft threshold nb = 1 beause each sample is tested
nc = 20*1000; nos = nc*nb; % no of channels (multiple of 100) and samples
s1 = 36;s0 = 1; 
SNR = (db2pow(s1^2/ s0^2))'; % SNR is mu1^2 / s0^2
ph1 = 0.5;ph0 = 1-ph1;% Priors
beta = ph0/ph1; % Current Beta

%% Test statistic
ts1 = (s1^2./((2*s0.^2).*(s0.^2+s1^2)));
ts2 = (0.5*log(s0.^2./(s0.^2+s1^2)));
%% Sig and noise
z = zeros(length(SNR),round(ph0*nos));y = zeros(length(SNR),round(ph1*nos));s=y;
for i = 1:length(SNR)
    z(i,:) = normrnd(0,s0(i),1,round(ph0*nos));% Generate noise. 
    s(i,:) = normrnd(0,s1,1,round(ph1*nos));% Rxd sig without noise
    y(i,:) = s(i,:)+ (normrnd(0,s0(i),1,round(ph1*nos)));% Rxd sig
end % Generate two disn
R = [z,y]; % Rxd sig and noise over complete duration
D = zeros(length(SNR),nc);
for k =1:nc
%    D(:,k) = sum((abs(R(:,((k-1)*nb)+1:k*nb))).^2,2);
%    D(:,k) = max((abs(R(:,((k-1)*nb)+1:k*nb))).^2,[],2);
    D(:,k) = mean((abs(R(:,((k-1)*nb)+1:k*nb))).^2,2);
end
%TS = (D .* repmat(ts1,1,nc)) + repmat(ts2,1,nc) + log(1/beta);
TS = (D * ts1) + ts2 + log(1/beta);
%% Analysis with changing Beta
res = 100000;PFD=zeros(1,length(res)-1); PLD= PFD;
for ctr = 1:res-1
    %% Generating Test Data R and decision Metric D
    pfa = ctr/res;
    ST = s0 * qfuncinv(pfa);
    %% PFD and PLD
    PFD(ctr) = sum(double(TS(:,1:round(ph0*nc)) >= ST),2)/(ph0*nc);
    PLD(ctr) = sum(double(TS(:,(round(ph0*nc))+1:end) <= ST),2)/(nc-(round(ph0*nc)));

end
loglog(PFD,PLD,'k','LineWidth',5);
title('ROC','fontsize',16,'fontweight','b');
xlabel('P_{FA}','fontsize',16,'fontweight','b');
ylabel('P_{MD}','fontsize',16,'fontweight','b');