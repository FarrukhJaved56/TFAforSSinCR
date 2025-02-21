clc; clear all; close all

%% Variables
a = 0.1;g = 1.5;sz2 = 1;mus=10;
ph0 = 0.5; ph1 = 1 - ph0;b= ph0/ph1;
nb = 1:.1:300;nb2 = sqrt(nb/2); % Prev Figs in thesis are on 300 nb
SNR_db = -24:1:25; SNR = (db2pow(SNR_db))'; K = 1./sqrt((2*SNR)+1);

%% NP and ANP
penp = zeros(length(SNR),length(nb));
peanp = zeros(length(SNR),length(nb));
for ctr= 1:length(SNR)
    AQnp = -K(ctr)*(qfuncinv(a)-(nb2*SNR(ctr)));% Arg Q Pe
    AQanp = -K(ctr)*((qfuncinv(a/ph0)./(2*nb2))-(nb2*(1+SNR(ctr))));% Arg Q Pe

    Qnp = 1-qfunc(-AQnp);
    Qanp = 1-qfunc(-AQanp);
    
    penp(ctr,:) = (a*ph0)+(ph1*Qnp);
    peanp(ctr,:) = (a*ph0)+(ph1*Qanp);
end

%% SEF
pesef = zeros(length(SNR),length(nb));
for ctr= 1:length(SNR)
    SQ = (((g-1)/2)*(SNR(ctr))^2) +(1/(K(ctr))^2);
    AQsef = nb2*SNR(ctr)/(1+sqrt(SQ));% Arg Q Pe
    pesef(ctr,:) = qfunc(AQsef);
end

%% PE
pepe = zeros(length(SNR),length(nb));
for ctr= 1:length(SNR)
    SNR1 = (1/sz2)*(1+(2/SNR(ctr)));
    SNR2 = SNR1^2;
    SNR3 = (nb.^2/2)*((sz2^2)*SNR(ctr));
    SNR4 = (nb*sz2^2)*(log((b/K(ctr))^2))*(2+(1/SNR(ctr)));
    STpe = SNR1+sqrt(SNR2+SNR3+SNR4);
    
    QA1 = STpe./(sz2*sqrt(2*nb));
    Q1 = qfunc(QA1-nb2);
    Q2 = qfunc(-K(ctr).*(QA1-(sqrt(nb/2)*(SNR(ctr)+1))));
    pepe(ctr,:) = (ph0*Q1)+(ph1*Q2);
end

%% Opt
peopt = zeros(length(SNR),length(nb));
for ctr= 1:length(SNR)
    SNR5 = 1./ ((2*SNR(ctr)) +(((g - 1)/2)*SNR(ctr).^2) );
    SNR6 = nb *(1-(SNR(ctr)*SNR5));
    SNR7 = sqrt(1+SNR5);
    SNR8 = sqrt((2*nb*log(b.^2))+(((nb*SNR(ctr)).^2) * SNR5));
    STopt = sz2*(SNR6+(SNR7*SNR8));
    
    QA1 = STopt./(sz2*sqrt(2*nb));
    Q1 = qfunc(QA1-nb2);
    Q2 = qfunc(-K(ctr).*(QA1-(sqrt(nb/2)*(SNR(ctr)+1))));
    peopt(ctr,:) = (ph0*Q1)+(ph1*Q2);
end

%% Soft
pesoft = zeros(length(SNR),length(nb));
for ctr= 1:length(SNR)
%    STsoft = mus^2*((nb*SNR(ctr))/(2*(1+SNR(ctr))))+(0.5/(1-SNR(ctr)))-log(b);
    STsoft = ((nb*SNR(ctr).^2)/(2*(1+SNR(ctr))))+(0.5/(1-SNR(ctr)))-log(b);
    
    QA1 = STsoft./(sz2*sqrt(2*nb));
    Q1 = qfunc(QA1-nb2);
    Q2 = qfunc(-K(ctr).*(QA1-(sqrt(nb/2)*(SNR(ctr)+1))));
    pesoft(ctr,:) = (ph0*Q1)+(ph1*Q2);
end
%% Figs
%figure(1)
subplot(3,1,1);mesh(penp);colorbar;xlabel('N_B');ylabel('SNR');title('P_e(\lambda_{NP})');
subplot(3,1,2);mesh(pesef);colorbar;xlabel('N_B');ylabel('SNR');title('P_e(\lambda_{SEF})');
subplot(3,1,3);mesh(pepe);colorbar;xlabel('N_B');ylabel('SNR');title('P_e(\lambda_{P_e})');
%figure(2);mesh(peopt);colorbar;xlabel('N_B');ylabel('SNR');title('P_e(\lambda_{Opt})');
figure(3)
subplot(2,1,1);mesh(pesoft);colorbar('fontsize',14,'fontweight','b');xlabel('N_B','fontsize',16,'fontweight','b');
ylabel('SNR','fontsize',16,'fontweight','b');title('P_e(\lambda_{soft})','fontsize',16,'fontweight','b');
subplot(2,1,2);mesh(peanp);colorbar('fontsize',14,'fontweight','b');xlabel('N_B','fontsize',16,'fontweight','b');
ylabel('SNR','fontsize',16,'fontweight','b');title('P_e(\lambda_{ANP})','fontsize',16,'fontweight','b');