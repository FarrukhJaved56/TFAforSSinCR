%% Confirm the correctness of prior estimation
clc
clear all
close all
%% Simulate a frame
f = 100;% No of frames
c = 10;% No of channels
toc = 20; % Time of occupation
Ph1_actual = zeros(11,toc);Ph0_actual= zeros(11,toc);
Ph1_calc = zeros(11,toc);Ph0_calc= zeros(11,toc);
%% Inc time of occupation
for l = 1:160
    ph1 = 0.005*(l-1); % Simulated prior probability of presence
    ph0 = 1 - ph1; % Simulated prior probability of presence
    sim = randsrc(c,f,[1,0; ph1,ph0]);
    simi=sim;simij=sim;
    % Generate random matrix with each row is a channel and colms are frames
    for k = 1:toc % Min time of occupation
        for i=1:c
            for j=1:f
                if sim(i,j) == 1 && j < f-(k-1)
                    simi(i,j:j+k) = 1;
                elseif sim(i,j) == 1 && j >= f-(k-1)
                    simi(i,j:end) = 1;
                end
                simij(i,j) = j*simi(i,j);
            end
        end
        Ph1_actual(l,k) = mean(sum(simi,2)/f);
        Ph0_actual(l,k) = 1-Ph1_actual(l,k);
        simi=zeros(size(sim));
        Ph1_calc(l,k) = mean(sum(simij,2)*(2/(f*(1+f))));
        Ph0_calc(l,k) = 1-Ph1_calc(l,k);
    end
end
%Error = abs(Ph1_calc-Ph1_actual);
Error = abs(Ph0_calc-Ph0_actual);
for i = 1:size(Error,2)
    Error(:,i) = smooth(smooth(smooth(Error(:,i))));
end
X = 0:1/(l-14):1-(1/l);Y=1:20; 
surf(Y,X,Error(15:end,:))
figure(2)
plot(0:1/160:1-1/160,Error)