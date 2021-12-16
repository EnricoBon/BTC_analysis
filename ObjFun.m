function [RMSE,r2,nRMSE,logRMSE,logr2,Pearson_r2,Pearson_logr2,KGE] = ObjFun(BTC_input,CC1)

N=length(BTC_input(:,1));       % number of observations
C_peak=max(BTC_input(:,2));     % Concentration peak

clear RMSE r2 nRMSE logRMSE logr2 Pearson_r2 Pearson_logr2 KGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get rid of negative values
for lll=1:1:length(CC1(:,1))
    if CC1(lll,1)<0
        CC1(lll,1)=0;
    end
end

clear lll

% Objective function -> RMSE and nRMSE
for k=1:1:length(CC1(1,:))
    for i=1:1:length(CC1(:,1))
        ssd(i,k)=(CC1(i,k)-BTC_input(i,2))^2;   % squared differences (simulated vs observed)
    end
end

for i=1:1:length(CC1(1,:))
    RMSE(1,i)=(sum(ssd(:,i))/N)^(1/2);    % Root mean squared error for every simulated BTC
    nRMSE(1,i)=RMSE(1,i)/C_peak;          % normalized RMSE
end

% Objective function -> logRMSE
% Note that some values of the BTC can be = 0. In this case the log(0)=-inf
% which does not allow us to obtain a finite value for logRMSE and log r^2
% to avoid the problem we can do two things
% 1) following the code from Ward et al (2017) we get rid of conc values=0.
% 2) We set 0values equal to 0.0001.

clear ssd
BTC_input_log(:,1)=BTC_input(:,1);  % keep the time
for i=1:1:length(CC1(:,1))
        if BTC_input(i,2)==0
           BTC_input_log(i,2)=0.0001; 
        else
           BTC_input_log(i,2)=BTC_input(i,2); 
        end
end
for k=1:1:length(CC1(1,:))
    for i=1:1:length(CC1(:,1))
        if CC1(i,k)==0
           CC1(i,k)=0.0001; 
        end
    end
end

% Evaluate logRMSE
for k=1:1:length(CC1(1,:))
    for i=1:1:length(CC1(:,1))
        ssd(i,k)=(log(CC1(i,k))-log(BTC_input_log(i,2)))^2;   % squared differences (log simulated vs log observed)
    end
end
for i=1:1:length(CC1(1,:))
    logRMSE(1,i)=(sum(ssd(:,i))/N)^(1/2);    % log RMSE
end

% Objective function -> r^2 -> NOTE THAT THIS r^2 is == to Nash Sutcliffe
% efficiency

MeanObs=mean(BTC_input(:,2));
for i=1:1:length(BTC_input(:,1))
    Diff2(i,1)=(BTC_input(i,2)-MeanObs(1,1))^2;
end
Denominator=sum(Diff2(:,1));

for k=1:1:length(CC1(1,:))
    for i=1:1:length(CC1(:,1))
        ssd(i,k)=(CC1(i,k)-BTC_input(i,2))^2;   % squared differences (simulated vs observed)
    end
    Numerator(1,k)=sum(ssd(:,k));
end

for i=1:1:length(CC1(1,:))
    r2(1,i)=1-Numerator(1,i)/Denominator;       % r^2 
end

% Objective function -> log r^2
% use again the "corrected" values (0.0001 instead of 0) to avoid "-inf"
% problems
clear Numerator Denominator Diff2 ssd

MeanObs=log(mean(BTC_input(:,2)));
for i=1:1:length(BTC_input(:,1))
    Diff2(i,1)=(log(BTC_input_log(i,2))-log(MeanObs(1,1)))^2;
end
Denominator=sum(Diff2(:,1));

for k=1:1:length(CC1(1,:))
    for i=1:1:length(CC1(:,1))
        ssd(i,k)=(log(CC1(i,k))-log(BTC_input_log(i,2)))^2;   % squared differences (simulated vs observed)
    end
    Numerator(1,k)=sum(ssd(:,k));
end

for i=1:1:length(CC1(1,:))
    logr2(1,i)=1-Numerator(1,i)/Denominator;       % log r^2 
end

% Objective function -> Pearson r^2 and log Pearson r^2
for k=1:1:length(CC1(1,:))
    R=corrcoef(CC1(:,k),BTC_input(:,2));
    Pearson_r2(1,k)=R(1,2);
    RR=corrcoef(log(CC1(:,k)),log(BTC_input_log(:,2)));
    Pearson_logr2(1,k)=RR(1,2);
end
clear R RR ssd 

% Objective function -> Kling Gupta
for k=1:1:length(Pearson_r2(1,:))
KGE(1,k)=1-((Pearson_r2(1,k)-1).^2 + (std(CC1(:,k))/std(BTC_input(:,2)) -1).^2 + (mean(CC1(:,k))/mean(BTC_input(:,2)) -1).^2);
end

end
