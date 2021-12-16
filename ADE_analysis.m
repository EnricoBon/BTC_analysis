function [ADE1, ADE2, ADE3, ADE3_limits,Hyperspace,M_g,Q1,time]=ADE_analysis(BTC_input,L, M, Qmeasured,Description1,Description2, n,BTC_result)

% In this function we calibrate 3 different Advection-dispersion equation
% and evaluate their performances with the observed BTC 

clear ADE1 ADE2 val idx C_fin1 C_fin2 D_fin1 D_fin2 CC1 logr2 logRMSE nRMSE 
clear Numerator j i k ssd RMSE r2 r2_1 A CC1 Diff2 RMSE_1 S ts v v1 M_lostPerc
clear M_recovered MeanObs R M_g Q1 Q B B_length Pearson_logr2 Pearson_r2

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                            _____  ______   __                       % %
% %                      /\   |  __ \|  ____| /_ |                      % %
% %                     /  \  | |  | | |__     | |                      % %
% %                    / /\ \ | |  | |  __|    | |                      % %
% %                   / ____ \| |__| | |____   | |                      % %
% %                  /_/    \_\_____/|______|  |_|                      % %
% %                                                                     % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% v is fixed and equal to L/t_peak - 
% Q is calculated from the BTC (Injected mass/sum of Concentration) ->
% dilution gauging method (Anderson and Burt, 1978 - The role of topography
% in controlling throughflow generation)
% A is calculated = Q/v 
% D is calibrated (step calibration).
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

Length=length(BTC_input(:,1));
ts=BTC_input(3,1)-BTC_input(2,1);       % Time step in hours

for i=1:1:Length
    CL_t(i,1)=BTC_input(i,2)*ts;   % (g*hours)/m^3
end

S=sum(CL_t);    % sum of concentration in (g*hours)/m^3

clear CL_t

M_g=M/1000;     % mass in gram (M is in mg)
Q=M_g/S;        % flow [m^3/hrs]
Q1=Q/3.6;       % flow [l/s]

v=L/BTC_result.tpeak;   %% streamflow velocity [m/h]
v1=v/3600;              %% streamflow velocity [m/s]

A=Q/v;                % m^2

% Advection dispersion equation
% we need time in seconds for the D
time(:,1)=3600.*BTC_input(:,1);          % from hours to seconds

% pre-allocating D
D=0.0001; % [m^2/s] advection-dispersion coefficient

% Deciding an interval of calibration fod D --> 0.0001 accuracy
% This interval (especially the upper limit) has to be chosen from the 
% modeller, depending on the investigated reach. 
B=[0.0001:0.0001:2];
B=B.';
B_length=length(B);
k=1;
% Solving ADE for each D
for i=1:1:B_length
    D=B(i,1);
    for j=1:1:Length
    CC1(k,i)=(M_g/((A*(4*3.14*D*time(j,1))^(1/2))))*...
        exp(-((L-v1*time(j,1))^2/(4*D*time(j,1))));     % g - m - s -> Conc in g/m3
    k=k+1;
    end
    clear k
    k=1;
end

clear k i

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% 
%%%% %%%% %%%% %%%% %%%% %%%% OBJECTIVE FUNCTION %%%% %%%% %%%% %%%% %%%% 
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%  

% Evaluate, for every combination, the performances using different
% objective functions. Note that every information of the initial BTC could
% potentially be used as objective function (skewness, holdback, apparent
% dispersion, c_peak and so on...) despite they could also be meaningless (Ward et al., 2017) 
% Here we will evaluate the performances via:
%
% > RMSE -> The most used objective function to assess BTC errors in stream 
% corridor works since it is equivalent form of RSS (used in OTIS).
%%% Runkel, R. L. (1998). One‐dimensional transport with inflow and storage (OTIS): A solute transport model for streams and Rivers 
%%% Marion, A., M. Zaramella, and A. Bottacin-Busolin (2008), Solute transport in rivers with multiple storage zones: The STIR model
%%% Kelleher, C., et al. (2013). Identifiability of transient storage model parameters along a mountain stream
%%% González-Pinzón, R., et al. (2015). A field comparison of multiple techniques to quantify groundwater–surface-water interactions
%%% Ward, A. S. et al. (2017). A software tool to assess uncertainty in transient‐storage model parameters using Monte Carlo simulations
%%% Kelleher, C. et al. (2019). Exploring tracer information and model framework trade‐offs to improve estimation of stream transient storage processes
%
% > r^2 -> RMSE is not a normalized error metric, so the coefficient of 
% determination r^2 assesses the fraction of variability in the observed 
% BTC that is described by variability in the simulated BTC -> Since it is
% normalized it can be used to compare different reaches
%%% Kelleher, C., et al. (2013). Identifiability of transient storage model parameters along a mountain stream parameters along a mountain stream  
%%% Definition of r^2 from --> Devore, J. L. (2000), Probability and Statistics for Engineering and the Sci- ences, Duxbury, Pacific Grove, Calif
% Note that this r^2 is different from Pearson's correlation (often called
% r^2 as well). r^2 correlation coefficient CAN BE HIGHLY NEGATIVE, on the 
% contrary than Pearson's correlation which is between -1 and +1
% r^2 <0 simply means the model fits worse than a linear model equal to the
% mean value C_mean of the BTC -> the simulated BTC is complete garbage.
% Consider that the mathematical formulation of r^2 is identical to the
% mathematical of Nash-Sutcliffe efficiency.
% To avoid misunderstandings I'll also evaluate Pearson's correlation
% coeffieicnt and call it PearsonR2

% For every BTC I'll also calculate: 
% > Kling–Gupta efficiency (KGE) -> Increasingly more used than r^2 (NSE)
% since it can be more informative than NSE (Knoben et al., 2019);
% > nRMSE -> RMSE normalized with the respect of peak concentration. It's
% useful when dealing with different peak concentrations -> different
% tracers -> not used in my application but could be useful for other
% studies
% > log-transformed RMSE and log-transformed r^2
% RMSE and r^2 evaluated after the observation and the simulation values of
% the contentration have been converted in logaritmic values
% log-transformed error metrics can be useful for estimating TSM parameters:
%%% Ward, A. S. et al. (2017). A software tool to assess uncertainty in transient‐storage model parameters using Monte Carlo simulations
%%% Wörman, A., & Wachniew, P. (2007). Reach scale and evaluation methods as limitations for transient storage properties in streams and rivers
%%% Kelleher, C. et al. (2019). Exploring tracer information and model framework trade‐offs to improve estimation of stream transient storage processes

%%%%% Okay, our theory is done, let's keep coding:
 
[RMSE,r2,nRMSE,logRMSE,logr2,Pearson_r2,Pearson_logr2,KGE] = ObjFun(BTC_input,CC1);

% Organize the results

ADE1.InjectedMass=M_g;      % injected Cl mass;
ADE1.L=L;                   % Length of the stream reach [m]
ADE1.Q=Q1;                  % Discharge l/s;
ADE1.A=A;                   % Wetted area A [m^2]
ADE1.v=v1;                  % velocity [m/s]

% Preallocating the ADE1.D
ADE1.D(1:(length(B(:,1))+1),1:8)={'a'};     % all the obj. fun.


% list of all the D-values used with the relative obj functions
ADE1.D(1,1)={'D'};          
ADE1.D(2:end,1)=num2cell(B);
ADE1.D(1,2)={'RMSE'};
ADE1.D(2:end,2)=num2cell(RMSE');
ADE1.D(1,3)={'r2 - Nash Sutcliff eff.'};
ADE1.D(2:end,3)=num2cell(r2');
ADE1.D(1,4)={'nRMSE'};
ADE1.D(2:end,4)=num2cell(nRMSE');
ADE1.D(1,5)={'logRMSE'};
ADE1.D(2:end,5)=num2cell(logRMSE');
ADE1.D(1,6)={'logr2'};
ADE1.D(2:end,6)=num2cell(logr2');
ADE1.D(1,7)={'Pearson_r2'};
ADE1.D(2:end,7)=num2cell(Pearson_r2');
ADE1.D(1,8)={'Pearson_logr2'};
ADE1.D(2:end,8)=num2cell(Pearson_logr2');
ADE1.D(1,9)={'Kling-Gupta eff.'};
ADE1.D(2:end,9)=num2cell(KGE');

% Regardless all the obj function, pick the best D for every objective function

ADE1.best_BTC(:,1)=BTC_input(:,1);        % time series 

[val, idx] = find(RMSE==min(RMSE(1,:)));
ADE1.Best(1,1)={'D'};                   % Best-fitting D [m^2/s] for min RMSE
ADE1.Best(2,1)=num2cell(B(idx,1));
ADE1.Best(1,2)={'RMSE'};                % min RMSE (the lower the better)
ADE1.Best(2,2)=num2cell(RMSE(val,idx));
ADE1.best_BTC(:,2)=CC1(:,idx);          % Save the conc. series for the best-fitting (RMSE) 
                                        % ADE1 curve

[val, idx] = find(r2==max(r2(1,:)));
ADE1.Best(1,3)={'D'};                   % Best-fitting D [m^2/s] for max r^2 (Nash-Sutcliff Eff)
ADE1.Best(2,3)=num2cell(B(idx,1));
ADE1.Best(1,4)={'r2 - Nash Sutcliff eff.'};                  % max r2 (the closer to 1 the better)
ADE1.Best(2,4)=num2cell(r2(val,idx));
ADE1.best_BTC(:,3)=CC1(:,idx);            % conc series of the best-fitting ADE1 curve

[val, idx] = find(nRMSE==min(nRMSE(1,:)));
ADE1.Best(1,5)={'D'};                      % Best-fitting D [m^2/s] for min nRMSE
ADE1.Best(2,5)=num2cell(B(idx,1));
ADE1.Best(1,6)={'nRMSE'};                  % min nRMSE (the lower the better)
ADE1.Best(2,6)=num2cell(nRMSE(val,idx));
ADE1.best_BTC(:,4)=CC1(:,idx);            % conc series of the best-fitting ADE1 curve

[val, idx] = find(logRMSE==min(logRMSE(1,:)));
ADE1.Best(1,7)={'D'};                        % Best-fitting D [m^2/s] for log RMSE
ADE1.Best(2,7)=num2cell(B(idx,1));
ADE1.Best(1,8)={'logRMSE'};                  % min logRMSE (the lower the better)
ADE1.Best(2,8)=num2cell(logRMSE(val,idx));
ADE1.best_BTC(:,5)=CC1(:,idx);            % conc series of the best-fitting ADE1 curve

[val, idx] = find(logr2==max(logr2(1,:)));
ADE1.Best(1,9)={'D'};                        % Best-fitting D [m^2/s] for log (r2 (or NSE))
ADE1.Best(2,9)=num2cell(B(idx,1));
ADE1.Best(1,10)={'logr2'};                    % max logr2 (the closer to 1 the better)
ADE1.Best(2,10)=num2cell(logr2(val,idx));
ADE1.best_BTC(:,6)=CC1(:,idx);            % conc series of the best-fitting ADE1 curve

[val, idx] = find(Pearson_r2==max(Pearson_r2(1,:)));
ADE1.Best(1,11)={'D'};                        % Best-fitting D [m^2/s] for Pearson
ADE1.Best(2,11)=num2cell(B(idx,1));
ADE1.Best(1,12)={'Pearson_r2'};               % max Pearson_r2 (the closer to 1 the better)
ADE1.Best(2,12)=num2cell(Pearson_r2(val,idx));
ADE1.best_BTC(:,7)=CC1(:,idx);            % conc series of the best-fitting ADE1 curve

[val, idx] = find(Pearson_logr2==max(Pearson_logr2(1,:)));
ADE1.Best(1,13)={'D'};                        % Best-fitting D [m^2/s] for log_Pearson
ADE1.Best(2,13)=num2cell(B(idx,1));
ADE1.Best(1,14)={'Pearson_logr2'};                    % max Pearson_logr2 (the closer to 1 the better)
ADE1.Best(2,14)=num2cell(Pearson_logr2(val,idx));
ADE1.best_BTC(:,8)=CC1(:,idx);            % conc series of the best-fitting ADE1 curve

[val, idx] = find(KGE==max(KGE(1,:)));
ADE1.Best(1,15)={'D'};                                % Best-fitting D [m^2/s] for ^max KGE
ADE1.Best(2,15)=num2cell(B(idx,1));
ADE1.Best(1,16)={'Kling-Gupta eff.'};                 % max KGE (the closer to 1 the better)
ADE1.Best(2,16)=num2cell(KGE(val,idx));
ADE1.best_BTC(:,9)=CC1(:,idx);            % conc series of the best-fitting ADE1 curve


%%%%%%%%%%%%%%%%

clear val idx C_fin1 C_fin2 D_fin1 D_fin2 Numerator j i k ssd RMSE r2 r2_1 A CC1 Diff2
clear B B_length D i k RMSE_1 RMSE r2 nRMSE logRMSE logr2 Pearson_r2 Pearson_logr2 KGE

%%%%%%%%%%%%%%%%
 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                             _____  ______   __                      % %
% %                       /\   |  __ \|  ____| |__ \                    % %
% %                      /  \  | |  | | |__       ) |                   % %
% %                     / /\ \ | |  | |  __|     / /                    % %
% %                    / ____ \| |__| | |____   / /_                    % %
% %                   /_/    \_\_____/|______| |____|                   % %
% %                                                                     % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% v is fixed and equal to L/t_peak 
% Q is measured by the v-notch upstream 
% A is calculated as (Q/v) 
% D is calibrated
% the consequence of this method is that Mass recovered from the injection 
% is =/= mass injected and it will be equal to Q_measured*sum of concentration

A=(Qmeasured/1000)/v1;    % m^2  Area derived using Qmeasured

M_recovered=Qmeasured*3.6*S;    % mass recovered [grams]
M_lostPerc=(M_g-M_recovered)*100/M_g; % Percentage of mass lost 
% We need to have mass balance! So the mass HAS TO be equal to the sum of
% the concentration times the measured discharge!

% pre-allocating D
D=0.0001; % [m^2/s] advection-dispersion coefficient

% Deciding an interval of calibration fod D -> again 0.0001:0.0001:2
B=[0.0001:0.0001:2];
B=B.';
B_length=length(B);
k=1;

% ADE using Mass recovered!!

for i=1:1:B_length
    D=B(i,1);
    for j=1:1:Length
    CC1(k,i)=(M_recovered/((A*(4*3.14*D*time(j,1))^(1/2))))*...
        exp(-((L-v1*time(j,1))^2/(4*D*time(j,1))));     % g - m - s -> Conc in g/m3
    k=k+1;
    end
    clear k
    k=1;
end

clear k i 

%%%%% Performances:
 
[RMSE,r2,nRMSE,logRMSE,logr2,Pearson_r2,Pearson_logr2,KGE] = ObjFun(BTC_input,CC1);

% Organize the results

ADE2.InjectedMass=M_g;      % injected Cl mass;
ADE2.RecoveredMass=M_recovered; % mass recovered;
ADE2.MassLostPerc=M_lostPerc; % Percentage of mass lost;
ADE2.L=L;           % Length of the stream reach [m]
ADE2.Qmeasured=Qmeasured;  % Discharge l/s;
ADE2.A=A;           % Wetted area A [m^2]
ADE2.v=v1;          % velocity [m/s]


% Preallocating the ADE2.D
ADE2.D(1:(length(B(:,1))+1),1:8)={'a'};     % all the obj. fun.

% list of all the D-values used with the relative obj functions
ADE2.D(1,1)={'D'};          
ADE2.D(2:end,1)=num2cell(B);
ADE2.D(1,2)={'RMSE'};
ADE2.D(2:end,2)=num2cell(RMSE');
ADE2.D(1,3)={'r2 - Nash Sutcliff eff.'};
ADE2.D(2:end,3)=num2cell(r2');
ADE2.D(1,4)={'nRMSE'};
ADE2.D(2:end,4)=num2cell(nRMSE');
ADE2.D(1,5)={'logRMSE'};
ADE2.D(2:end,5)=num2cell(logRMSE');
ADE2.D(1,6)={'logr2'};
ADE2.D(2:end,6)=num2cell(logr2');
ADE2.D(1,7)={'Pearson_r2'};
ADE2.D(2:end,7)=num2cell(Pearson_r2');
ADE2.D(1,8)={'Pearson_logr2'};
ADE2.D(2:end,8)=num2cell(Pearson_logr2');
ADE2.D(1,9)={'Kling-Gupta eff.'};
ADE2.D(2:end,9)=num2cell(KGE');

% regardless all the obj function, pick the best D for every objective function
ADE2.best_BTC(:,1)=BTC_input(:,1);        % time series of the best-fitting ADE2 curve 
    
[val, idx] = find(RMSE==min(RMSE(1,:)));
ADE2.Best(1,1)={'D'};                   % Best-fitting D [m^2/s] for min RMSE
ADE2.Best(2,1)=num2cell(B(idx,1));
ADE2.Best(1,2)={'RMSE'};                % min RMSE (the lower the better)
ADE2.Best(2,2)=num2cell(RMSE(val,idx));
ADE2.best_BTC(:,2)=CC1(:,idx);          % conc series of the best-fitting ADE2 curve 

[val, idx] = find(r2==max(r2(1,:)));
ADE2.Best(1,3)={'D'};                   % Best-fitting D [m^2/s] for max r^2 (Nash-Sutcliff Eff)
ADE2.Best(2,3)=num2cell(B(idx,1));
ADE2.Best(1,4)={'r2 - Nash Sutcliff eff.'};                  % max r2 (the closer to 1 the better)
ADE2.Best(2,4)=num2cell(r2(val,idx));
ADE2.best_BTC(:,3)=CC1(:,idx);            % conc series of the best-fitting ADE2 curve 

[val, idx] = find(nRMSE==min(nRMSE(1,:)));
ADE2.Best(1,5)={'D'};                      % Best-fitting D [m^2/s] for min nRMSE
ADE2.Best(2,5)=num2cell(B(idx,1));
ADE2.Best(1,6)={'nRMSE'};                  % min nRMSE (the lower the better)
ADE2.Best(2,6)=num2cell(nRMSE(val,idx));
ADE2.best_BTC(:,4)=CC1(:,idx);            % conc series of the best-fitting ADE2 curve 

[val, idx] = find(logRMSE==min(logRMSE(1,:)));
ADE2.Best(1,7)={'D'};                        % Best-fitting D [m^2/s] for log RMSE
ADE2.Best(2,7)=num2cell(B(idx,1));
ADE2.Best(1,8)={'logRMSE'};                  % min logRMSE (the lower the better)
ADE2.Best(2,8)=num2cell(logRMSE(val,idx));
ADE2.best_BTC(:,5)=CC1(:,idx);            % conc series of the best-fitting ADE2 curve 

[val, idx] = find(logr2==max(logr2(1,:)));
ADE2.Best(1,9)={'D'};                        % Best-fitting D [m^2/s] for log (r2 (or NSE))
ADE2.Best(2,9)=num2cell(B(idx,1));
ADE2.Best(1,10)={'logr2'};                    % max logr2 (the closer to 1 the better)
ADE2.Best(2,10)=num2cell(logr2(val,idx));
ADE2.best_BTC(:,6)=CC1(:,idx);            % conc series of the best-fitting ADE2 curve 

[val, idx] = find(Pearson_r2==max(Pearson_r2(1,:)));
ADE2.Best(1,11)={'D'};                        % Best-fitting D [m^2/s] for Pearson
ADE2.Best(2,11)=num2cell(B(idx,1));
ADE2.Best(1,12)={'Pearson_r2'};               % max Pearson_r2 (the closer to 1 the better)
ADE2.Best(2,12)=num2cell(Pearson_r2(val,idx));
ADE2.best_BTC(:,7)=CC1(:,idx);            % conc series of the best-fitting ADE2 curve 

[val, idx] = find(Pearson_logr2==max(Pearson_logr2(1,:)));
ADE2.Best(1,13)={'D'};                        % Best-fitting D [m^2/s] for log_Pearson
ADE2.Best(2,13)=num2cell(B(idx,1));
ADE2.Best(1,14)={'Pearson_logr2'};                    % max Pearson_logr2 (the closer to 1 the better)
ADE2.Best(2,14)=num2cell(Pearson_logr2(val,idx));
ADE2.best_BTC(:,8)=CC1(:,idx);            % conc series of the best-fitting ADE2 curve 

[val, idx] = find(KGE==max(KGE(1,:)));
ADE2.Best(1,15)={'D'};                                % Best-fitting D [m^2/s] for ^max KGE
ADE2.Best(2,15)=num2cell(B(idx,1));
ADE2.Best(1,16)={'Kling-Gupta eff.'};                 % max KGE (the closer to 1 the better)
ADE2.Best(2,16)=num2cell(KGE(val,idx));
ADE2.best_BTC(:,9)=CC1(:,idx);            % conc series of the best-fitting ADE2 curve 


clear val idx C_fin1 C_fin2 D_fin1 D_fin2 j i k ssd RMSE r2 r2_1 A CC1 RMSE_1 v1 D
clear TOP1 TOP01 TOP10 X_normalized X_scaled Lim01 Lim1 Lim10
clear Interp_curve1 Interp_curve01 Interp_curve10 Interp_curve_Hyper2 Interp_curve_Latin
clear k i LatinSample max_ranges_p min_ranges_p OFFSET offset SLOPE slope str
clear val idx C_fin1 C_fin2 D_fin1 D_fin2 Numerator j i k ssd RMSE r2 r2_1 A CC1 Diff2
clear B B_length D i k RMSE_1 RMSE r2 nRMSE logRMSE logr2 Pearson_r2 Pearson_logr2 KGE

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %                             _____  ______   ____                    % %
% %                       /\   |  __ \|  ____| |___ \                   % % 
% %                      /  \  | |  | | |__      __) |                  % %
% %                     / /\ \ | |  | |  __|    |__ <                   % %
% %                    / ____ \| |__| | |____   ___) |                  % %
% %                   /_/    \_\_____/|______| |____/                   % %
% %                                                                     % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% v is calibrated; A is calibrated; D is calibrated
% this option uses Monte Carlo or Latin Hypercube sampling and considers as reliable
% boundaries the one above a certain threshold for every objective function  

% PART 1 -> First Monte Carlo
% First parameter -> v
% lower limit: 20% less than the value obtained from v = L/t_peak

v1=ADE1.v;
v1=round(v1,4);             % Use the m/s
vmin=round(v1*0.8,4);       % Define the lower limit

vstep=0.0001;               % Define the step -> Useful only for Monte Carlo

% NEW 
% Define max flow velocity as calculated by the first detection time ...
% and linear distance between injection and detection points (Zhang et al.,
% 2020 - https://www.mdpi.com/1660-4601/17/19/7219) 

for i=1:1:length(BTC_input(:,1))
    if BTC_input(i,2)>0
        iii=i;
        break
    end
end
vmax=L/BTC_input(iii,1); % m/h
vmax=vmax/3600;
% END NEW 

% Second parameter -> A

A=round(ADE1.A,4);           % Use the m^2
Amin=A-0.2*A;                % Define the lower limit (-20%)
Amax=A*0.2+A;                % Define the upper limit (+20%)
Astep=0.0001;              % Define the step -> Only for Monte Carlo

% Third parameter -> D
D=cell2mat(ADE1.Best(2,1)); % Use in m^2/s
Dmin=0.0001;                % Define lower limit (my suggestion -> set it much lower than 0.8*D)
Dmax=D*0.2+D;               % Define upper limit (20% higher from the best D (RMSE) in ADE1)
Dstep=0.0001;             % An accuracy of 0.0001 -> Only for Monte Carlo
clear v1 A D B B_length

% The definition of the step (accuracy) is useful only for the Montecarlo,
% the script I wrote for the latin Hypercube just needs upper and lower
% limits, BUT it requires "lhsdesign" function, so an installed toolbox with
% this function (eg. Statistics and Machine Learning Toolbox)

[ADE3,ADE3_limits,Hyperspace]=MonteCarlo(n, vmin,vmax,vstep,Amin,Amax,Astep,Dmin,Dmax,...
    Dstep,M_g,time,L,BTC_input,ADE1, ADE2,Description1,Description2,Length,Q1);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

end

