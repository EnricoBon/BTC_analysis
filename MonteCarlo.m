function [ADE3,ADE3_limits,Hyperspace] = MonteCarlo(n, vmin,vmax,vstep,Amin,Amax,Astep,Dmin,...
    Dmax,Dstep,M_g,time,L,BTC_input, ADE1, ADE2,Description1,Description2, Length,Q1)

% NOTE: to run properly the Latin Hypercube sampling the MATLAB version
% should have installed a package with "lhsdesign" function (eg. Statistics and
% Machine Learning Toolbox).

% OPTION 1 -> UNIFORMALLY DISTRIBUTED SAMPLES VIA MONTE CARLO SAMPLING
% If you choose this method please note that you need the Dstep, vstep, and Astep
% that I've originally removed at the end of ADE_analysis, before calling
% this (MonteCarlo) function

% D=randsample([Dmin:Dstep:Dmax],n,true);         % replacement are allowed 
% D=D';
% A=randsample([Amin:Astep:Amax],n,true);
% A=A';
% v=randsample([vmin:vstep:vmax],n,true);
% v=v';
%%%%% END OPTION 1

% OPTION 2 -> BUILD A LATIN HYPERCUBE SAMPLING
% Build 2 vectors with 1 column -> each row of the column has the min and
% max values for the parameters
min_ranges= [vmin; Amin; Dmin];     % Min values acceptable for v, A and D
max_ranges= [vmax; Amax; Dmax];     % Max values acceptable for v, A and D
p=length(min_ranges);     

slope=max_ranges-min_ranges;
offset=min_ranges;
SLOPE=ones(n,p);
OFFSET=ones(n,p);
for i=1:p
    SLOPE(:,i)=ones(n,1).*slope(i);
    OFFSET(:,i)=ones(n,1).*offset(i);
end
X_normalized = lhsdesign(n,p);        % this one builds latin hypercube in [0 1] interval
X_scaled=SLOPE.*X_normalized+OFFSET;  % this one modify the [0 1] latin hypercube in the
                                      % latin hypercube of our parameters
v=X_scaled(:,1);
A=X_scaled(:,2);
D=X_scaled(:,3);
%%%%% END OPTION 2

%%%% NOTE -> take all the combinations or not...?
% Choice for the modeller: the random v * the random A has to be between the
% reliable values of discharge: +- 10% of error for Dilution gauging 
% remember that Q1=l/s
Q_ref=Q1/1000;  % m^3/s
Q_ref_min=Q_ref*0.9;
Q_ref_max=Q_ref*1.1;
% +-10% of Q from dilution gauging method coming from: Schmadel, N. M., 
% Neilson, B. T., and Stevens, D. K.: Approaches to estimate uncertainty in
% longitudinal channel water balances, J. Hydrol., 394, 357–369, 2010
%%%% END of the note
clear Input

Input(:,1)=v; % m/s
Input(:,2)=A; % m^2
Input(:,3)=D; % m^2/s
Input(:,4)=Input(:,1).*Input(:,2);

clear i
% 
k=1;
wbhandle=waitbar(0,'Creating Latin Hypercube set for ADE - Removing unrealistic discharge values.');
 for i=1:1:n
    if Input(i,4)>Q_ref_min && Input(i,4)<Q_ref_max
	Input_fin(k,:)=Input(i,:); 
	k=k+1;
    end
waitbar(i / n)
 end

clear v A D Input
Input=Input_fin;
clear Input_fin

N=length(Input(:,1)); % new length with removed unrealistic values
v(1:N,1)=Input(1:N,1);
A(1:N,1)=Input(1:N,2);
D(1:N,1)=Input(1:N,3);

close(wbhandle)
% Note for the modeller -> consider that a great amount of parameter sets can
% be removed by setting this discharge condition. Therefore, if you need a 
% NET number of simulation = n consider setting n=5*n at the beginning of 
% the latin hypercube for the ADE.
% However, due to the strong identifiability of A, D and v in the ADE even
% 10000 simulation will give a strong identifiability for the parameters.

% Build the string for the figures -> so we know the extremes of the
% sampling immediately once we see the first figure

formatSpec1="v_m_i_n=%g m/s";
formatSpec2="v_m_a_x=%g m/s";
formatSpec3="A_m_i_n=%g m^2";
formatSpec4="A_m_a_x=%g m^2";
formatSpec5="D_m_i_n=%g m^2/s"; 
formatSpec6="D_m_a_x=%g m^2/s";
str(1,1)=sprintf(formatSpec1,vmin);
str(2,1)=sprintf(formatSpec2,vmax);
str(3,1)=sprintf(formatSpec3,Amin);
str(4,1)=sprintf(formatSpec4,Amax);
str(5,1)=sprintf(formatSpec5,Dmin);
str(6,1)=sprintf(formatSpec6,Dmax);

%%%% Let's clear the workspace a bit...
clear X_scaled X_normalized offset min_ranges_p SLOPE OFFSET max_ranges_p p
clear Dmin Dstep Dmax Amin Astep Amax vmin vstep vmax ssd RMSE i k
% values we need for the ADE: Mass in grams, length of the stream reach
% let's go with our ADE model run

clear Hyperspace
wbhandle=waitbar(0,'Completing LatinHypercube ADE simulations. Please wait.');

for k=1:1:N
    for i=1:1:Length 
        CC1(i,1)=(M_g/((A(k,1)*(4*3.14*D(k,1)*time(i,1))^(1/2))))*...   % g - m - s -> CC1 in g/m3
            exp(-((L-v(k,1)*time(i,1))^2/(4*D(k,1)*time(i,1))));
    end
    clear i

% % % For every BTC let's calculate the fitting with the ObjFunc code 
[RMSE,r2,nRMSE,logRMSE,logr2,Pearson_r2,Pearson_logr2,KGE] = ObjFun(BTC_input,CC1);
% Store the results in the Hyperspace
    Hyperspace(k,1)=v(k,1);
    Hyperspace(k,2)=A(k,1);
    Hyperspace(k,3)=D(k,1);
    Hyperspace(k,4)=RMSE(1,1);
    Hyperspace(k,5)=r2(1,1);
    Hyperspace(k,6)=nRMSE(1,1);
    Hyperspace(k,7)=logRMSE(1,1);
    Hyperspace(k,8)=logr2(1,1);
    Hyperspace(k,9)=Pearson_r2(1,1);
    Hyperspace(k,10)=Pearson_logr2(1,1);
    Hyperspace(k,11)=KGE(1,1);
    
ADE_temp(:,k)= CC1(:,1);
clear CC1 RMSE r2 nRMSE logRMSE logr2 Pearson_r2 Pearson_logr2 KGE
waitbar(k / N)
end
close(wbhandle)

% % % Our Hyperspace matrix contain all the MonteCarlo/Latin Hypercube runs 
% with the parameter values, their efficiency in term of RMSE, r2, etc...
% Let's pick some thresholds for a better representation
N=length(Hyperspace(:,1));   % set the limit for the actual simulation considered

top20=round(N*0.2);        % Limit for top 20% of the results
top10=round(N*0.1);        % Limit for top 10% of the results
top1=round(N*0.01);        % Limit for top 1% of the results
top01=round(N*0.001);      % Limit for top 0.1% of the results
Hyperspace_Order=sortrows(Hyperspace,4);   % Hyperspace sorted depending on RMSE value

for i=1:1:length(Hyperspace(:,1))
[val,idx]=find(Hyperspace_Order(i,4)==Hyperspace(:,4));
ORDER(i,1)=val;
end
% Order ADEs accordingly
for i=1:1:length(Hyperspace(:,1))
ADE_temp_Order(:,i)=ADE_temp(:,ORDER(i,1));
end

% Evaluate the curve where, for each time value, we have the max and min
% concentration among all the possible ADE, so we can represent the grey 
% area of the possibilities without affecting the plotting computational time
clear i

for i=1:1:Length
    Interp_curve_All(i,1)=time(i,1);
    Interp_curve_All(i,2)=min(ADE_temp_Order(i,:));
    Interp_curve_All(i,3)=max(ADE_temp_Order(i,:));
    
    Interp_curve_20(i,1)=time(i,1);
    Interp_curve_20(i,2)=min(ADE_temp_Order(i,1:top20));
    Interp_curve_20(i,3)=max(ADE_temp_Order(i,1:top20));
    
    Interp_curve_10(i,1)=time(i,1);
    Interp_curve_10(i,2)=min(ADE_temp_Order(i,1:top10));
    Interp_curve_10(i,3)=max(ADE_temp_Order(i,1:top10));
    
    Interp_curve_1(i,1)=time(i,1);
    Interp_curve_1(i,2)=min(ADE_temp_Order(i,1:top1));
    Interp_curve_1(i,3)=max(ADE_temp_Order(i,1:top1));
    
    Interp_curve_01(i,1)=time(i,1);
    Interp_curve_01(i,2)=min(ADE_temp_Order(i,1:top01));
    Interp_curve_01(i,3)=max(ADE_temp_Order(i,1:top01));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %  ______ _      __            _    _                                 _          
% % |  ____(_)    /_ |          | |  | |                               | |         
% % | |__   _  __ _| |  ______  | |__| |_   _ _ __   ___ _ __ ___ _   _| |__   ___ 
% % |  __| | |/ _` | | |______| |  __  | | | | '_ \ / _ \ '__/ __| | | | '_ \ / _ \
% % | |    | | (_| | |          | |  | | |_| | |_) |  __/ | | (__| |_| | |_) |  __/
% % |_|    |_|\__, |_|          |_|  |_|\__, | .__/ \___|_|  \___|\__,_|_.__/ \___|
% %            __/ |                     __/ | |                                   
% %           |___/                     |___/|_|                                   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% FIGURE #1 - All the results of the hypersphere
% In this figure we show all the results relative to one objective function
% (RMSE) and we plot the RMSE vs the three variables (v, A and D). Since
% the Hyperspace is made by only 3 variables we can easily represent it as 
% a cube (each variable for each dimension of the cube).

f=figure;
subplot(2,4,1)
plot(Hyperspace(:,1),Hyperspace(:,4),'.k')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Hyperspace(:,1)) max(Hyperspace(:,1))])
ylim([0 max(Hyperspace(:,4))])

subplot(2,4,2)
plot(Hyperspace(:,2),Hyperspace(:,4),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
xlim([min(Hyperspace(:,2)) max(Hyperspace(:,2))])
ylim([0 max(Hyperspace(:,4))])

subplot(2,4,5)
plot(Hyperspace(:,3),Hyperspace(:,4),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Hyperspace(:,3)) max(Hyperspace(:,3))])
ylim([0 max(Hyperspace(:,4))])

subplot(2,4,6)
plot(Interp_curve_All(:,1), Interp_curve_All(:,2), 'k', 'LineWidth', 1,'HandleVisibility','off');
hold on;
plot(Interp_curve_All(:,1), Interp_curve_All(:,3), 'k', 'LineWidth', 1,'HandleVisibility','off');
x2 = [Interp_curve_All(:,1)', fliplr(Interp_curve_All(:,1)')];
inBetween = [Interp_curve_All(:,2)', fliplr(Interp_curve_All(:,3)')];
fill(x2, inBetween, [0.86,0.86,0.86]);
hold on
plot(time(:,1),BTC_input(:,2),'-r','LineWidth',2)
xlabel ('t [s]');
ylabel ('Conc [mg/l]');
legend('Simulation','Observed BTC')
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String', "Latin Hypercube sampling",...
    'FontSize',12,'LineStyle','none','FitBoxToText','on')
annotation('textbox', [0.01, 0.8, 0.1, 0.1], 'String',str,'FontSize',11,'LineStyle','none')

subplot(2,4,[3,4,7,8])
plot3(Hyperspace(:,1),Hyperspace(:,2),Hyperspace(:,3),'.k')
hold on
plot3(Hyperspace_Order(1:top10,1),Hyperspace_Order(1:top10,2),Hyperspace_Order(1:top10,3),'.r','MarkerSize',10)
hold on
plot3(Hyperspace_Order(1:top1,1),Hyperspace_Order(1:top1,2),Hyperspace_Order(1:top1,3),'.g','MarkerSize',10)
legend('All Simulations','Top 10%','Top 1%')
legend('Location','northeast')
xlim([min(Hyperspace(:,1)) max(Hyperspace(:,1))])
ylim([min(Hyperspace(:,2)) max(Hyperspace(:,2))])
zlim([min(Hyperspace(:,3)) max(Hyperspace(:,3))])
xlabel ('v [m/s]');
ylabel ('A [m^2]');
zlabel ('D [m^2/s]');

sgtitle({Description1;Description2},'FontSize',14);
set(gcf, 'WindowState', 'maximized');
saveas(f,[pwd '/Output_files/Hypercube.fig']);
saveas(f,[pwd '/Output_files/Hypercube.tif']);
close(f)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %  ______ _       ___             _______            _____  _       _       
% % |  ____(_)     |__ \           |__   __|          |  __ \| |     | |      
% % | |__   _  __ _   ) |  ______     | | ___  _ __   | |__) | | ___ | |_ ___ 
% % |  __| | |/ _` | / /  |______|    | |/ _ \| '_ \  |  ___/| |/ _ \| __/ __|
% % | |    | | (_| |/ /_              | | (_) | |_) | | |    | | (_) | |_\__ \
% % |_|    |_|\__, |____|             |_|\___/| .__/  |_|    |_|\___/ \__|___/
% %            __/ |                          | |                             
% %           |___/                           |_|                             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
% FIGURE #2 - We plot only the top 20 - 10 - and 1% of the results of the 
% hypersphere with the performances of the 3 parameters (v, A and D) and
% the simulations (in grey area) of tracer-solute concentrations
% corresponding to the set of v, A and D

g=figure;
subplot(3,4,1)
plot(Hyperspace_Order(1:top20,1),Hyperspace_Order(1:top20,4),'.k')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Hyperspace_Order(1:top20,1)) max(Hyperspace_Order(1:top20,1))])
ylim([0 max(Hyperspace_Order(1:top20,4))])

subplot(3,4,2)
plot(Hyperspace_Order(1:top20,2),Hyperspace_Order(1:top20,4),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
xlim([min(Hyperspace_Order(1:top20,2)) max(Hyperspace_Order(1:top20,2))])
ylim([0 max(Hyperspace_Order(1:top20,4))])

subplot(3,4,3)
plot(Hyperspace_Order(1:top20,3),Hyperspace_Order(1:top20,4),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Hyperspace_Order(1:top20,3)) max(Hyperspace_Order(1:top20,3))])
ylim([0 max(Hyperspace_Order(1:top20,4))])

subplot(3,4,4)
plot(Interp_curve_20(:,1), Interp_curve_20(:,2), 'k', 'LineWidth', 1,'HandleVisibility','off');
hold on;
plot(Interp_curve_20(:,1), Interp_curve_20(:,3), 'k', 'LineWidth', 1,'HandleVisibility','off');
x2 = [Interp_curve_20(:,1)', fliplr(Interp_curve_20(:,1)')];
inBetween = [Interp_curve_20(:,2)', fliplr(Interp_curve_20(:,3)')];
fill(x2, inBetween, [0.86,0.86,0.86]);
hold on
plot(time(:,1),BTC_input(:,2),'-r','LineWidth',2)
xlabel ('t [s]');
ylabel ('Conc [mg/l]');
legend('Simulation','Observed BTC')
annotation('textbox', [0.12, 0.85, 0.1, 0.1], 'String', "Latin Hypercube sampling - top 20%",...
    'FontSize',12,'LineStyle','none','FitBoxToText','on')
annotation('textbox', [0.01, 0.8, 0.1, 0.1], 'String',str,'FontSize',11,'LineStyle','none')
set(gcf, 'WindowState', 'maximized');

subplot(3,4,5)
plot(Hyperspace_Order(1:top10,1),Hyperspace_Order(1:top10,4),'.k')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Hyperspace_Order(1:top10,1)) max(Hyperspace_Order(10:top1,1))])
ylim([0 max(Hyperspace_Order(1:top10,4))])

subplot(3,4,6)
plot(Hyperspace_Order(1:top10,2),Hyperspace_Order(1:top10,4),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
xlim([min(Hyperspace_Order(1:top10,2)) max(Hyperspace_Order(1:top10,2))])
ylim([0 max(Hyperspace_Order(1:top10,4))])

subplot(3,4,7)
plot(Hyperspace_Order(1:top10,3),Hyperspace_Order(1:top10,4),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Hyperspace_Order(1:top10,3)) max(Hyperspace_Order(1:top10,3))])
ylim([0 max(Hyperspace_Order(1:top10,4))])

subplot(3,4,8)
plot(Interp_curve_10(:,1), Interp_curve_10(:,2), 'k', 'LineWidth', 1,'HandleVisibility','off');
hold on;
plot(Interp_curve_10(:,1), Interp_curve_10(:,3), 'k', 'LineWidth', 1,'HandleVisibility','off');
x2 = [Interp_curve_10(:,1)', fliplr(Interp_curve_10(:,1)')];
inBetween = [Interp_curve_10(:,2)', fliplr(Interp_curve_10(:,3)')];
fill(x2, inBetween, [0.86,0.86,0.86]);
hold on
plot(time(:,1),BTC_input(:,2),'-r','LineWidth',2)
xlabel ('t [s]');
ylabel ('Conc [mg/l]');
legend('Simulation','Observed BTC')
annotation('textbox', [0.12, 0.55, 0.1, 0.1], 'String', "Latin Hypercube sampling - top 10%",...
    'FontSize',12,'LineStyle','none','FitBoxToText','on')
set(gcf, 'WindowState', 'maximized');

subplot(3,4,9)
plot(Hyperspace_Order(1:top1,1),Hyperspace_Order(1:top1,4),'.k')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Hyperspace_Order(1:top1,1)) max(Hyperspace_Order(1:top1,1))])
ylim([0 max(Hyperspace_Order(1:top1,4))])

subplot(3,4,10)
plot(Hyperspace_Order(1:top1,2),Hyperspace_Order(1:top1,4),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
xlim([min(Hyperspace_Order(1:top1,2)) max(Hyperspace_Order(1:top1,2))])
ylim([0 max(Hyperspace_Order(1:top1,4))])

subplot(3,4,11)
plot(Hyperspace_Order(1:top1,3),Hyperspace_Order(1:top1,4),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Hyperspace_Order(1:top1,3)) max(Hyperspace_Order(1:top1,3))])
ylim([0 max(Hyperspace_Order(1:top1,4))])

subplot(3,4,12)
plot(Interp_curve_1(:,1), Interp_curve_1(:,2), 'k', 'LineWidth', 1,'HandleVisibility','off');
hold on;
plot(Interp_curve_1(:,1), Interp_curve_1(:,3), 'k', 'LineWidth', 1,'HandleVisibility','off');
x2 = [Interp_curve_1(:,1)', fliplr(Interp_curve_1(:,1)')];
inBetween = [Interp_curve_1(:,2)', fliplr(Interp_curve_1(:,3)')];
fill(x2, inBetween, [0.86,0.86,0.86]);
hold on
plot(time(:,1),BTC_input(:,2),'-r','LineWidth',2)
xlabel ('t [s]');
ylabel ('Conc [mg/l]');
legend('Simulation','Observed BTC')
annotation('textbox', [0.12, 0.255, 0.1, 0.1], 'String', "Latin Hypercube sampling - top 1%",...
    'FontSize',12,'LineStyle','none','FitBoxToText','on')
set(gcf, 'WindowState', 'maximized');

sgtitle({Description1,Description2},'FontSize',14);

saveas(g,[pwd '/Output_files/Hypercube_topResults.fig']);
saveas(g,[pwd '/Output_files/Hypercube_topResults.tif']);
close (g);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % %                  _____ _        _      __ _       
% % %                 / ____| |      | |    / _(_)      
% % %                | (___ | |_ __ _| |_  | |_ _  __ _ 
% % %                 \___ \| __/ _` | __| |  _| |/ _` |
% % %                 ____) | || (_| | |_  | | | | (_| |
% % %                |_____/ \__\__,_|\__| |_| |_|\__, |
% % %                                              __/ |
% % %                                             |___/ 
% % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

[stat1, stat2,Range_param_RMSE,Likelihood_RMSE] = StatFigures(Hyperspace,Description1,Description2);

ADE3.Range_param_RMSE=Range_param_RMSE; % Save the range of top 20%-0.1% parameters
ADE3.Likelihood_RMSE=Likelihood_RMSE;   % Save the range parameters from 10 to 1%
% with a step of 1% where the parameters are listed with their likelyhood
% Likelihood is evaluated as:
% of=Hyperspace_best(:,4);       % Select objective function RMSE
% of=of./max(of);                % normalised 
% of=1-of;                       % likelihood (high values indicate more likely [probable] models)
% if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods
% According to Ward et al., 2017;
% So Likelihood_RMSE.v.ten will have the first column with all the velocity
% values in the 10th interval of Hyperspace ordered according to RMSE, and
% their associated likelihood in the second column

% stat 1 = figure containing v, A and D data vs RMSE (1st row);
% Frequency plots (2nd row); normalized cumulative distribution of the
% parameters (3rd row)
% stat 2 = figure containing v, A and D data vs RMSE (1st row);
% A posteriori parameter distribution; (2nd row); parameter
% identifiability (3rd row)

saveas(stat1,[pwd '/Output_files/Stat_Results1_RMSE.fig']);
saveas(stat1,[pwd '/Output_files/Stat_Results1_RMSE.tif']);
saveas(stat2,[pwd '/Output_files/Stat_Results2_RMSE.fig']);
saveas(stat2,[pwd '/Output_files/Stat_Results2_RMSE.tif']);

close (stat1, stat2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If wanted, it's possible to include StatFigures_var, which gives the
% statistic plots for r^2 (NSE), logRMSE and KGE. These plots might help the
% modeller to understand there is nosubstantial difference between using
% Nash-Sutcliff and RMSE (both of them give the same interval for
% performance thresholds), the less-efficiency of logRMSE; and the KGE that
% indicates large parameter ranges for equally good performances, and thus does not
% help to target identifiability.

[stat1_r2, stat2_r2, stat1_logRMSE, stat2_logRMSE,stat1_KGE,stat2_KGE,...
Range_param_r2,Range_param_logRMSE,Range_param_KGE,Range_param_nRMSE] = StatFigures_var(Hyperspace,Description1,Description2);

% Save the range of the top 20%-0.1% of the parameters for different Obj.
% Function
ADE3.Range_param_r2=Range_param_r2; % r2
ADE3.Range_param_logRMSE=Range_param_logRMSE; % logRMSE
ADE3.Range_param_KGE=Range_param_KGE; % KGE
ADE3.Range_param_nRMSE=Range_param_nRMSE; % nRMSE

% r^2 - NSE
saveas(stat1_r2,[pwd '/Output_files/Stat_Results1_r2.fig']);
saveas(stat1_r2,[pwd '/Output_files/Stat_Results1_r2.tif']);
saveas(stat2_r2,[pwd '/Output_files/Stat_Results2_r2.fig']);
saveas(stat2_r2,[pwd '/Output_files/Stat_Results2_r2.tif']);
% logRMSE
saveas(stat1_logRMSE,[pwd '/Output_files/Stat_Results1_logRMSE.fig']);
saveas(stat1_logRMSE,[pwd '/Output_files/Stat_Results1_logRMSE.tif']);
saveas(stat2_logRMSE,[pwd '/Output_files/Stat_Results2_logRMSE.fig']);
saveas(stat2_logRMSE,[pwd '/Output_files/Stat_Results2_logRMSE.tif']);
% KGE
saveas(stat1_KGE,[pwd '/Output_files/Stat_Results1_KGE.fig']);
saveas(stat1_KGE,[pwd '/Output_files/Stat_Results1_KGE.tif']);
saveas(stat2_KGE,[pwd '/Output_files/Stat_Results2_KGE.fig']);
saveas(stat2_KGE,[pwd '/Output_files/Stat_Results2_KGE.tif']);

close (stat1_r2, stat2_r2, stat1_logRMSE, stat2_logRMSE,stat1_KGE,stat2_KGE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Our analysis is terminated, we need to wrap up the results of the latin
% hypercube analysis into a unique matrix

ADE3.InjectedMass=M_g;      % injected Cl mass;
ADE3.L=L;                   % Length of the stream reach [m]
% preallocating the hyperspace
ADE3.Hyperspace(1:(length(Hyperspace(:,1))+1),1:11)={'a'};
% This struct contains all the results 
ADE3.Hyperspace(1,1)={'v'};          
ADE3.Hyperspace(2:end,1)=num2cell(Hyperspace(:,1));
ADE3.Hyperspace(1,2)={'A'};          
ADE3.Hyperspace(2:end,2)=num2cell(Hyperspace(:,2));
ADE3.Hyperspace(1,3)={'D'};          
ADE3.Hyperspace(2:end,3)=num2cell(Hyperspace(:,3));
ADE3.Hyperspace(1,4)={'RMSE'};          
ADE3.Hyperspace(2:end,4)=num2cell(Hyperspace(:,4));
ADE3.Hyperspace(1,5)={'r2 - Nash Sutcliff eff.'};          
ADE3.Hyperspace(2:end,5)=num2cell(Hyperspace(:,5));
ADE3.Hyperspace(1,6)={'nRMSE'};          
ADE3.Hyperspace(2:end,6)=num2cell(Hyperspace(:,6));
ADE3.Hyperspace(1,7)={'logRMSE'};          
ADE3.Hyperspace(2:end,7)=num2cell(Hyperspace(:,7));
ADE3.Hyperspace(1,8)={'logr2'};          
ADE3.Hyperspace(2:end,8)=num2cell(Hyperspace(:,8));
ADE3.Hyperspace(1,9)={'Pearson_r2'};          
ADE3.Hyperspace(2:end,9)=num2cell(Hyperspace(:,9));
ADE3.Hyperspace(1,10)={'Pearson_logr2'};          
ADE3.Hyperspace(2:end,10)=num2cell(Hyperspace(:,10));
ADE3.Hyperspace(1,11)={'Kling-Gupta eff.'};          
ADE3.Hyperspace(2:end,11)=num2cell(Hyperspace(:,11));

% It's worth saving the range of the top 10% of the results to give a range
% of the parameters and report the performances associated
Hyperspace_10=sortrows(Hyperspace,4);
% Take the 10%
Hyperspace_10=Hyperspace_10(1:top10,:);

% preallocating the Hyperspace_10
ADE3.Hyperspace_10(1:length(Hyperspace_10(:,1))+1,1:4)={'a'};
%
ADE3.Hyperspace_10(1,1)={'v'};          
ADE3.Hyperspace_10(2:end,1)=num2cell(Hyperspace_10(:,1));
ADE3.Hyperspace_10(1,2)={'A'};          
ADE3.Hyperspace_10(2:end,2)=num2cell(Hyperspace_10(:,2));
ADE3.Hyperspace_10(1,3)={'D'};          
ADE3.Hyperspace_10(2:end,3)=num2cell(Hyperspace_10(:,3));
ADE3.Hyperspace_10(1,4)={'RMSE'};          
ADE3.Hyperspace_10(2:end,4)=num2cell(Hyperspace_10(:,4));

% Parameter boundary and RMSE range
% preallocating the Hyperspace_10_Limits
ADE3.Hyperspace_10_Limits(1:2,1:10)={'a'};
%

ADE3.Hyperspace_10_Limits(1,1)={'v min'};          
ADE3.Hyperspace_10_Limits(2,1)=num2cell(min(Hyperspace_10(:,1)));
ADE3.Hyperspace_10_Limits(1,2)={'v max'};          
ADE3.Hyperspace_10_Limits(2,2)=num2cell(max(Hyperspace_10(:,1)));
ADE3.Hyperspace_10_Limits(1,3)={'A min'};          
ADE3.Hyperspace_10_Limits(2,3)=num2cell(min(Hyperspace_10(:,2)));
ADE3.Hyperspace_10_Limits(1,4)={'A max'};          
ADE3.Hyperspace_10_Limits(2,4)=num2cell(max(Hyperspace_10(:,2)));
ADE3.Hyperspace_10_Limits(1,5)={'D min'};          
ADE3.Hyperspace_10_Limits(2,5)=num2cell(min(Hyperspace_10(:,3)));
ADE3.Hyperspace_10_Limits(1,6)={'D max'};          
ADE3.Hyperspace_10_Limits(2,6)=num2cell(max(Hyperspace_10(:,3)));
ADE3.Hyperspace_10_Limits(1,7)={'RMSE min'};          
ADE3.Hyperspace_10_Limits(2,7)=num2cell(min(Hyperspace_10(:,4)));
ADE3.Hyperspace_10_Limits(1,8)={'RMSE max'};          
ADE3.Hyperspace_10_Limits(2,8)=num2cell(max(Hyperspace_10(:,4)));
ADE3.Hyperspace_10_Limits(1,9)={'r^2 (NSE) min'};          
ADE3.Hyperspace_10_Limits(2,9)=num2cell(min(Hyperspace_10(:,5)));
ADE3.Hyperspace_10_Limits(1,10)={'r^2 (NSE) max'};          
ADE3.Hyperspace_10_Limits(2,10)=num2cell(max(Hyperspace_10(:,5)))

% Now it's worth to consider just the top 0.1% of the result to give a range
% of the parameters v, A and D and report the performances associated ->
% RMSE chose as preferred objective function
Hyperspace_01=sortrows(Hyperspace,4);
% Take the 0.1%
Hyperspace_01=Hyperspace_01(1:top01,:);

% preallocating the Hyperspace_01
ADE3.Hyperspace_01(1:length(Hyperspace_01(:,1))+1,1:4)={'a'};
%
ADE3.Hyperspace_01(1,1)={'v'};          
ADE3.Hyperspace_01(2:end,1)=num2cell(Hyperspace_01(:,1));
ADE3.Hyperspace_01(1,2)={'A'};          
ADE3.Hyperspace_01(2:end,2)=num2cell(Hyperspace_01(:,2));
ADE3.Hyperspace_01(1,3)={'D'};          
ADE3.Hyperspace_01(2:end,3)=num2cell(Hyperspace_01(:,3));
ADE3.Hyperspace_01(1,4)={'RMSE'};          
ADE3.Hyperspace_01(2:end,4)=num2cell(Hyperspace_01(:,4));

% Parameter boundary and RMSE range
% preallocating the Hyperspace_01_Limits
ADE3.Hyperspace_01_Limits(1:2,1:10)={'a'};
%

ADE3.Hyperspace_01_Limits(1,1)={'v min'};          
ADE3.Hyperspace_01_Limits(2,1)=num2cell(min(Hyperspace_01(:,1)));
ADE3.Hyperspace_01_Limits(1,2)={'v max'};          
ADE3.Hyperspace_01_Limits(2,2)=num2cell(max(Hyperspace_01(:,1)));
ADE3.Hyperspace_01_Limits(1,3)={'A min'};          
ADE3.Hyperspace_01_Limits(2,3)=num2cell(min(Hyperspace_01(:,2)));
ADE3.Hyperspace_01_Limits(1,4)={'A max'};          
ADE3.Hyperspace_01_Limits(2,4)=num2cell(max(Hyperspace_01(:,2)));
ADE3.Hyperspace_01_Limits(1,5)={'D min'};          
ADE3.Hyperspace_01_Limits(2,5)=num2cell(min(Hyperspace_01(:,3)));
ADE3.Hyperspace_01_Limits(1,6)={'D max'};          
ADE3.Hyperspace_01_Limits(2,6)=num2cell(max(Hyperspace_01(:,3)));
ADE3.Hyperspace_01_Limits(1,7)={'RMSE min'};          
ADE3.Hyperspace_01_Limits(2,7)=num2cell(min(Hyperspace_01(:,4)));
ADE3.Hyperspace_01_Limits(1,8)={'RMSE max'};          
ADE3.Hyperspace_01_Limits(2,8)=num2cell(max(Hyperspace_01(:,4)));
ADE3.Hyperspace_01_Limits(1,9)={'r^2 (NSE) min'};          
ADE3.Hyperspace_01_Limits(2,9)=num2cell(min(Hyperspace_01(:,5)));
ADE3.Hyperspace_01_Limits(1,10)={'r^2 (NSE) max'};          
ADE3.Hyperspace_01_Limits(2,10)=num2cell(max(Hyperspace_01(:,5)));

% It's worth saving the range of the top 20% of the results to give a range
% of the parameters and report the performances associated
Hyperspace_20=sortrows(Hyperspace,4);
% Take the 20%
Hyperspace_20=Hyperspace_20(1:top20,:);

% preallocating the Hyperspace_20
ADE3.Hyperspace_20(1:length(Hyperspace_20(:,1))+1,1:4)={'a'};
%
ADE3.Hyperspace_20(1,1)={'v'};          
ADE3.Hyperspace_20(2:end,1)=num2cell(Hyperspace_20(:,1));
ADE3.Hyperspace_20(1,2)={'A'};          
ADE3.Hyperspace_20(2:end,2)=num2cell(Hyperspace_20(:,2));
ADE3.Hyperspace_20(1,3)={'D'};          
ADE3.Hyperspace_20(2:end,3)=num2cell(Hyperspace_20(:,3));
ADE3.Hyperspace_20(1,4)={'RMSE'};          
ADE3.Hyperspace_20(2:end,4)=num2cell(Hyperspace_20(:,4));

% Parameter boundary and RMSE range
% preallocating the Hyperspace_20_Limits
ADE3.Hyperspace_20_Limits(1:2,1:20)={'a'};
%

ADE3.Hyperspace_20_Limits(1,1)={'v min'};          
ADE3.Hyperspace_20_Limits(2,1)=num2cell(min(Hyperspace_20(:,1)));
ADE3.Hyperspace_20_Limits(1,2)={'v max'};          
ADE3.Hyperspace_20_Limits(2,2)=num2cell(max(Hyperspace_20(:,1)));
ADE3.Hyperspace_20_Limits(1,3)={'A min'};          
ADE3.Hyperspace_20_Limits(2,3)=num2cell(min(Hyperspace_20(:,2)));
ADE3.Hyperspace_20_Limits(1,4)={'A max'};          
ADE3.Hyperspace_20_Limits(2,4)=num2cell(max(Hyperspace_20(:,2)));
ADE3.Hyperspace_20_Limits(1,5)={'D min'};          
ADE3.Hyperspace_20_Limits(2,5)=num2cell(min(Hyperspace_20(:,3)));
ADE3.Hyperspace_20_Limits(1,6)={'D max'};          
ADE3.Hyperspace_20_Limits(2,6)=num2cell(max(Hyperspace_20(:,3)));
ADE3.Hyperspace_20_Limits(1,7)={'RMSE min'};          
ADE3.Hyperspace_20_Limits(2,7)=num2cell(min(Hyperspace_20(:,4)));
ADE3.Hyperspace_20_Limits(1,8)={'RMSE max'};          
ADE3.Hyperspace_20_Limits(2,8)=num2cell(max(Hyperspace_20(:,4)));
ADE3.Hyperspace_20_Limits(1,9)={'r^2 (NSE) min'};          
ADE3.Hyperspace_20_Limits(2,9)=num2cell(min(Hyperspace_20(:,5)));
ADE3.Hyperspace_20_Limits(1,20)={'r^2 (NSE) max'};          
ADE3.Hyperspace_20_Limits(2,20)=num2cell(max(Hyperspace_20(:,5)));


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %                           _____  ______     
% %                     /\   |  __ \|  ____|    
% %                    /  \  | |  | | |__   ___ 
% %                   / /\ \ | |  | |  __| / __|
% %                  / ____ \| |__| | |____\__ \
% %                 /_/    \_\_____/|______|___/
% % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Let's end with a figure which considers ADE1, ADE2, and top 0.1% of the
% results of the Latin hypercube sampling

[ADE_1] = ADE_collection1(ADE1,ADE2,ADE3,Interp_curve_01,time,BTC_input,Description1,Description2);

saveas(ADE_1,[pwd '/Output_files/ADE_1.fig']);
saveas(ADE_1,[pwd '/Output_files/ADE_1.tif']);
close (ADE_1)

[ADE_2,ADE3_limits] = ADE_collection2(Hyperspace, time, Length, BTC_input,M_g,L);
saveas(ADE_2,[pwd '/Output_files/ADE_2.fig']);
saveas(ADE_2,[pwd '/Output_files/ADE_2.tif']);
close (ADE_2)

end

