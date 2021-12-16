function [Working_matrix,Order_Working_matrix,Order_Working_matrix_20,Order_Working_matrix_10, TEMP_BTC,Interp_curve] = figures_OTIS(OTIS_hypercube_input,...
    Description1, Description2,Data,Instate,OSFLAG,Not_used_param,BTC_input,str)

clear Working_matrix
Working_matrix=OTIS_hypercube_input;
    % Working_matrix(:,1) - v
    % Working_matrix(:,2) - Area
    % Working_matrix(:,3) - D
    % Working_matrix(:,4) - Alpha
    % Working_matrix(:,5) - Area transient storage
    % Working_matrix(:,6) - Discharge
    % Working_matrix(:,7) - Initial Concentration

    Working_matrix(:,8)=Data.RMSE(:,1);         % RMSE
    Working_matrix(:,9)=Data.r2(:,1);           % r2
    Working_matrix(:,10)=Data.nRMSE(:,1);       % nRMSE
    Working_matrix(:,11)=Data.logRMSE(:,1);     % logRMSE
    Working_matrix(:,12)=Data.logr2(:,1);       % logr2
    Working_matrix(:,13)=Data.Pearson_r2(:,1);  % Pearson r2
    Working_matrix(:,14)=Data.Pearson_logr2(:,1); % Pearson log r2
    Working_matrix(:,15)=Data.KGE(:,1);         % Kling gupta eff

Order_Working_matrix=sortrows(Working_matrix,8);   % Hyperspace sorted depending on RMSE value

% Create the BTC curves we need 
% Curves evaluated from the top 20%

top20=length(Working_matrix(:,1))*0.2;        % Limit for top 20% of the results
Order_Working_matrix_20=Order_Working_matrix(1:top20,:);

% % Build the top 100 BTC -> 

wbhandle=waitbar(0,'Building the top 100 OTIS BTC. Please wait.');
for i = 1:100             % define every run
%
%Run the TSM for i-th parameter set
Model=OTIS_run(Instate,i,Order_Working_matrix,OSFLAG,Not_used_param);
Sim = interp1(Model.ttime,Model.conc_Channel,BTC_input(:,1));
TEMP_BTC(:,i)=Sim;
clear Sim
% progress with the waitbar
waitbar(i / 100)
end

close(wbhandle)
% build the Area of the top results 
% top10=length(Working_matrix(:,1))*0.1;        % Limit for top 10% of the results
% top1=length(Working_matrix(:,1))*0.01;        % Limit for top 1% of the results
% top01=length(Working_matrix(:,1))*0.001;      % Limit for top 0.1% of the results
% 
% for i=1:1:length(BTC_input(:,1)) % OLD OPTION FOR TOP 20%-10%-1% OF THE
% % BTC. Discarded because massively time demanding to build all the BTC of the top 20% of the results
%     
%     Interp_curve_20(i,1)=BTC_input(i,1);
%     Interp_curve_20(i,2)=min(TEMP_BTC(i,1:top20));
%     Interp_curve_20(i,3)=max(TEMP_BTC(i,1:top20));
%     
%     Interp_curve_10(i,1)=BTC_input(i,1);
%     Interp_curve_10(i,2)=min(TEMP_BTC(i,1:top10));
%     Interp_curve_10(i,3)=max(TEMP_BTC(i,1:top10));
%     
%     Interp_curve_1(i,1)=BTC_input(i,1);
%     Interp_curve_1(i,2)=min(TEMP_BTC(i,1:top1));
%     Interp_curve_1(i,3)=max(TEMP_BTC(i,1:top1));
%     
%     Interp_curve_01(i,1)=BTC_input(i,1);
%     Interp_curve_01(i,2)=min(TEMP_BTC(i,1:top01));
%     Interp_curve_01(i,3)=max(TEMP_BTC(i,1:top01));
% end

 for i=1:1:length(BTC_input(:,1)) % I leve the option for the top 100 BTC
    Interp_curve(i,1)=BTC_input(i,1);
    Interp_curve(i,2)=min(TEMP_BTC(i,1:100));
    Interp_curve(i,3)=max(TEMP_BTC(i,1:100));     
 end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %              ______ _      __                     _ _ 
% %             |  ____(_)    /_ |              /\   | | |
% %             | |__   _  __ _| |  ______     /  \  | | |
% %             |  __| | |/ _` | | |______|   / /\ \ | | |
% %             | |    | | (_| | |           / ____ \| | |
% %             |_|    |_|\__, |_|          /_/    \_\_|_|
% %                        __/ |                          
% %                       |___/                           
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Plot all the variables (v, A, D, Alpha and A_TS) vs chosen objective
% functions (First row -> RMSE; Second row -> Nash Sutcliffe; Third row ->
% Kling-Gupta Eff.

f=figure;
% RMSE
subplot(3,5,1)
plot(Order_Working_matrix(:,1),Order_Working_matrix(:,8),'.k')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,1)) max(Order_Working_matrix(:,1))])
ylim([0 max(Order_Working_matrix(:,8))])

subplot(3,5,2)
plot(Order_Working_matrix(:,2),Order_Working_matrix(:,8),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,2)) max(Order_Working_matrix(:,2))])
ylim([0 max(Order_Working_matrix(:,8))])

subplot(3,5,3)
plot(Order_Working_matrix(:,3),Order_Working_matrix(:,8),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,3)) max(Order_Working_matrix(:,3))])
ylim([0 max(Order_Working_matrix(:,8))])

subplot(3,5,4)
plot(Order_Working_matrix(:,4),Order_Working_matrix(:,8),'.k')
xlabel ('Alpha [1/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,4)) max(Order_Working_matrix(:,4))])
ylim([0 max(Order_Working_matrix(:,8))])

subplot(3,5,5)
plot(Order_Working_matrix(:,5),Order_Working_matrix(:,8),'.k')
xlabel ('A_T_S [m^2]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,5)) max(Order_Working_matrix(:,5))])
ylim([0 max(Order_Working_matrix(:,8))])

%r^2
subplot(3,5,6)
plot(Order_Working_matrix(:,1),Order_Working_matrix(:,9),'.k')
xlabel ('v [m/s]');
ylabel ('r^2');
xlim([min(Order_Working_matrix(:,1)) max(Order_Working_matrix(:,1))])
ylim([0 1])

subplot(3,5,7)
plot(Order_Working_matrix(:,2),Order_Working_matrix(:,9),'.k')
xlabel ('A [m^2]');
ylabel ('r^2');
xlim([min(Order_Working_matrix(:,2)) max(Order_Working_matrix(:,2))])
ylim([0 1])

subplot(3,5,8)
plot(Order_Working_matrix(:,3),Order_Working_matrix(:,9),'.k')
xlabel ('D [m^2/s]');
ylabel ('r^2');
xlim([min(Order_Working_matrix(:,3)) max(Order_Working_matrix(:,3))])
ylim([0 1])

subplot(3,5,9)
plot(Order_Working_matrix(:,4),Order_Working_matrix(:,9),'.k')
xlabel ('Alpha [1/s]');
ylabel ('r^2');
xlim([min(Order_Working_matrix(:,4)) max(Order_Working_matrix(:,4))])
ylim([0 1])

subplot(3,5,10)
plot(Order_Working_matrix(:,5),Order_Working_matrix(:,9),'.k')
xlabel ('A_T_S [m^2]');
ylabel ('r^2');
xlim([min(Order_Working_matrix(:,5)) max(Order_Working_matrix(:,5))])
ylim([0 1])

% KGE
subplot(3,5,11)
plot(Order_Working_matrix(:,1),Order_Working_matrix(:,15),'.k')
xlabel ('v [m/s]');
ylabel ('KGE');
xlim([min(Order_Working_matrix(:,1)) max(Order_Working_matrix(:,1))])
ylim([0 1])

subplot(3,5,12)
plot(Order_Working_matrix(:,2),Order_Working_matrix(:,15),'.k')
xlabel ('A [m^2]');
ylabel ('KGE');
xlim([min(Order_Working_matrix(:,2)) max(Order_Working_matrix(:,2))])
ylim([0 1])

subplot(3,5,13)
plot(Order_Working_matrix(:,3),Order_Working_matrix(:,15),'.k')
xlabel ('D [m^2/s]');
ylabel ('KGE');
xlim([min(Order_Working_matrix(:,3)) max(Order_Working_matrix(:,3))])
ylim([0 1])

subplot(3,5,14)
plot(Order_Working_matrix(:,4),Order_Working_matrix(:,15),'.k')
xlabel ('Alpha [1/s]');
ylabel ('KGE');
xlim([min(Order_Working_matrix(:,4)) max(Order_Working_matrix(:,4))])
ylim([0 1])

subplot(3,5,15)
plot(Order_Working_matrix(:,5),Order_Working_matrix(:,15),'.k')
xlabel ('A_T_S [m^2]');
ylabel ('KGE');
xlim([min(Order_Working_matrix(:,5)) max(Order_Working_matrix(:,5))])
ylim([0 1])

annotation('textbox', [0.01, 0.8, 0.1, 0.1], 'String',str,'FontSize',11,'LineStyle','none')
set(gcf, 'WindowState', 'maximized');
sgtitle({'OTIS latin hypercube';Description1;Description2},'FontSize',14);

saveas(f,[pwd '/Output_files_OTIS/OTIS_perf_all.fig']);
saveas(f,[pwd '/Output_files_OTIS/OTIS_perf_all.tif']);
close(f)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %          ______         ___                      _ _ ___  
% %         |  ____(_)     |__ \               /\   | | |__ \ 
% %         | |__   _  __ _   ) |  ______     /  \  | | |  ) |
% %         |  __| | |/ _` | / /  |______|   / /\ \ | | | / / 
% %         | |    | | (_| |/ /_            / ____ \| | |/ /_ 
% %         |_|    |_|\__, |____|          /_/    \_\_|_|____|
% %                    __/ |                                  
% %                   |___/                                   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

f=figure;           % identical to figure 1, we just plot the top 20%
% RMSE
subplot(3,5,1)
plot(Order_Working_matrix_20(:,1),Order_Working_matrix_20(:,8),'.k')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_20(:,1)) max(Order_Working_matrix_20(:,1))])
ylim([0 max(Order_Working_matrix_20(:,8))])

subplot(3,5,2)
plot(Order_Working_matrix_20(:,2),Order_Working_matrix_20(:,8),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_20(:,2)) max(Order_Working_matrix_20(:,2))])
ylim([0 max(Order_Working_matrix_20(:,8))])

subplot(3,5,3)
plot(Order_Working_matrix_20(:,3),Order_Working_matrix_20(:,8),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_20(:,3)) max(Order_Working_matrix_20(:,3))])
ylim([0 max(Order_Working_matrix_20(:,8))])

subplot(3,5,4)
plot(Order_Working_matrix_20(:,4),Order_Working_matrix_20(:,8),'.k')
xlabel ('Alpha [1/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_20(:,4)) max(Order_Working_matrix_20(:,4))])
ylim([0 max(Order_Working_matrix_20(:,8))])

subplot(3,5,5)
plot(Order_Working_matrix_20(:,5),Order_Working_matrix_20(:,8),'.k')
xlabel ('A_T_S [m^2]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_20(:,5)) max(Order_Working_matrix_20(:,5))])
ylim([0 max(Order_Working_matrix_20(:,8))])

%r^2
subplot(3,5,6)
plot(Order_Working_matrix_20(:,1),Order_Working_matrix_20(:,9),'.k')
xlabel ('v [m/s]');
ylabel ('r^2');
xlim([min(Order_Working_matrix_20(:,1)) max(Order_Working_matrix_20(:,1))])
ylim([min(Order_Working_matrix_20(:,9)) 1])

subplot(3,5,7)
plot(Order_Working_matrix_20(:,2),Order_Working_matrix_20(:,9),'.k')
xlabel ('A [m^2]');
ylabel ('r^2');
xlim([min(Order_Working_matrix_20(:,2)) max(Order_Working_matrix_20(:,2))])
ylim([min(Order_Working_matrix_20(:,9)) 1])

subplot(3,5,8)
plot(Order_Working_matrix_20(:,3),Order_Working_matrix_20(:,9),'.k')
xlabel ('D [m^2/s]');
ylabel ('r^2');
xlim([min(Order_Working_matrix_20(:,3)) max(Order_Working_matrix_20(:,3))])
ylim([min(Order_Working_matrix_20(:,9)) 1])

subplot(3,5,9)
plot(Order_Working_matrix_20(:,4),Order_Working_matrix_20(:,9),'.k')
xlabel ('Alpha [1/s]');
ylabel ('r^2');
xlim([min(Order_Working_matrix_20(:,4)) max(Order_Working_matrix_20(:,4))])
ylim([min(Order_Working_matrix_20(:,9)) 1])

subplot(3,5,10)
plot(Order_Working_matrix_20(:,5),Order_Working_matrix_20(:,9),'.k')
xlabel ('A_T_S [m^2]');
ylabel ('r^2');
xlim([min(Order_Working_matrix_20(:,5)) max(Order_Working_matrix_20(:,5))])
ylim([min(Order_Working_matrix_20(:,9)) 1])

% KGE
subplot(3,5,11)
plot(Order_Working_matrix_20(:,1),Order_Working_matrix_20(:,15),'.k')
xlabel ('v [m/s]');
ylabel ('KGE');
xlim([min(Order_Working_matrix_20(:,1)) max(Order_Working_matrix_20(:,1))])
ylim([min(Order_Working_matrix_20(:,15)) 1])

subplot(3,5,12)
plot(Order_Working_matrix_20(:,2),Order_Working_matrix_20(:,15),'.k')
xlabel ('A [m^2]');
ylabel ('KGE');
xlim([min(Order_Working_matrix_20(:,2)) max(Order_Working_matrix_20(:,2))])
ylim([min(Order_Working_matrix_20(:,15)) 1])

subplot(3,5,13)
plot(Order_Working_matrix_20(:,3),Order_Working_matrix_20(:,15),'.k')
xlabel ('D [m^2/s]');
ylabel ('KGE');
xlim([min(Order_Working_matrix_20(:,3)) max(Order_Working_matrix_20(:,3))])
ylim([min(Order_Working_matrix_20(:,15)) 1])

subplot(3,5,14)
plot(Order_Working_matrix_20(:,4),Order_Working_matrix_20(:,15),'.k')
xlabel ('Alpha [1/s]');
ylabel ('KGE');
xlim([min(Order_Working_matrix_20(:,4)) max(Order_Working_matrix_20(:,4))])
ylim([min(Order_Working_matrix_20(:,15)) 1])

subplot(3,5,15)
plot(Order_Working_matrix_20(:,5),Order_Working_matrix_20(:,15),'.k')
xlabel ('A_T_S [m^2]');
ylabel ('KGE');
xlim([min(Order_Working_matrix_20(:,5)) max(Order_Working_matrix_20(:,5))])
ylim([min(Order_Working_matrix_20(:,15)) 1])

set(gcf, 'WindowState', 'maximized');
sgtitle({'OTIS latin hypercube - top20%';Description1;Description2},'FontSize',14);
saveas(f,[pwd '/Output_files_OTIS/OTIS_perf_20.fig']);
saveas(f,[pwd '/Output_files_OTIS/OTIS_perf_20.tif']);

close (f)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %  ______ _       ____             _                      _       _       
% % |  ____(_)     |___ \           | |                    | |     | |      
% % | |__   _  __ _  __) |  ______  | |_ ___  _ __    _ __ | | ___ | |_ ___ 
% % |  __| | |/ _` ||__ <  |______| | __/ _ \| '_ \  | '_ \| |/ _ \| __/ __|
% % | |    | | (_| |___) |          | || (_) | |_) | | |_) | | (_) | |_\__ \
% % |_|    |_|\__, |____/            \__\___/| .__/  | .__/|_|\___/ \__|___/
% %            __/ |                         | |     | |                    
% %           |___/                          |_|     |_|                    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

% % % In this option we just build the top 100 BTC 
 formatSpec1="v_m_i_n=%g m/s";
 formatSpec2="v_m_a_x=%g m/s";
 formatSpec3="A_m_i_n=%g m^2";
 formatSpec4="A_m_a_x=%g m^2";
 formatSpec5="D_m_i_n=%g m^2/s"; 
 formatSpec6="D_m_a_x=%g m^2/s";
 formatSpec7="Alpha_m_i_n=%g 1/s";
 formatSpec8="Alpha_m_a_x=%g 1/s";
 formatSpec9="A_T_S_m_i_n=%g m^2";
 formatSpec10="A_T_S_m_a_x=%g m^2";
 formatSpec11="RMSE_m_i_n=%g m^2";
 formatSpec12="RMSE_m_a_x=%g m^2";
% 
 ssstr(1,1)=sprintf(formatSpec1,min(Order_Working_matrix(1:100,1)));
 ssstr(2,1)=sprintf(formatSpec2,max(Order_Working_matrix(1:100,1)));
 ssstr(3,1)=sprintf(formatSpec3,min(Order_Working_matrix(1:100,2)));
 ssstr(4,1)=sprintf(formatSpec4,max(Order_Working_matrix(1:100,2)));
 ssstr(5,1)=sprintf(formatSpec5,min(Order_Working_matrix(1:100,3)));
 ssstr(6,1)=sprintf(formatSpec6,max(Order_Working_matrix(1:100,3)));
 ssstr(7,1)=sprintf(formatSpec7,min(Order_Working_matrix(1:100,4)));
 ssstr(8,1)=sprintf(formatSpec8,max(Order_Working_matrix(1:100,4)));
 ssstr(9,1)=sprintf(formatSpec9,min(Order_Working_matrix(1:100,5)));
 ssstr(10,1)=sprintf(formatSpec10,max(Order_Working_matrix(1:100,5)));
 ssstr(11,1)=sprintf(formatSpec11,min(Order_Working_matrix(1:100,8)));
 ssstr(12,1)=sprintf(formatSpec12,max(Order_Working_matrix(1:100,8)));

clear formatSpec1 formatSpec2 formatSpec3 formatSpec4 formatSpec4 formatSpec6...
    formatSpec7 formatSpec8 formatSpec9 formatSpec10 formatSpec11 formatSpec12
 
f=figure;

subplot(2,3,1)
plot(Order_Working_matrix(1:100,1),Order_Working_matrix(1:100,8),'.k')    % v VS RMSE
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(1:100,1)) max(Order_Working_matrix(1:100,1))])
ylim([0 max(Order_Working_matrix(1:100,8))])
% 
subplot(2,3,2)
plot(Order_Working_matrix(1:100,2),Order_Working_matrix(1:100,8),'.k')    % A VS RMSE
xlabel ('A [m^2]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(1:100,2)) max(Order_Working_matrix(1:100,2))])
ylim([0 max(Order_Working_matrix(1:100,8))])
% 
subplot(2,3,3)
plot(Order_Working_matrix(1:100,3),Order_Working_matrix(1:100,8),'.k')    % D VS RMSE
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(1:100,3)) max(Order_Working_matrix(1:100,3))])
ylim([0 max(Order_Working_matrix(1:100,8))])
% 
subplot(2,3,4)
plot(Order_Working_matrix(1:100,4),Order_Working_matrix(1:100,8),'.k')    % Alpha VS RMSE
xlabel ('Alpha [1/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(1:100,4)) max(Order_Working_matrix(1:100,4))])
ylim([0 max(Order_Working_matrix(1:100,8))])
% 
subplot(2,3,5)
plot(Order_Working_matrix(1:100,5),Order_Working_matrix(1:100,8),'.k')    % A_TS VS RMSE
xlabel ('A_T_S [m^2]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(1:100,5)) max(Order_Working_matrix(1:100,5))])
ylim([0 max(Order_Working_matrix(1:100,8))])
%
subplot(2,3,6)
plot(Interp_curve(:,1), Interp_curve(:,2), 'k', 'LineWidth', 1,'HandleVisibility','off');
hold on;
plot(Interp_curve(:,1), Interp_curve(:,3), 'k', 'LineWidth', 1,'HandleVisibility','off');
x2 = [Interp_curve(:,1)', fliplr(Interp_curve(:,1)')];
inBetween = [Interp_curve(:,2)', fliplr(Interp_curve(:,3)')];
fill(x2, inBetween, [0.86,0.86,0.86]);
hold on
plot(BTC_input(:,1),TEMP_BTC(:,1),'--b','LineWidth',1)
hold on
plot(BTC_input(:,1),BTC_input(:,2),'-r','LineWidth',2)
xlabel ('t [s]');
ylabel ('Conc [mg/l]');
legend('Simulation range','Best-fitting OTIS BTC','Observed BTC')
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String', "Latin Hypercube sampling - top 100 results",...
    'FontSize',12,'LineStyle','none','FitBoxToText','on')
annotation('textbox', [0.01, 0.8, 0.1, 0.1], 'String',str,'FontSize',11,'LineStyle','none')
set(gcf, 'WindowState', 'maximized');

sgtitle({'OTIS latin hypercube';Description1;Description2},'FontSize',14);
% 
saveas(f,[pwd '/Output_files_OTIS/OTIS_top100.fig']);
saveas(f,[pwd '/Output_files_OTIS/OTIS_top100.tif']);
% 
close (f)

% Let's compute the statistics for only the best 10% of the simulations
% based on RMSE values -> Following Ward et al., 2017

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %                    _____ _        _     __ 
% %                   / ____| |      | |   /_ |
% %                  | (___ | |_ __ _| |_   | |
% %                   \___ \| __/ _` | __|  | |
% %                   ____) | || (_| | |_   | |
% %                  |_____/ \__\__,_|\__|  |_|
% %                           
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Compute the statistics only on the top 10% of the results
% TAKE ONLY THE TOP 10% OF THE DISCENDING (in RMSE) WORKING MATRIX

Order_Working_matrix_10=sortrows(Working_matrix,8);
% Take the 10%
rem=length(Working_matrix)*0.1;
Order_Working_matrix_10=Order_Working_matrix_10(1:rem,:);
clear rem

%%%
[a,b] = sort(Order_Working_matrix_10(:,8)); % Sort depending on RMSE

% vector "a" puts the RMSE in crescent order while vector "b" returns a 
% the location of every RMSE value in the "a" vector

n = length(a);
nbins = 25;      % number of intervals to create the frequency plots

% thresholds corresponding to the top 20%, 10%, 5%, 1%, 0.5% and 0.01% 
% of the RMSE values --> if you want, you can change them
m = [0.2 0.1 0.05 0.01 0.005 0.001];   

% y-axis ranges
minhist = 0;
maxhist = 80;     
MAP = jet(256);

%%%%%%%%%%%%%%%%%
% Figure 1 = v, A, D, As and Alpha data vs RMSE; 
% Frequency plots and normalized cumulative distribution

stat1=figure;

subplot(3,5,1)
plot(Order_Working_matrix_10(:,1),Order_Working_matrix_10(:,8),'.k')
hold on
plot(Order_Working_matrix_10(1,1),Order_Working_matrix_10(1,8),'ro')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_10(:,1)) max(Order_Working_matrix_10(:,1))])
ylim([0 max(Order_Working_matrix_10(:,8))])

subplot(3,5,2)
plot(Order_Working_matrix_10(:,2),Order_Working_matrix_10(:,8),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
hold on
plot(Order_Working_matrix_10(1,2),Order_Working_matrix_10(1,8),'ro')
xlim([min(Order_Working_matrix_10(:,2)) max(Order_Working_matrix_10(:,2))])
ylim([0 max(Order_Working_matrix_10(:,8))])

subplot(3,5,3)
plot(Order_Working_matrix_10(:,3),Order_Working_matrix_10(:,8),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_10(:,3)) max(Order_Working_matrix_10(:,3))])
ylim([0 max(Order_Working_matrix_10(:,8))])
hold on
plot(Order_Working_matrix_10(1,3),Order_Working_matrix_10(1,8),'ro')

subplot(3,5,4)
plot(Order_Working_matrix_10(:,4),Order_Working_matrix_10(:,8),'.k')
xlabel ('Alpha [1/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_10(:,4)) max(Order_Working_matrix_10(:,4))])
ylim([0 max(Order_Working_matrix_10(:,8))])
hold on
plot(Order_Working_matrix_10(1,4),Order_Working_matrix_10(1,8),'ro')

subplot(3,5,5)
plot(Order_Working_matrix_10(:,5),Order_Working_matrix_10(:,8),'.k')
xlabel ('A_T_S [m^2]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_10(:,5)) max(Order_Working_matrix_10(:,5))])
ylim([0 max(Order_Working_matrix_10(:,8))])
hold on
plot(Order_Working_matrix_10(1,5),Order_Working_matrix_10(1,8),'ro')

% Second row - frequency curves 
% this is parameter value probability distribution functions that display
% the distribution of the parameter values CORRESPONDING to a certin
% threshold range. If there is narrow PDF around a certain range/value -> the
% parameter is certain and the peak is around the best value
clear parameter temp p1 c1 d1 

subplot (3,5,6)
parameter=Order_Working_matrix_10(:,1);     % Select velocity 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end
plot([parameter(b(1)) parameter(b(1))],[0 100*max(c1)./(sum(c1))+10],'Color',[0.5 0.5 0.5],'LineWidth',1);
if sum(c1)==0
        axis([min(parameter) max(parameter) 0 100]);
    else
        axis([min(parameter) max(parameter) 0 100*max(c1)./(sum(c1))+10]);
end 
legend(num2str(m'))
xlabel ('v [m/s]');           
ylabel('Frequency')   
clear parameter temp p1 c1 d1 

subplot (3,5,7)
parameter=Order_Working_matrix_10(:,2);     % Select Area 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end
plot([parameter(b(1)) parameter(b(1))],[0 100*max(c1)./(sum(c1))+10],'Color',[0.5 0.5 0.5],'LineWidth',1);
if sum(c1)==0
        axis([min(parameter) max(parameter) 0 100]);
    else
        axis([min(parameter) max(parameter) 0 100*max(c1)./(sum(c1))+10]);
end 
legend(num2str(m'))
xlabel ('A [m^2]');           
ylabel('Frequency')  
clear parameter temp p1 c1 d1 

subplot (3,5,8)
parameter=Order_Working_matrix_10(:,3);     % Select Dispersion 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end
plot([parameter(b(1)) parameter(b(1))],[0 100*max(c1)./(sum(c1))+10],'Color',[0.5 0.5 0.5],'LineWidth',1);
if sum(c1)==0
        axis([min(parameter) max(parameter) 0 100]);
    else
        axis([min(parameter) max(parameter) 0 100*max(c1)./(sum(c1))+10]);
end 
legend(num2str(m'))
xlabel ('D [m^2/s]');           
ylabel('Frequency')  
clear parameter temp p1 c1 d1 

subplot (3,5,9)
parameter=Order_Working_matrix_10(:,4);     % Select Alpha 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end
plot([parameter(b(1)) parameter(b(1))],[0 100*max(c1)./(sum(c1))+10],'Color',[0.5 0.5 0.5],'LineWidth',1);
if sum(c1)==0
        axis([min(parameter) max(parameter) 0 100]);
    else
        axis([min(parameter) max(parameter) 0 100*max(c1)./(sum(c1))+10]);
end 
legend(num2str(m'))
xlabel ('Alpha [1/s]');           
ylabel('Frequency')  
clear parameter temp p1 c1 d1 

subplot (3,5,10)
parameter=Order_Working_matrix_10(:,5);     % Select Area transient storage 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end
plot([parameter(b(1)) parameter(b(1))],[0 100*max(c1)./(sum(c1))+10],'Color',[0.5 0.5 0.5],'LineWidth',1);
if sum(c1)==0
        axis([min(parameter) max(parameter) 0 100]);
    else
        axis([min(parameter) max(parameter) 0 100*max(c1)./(sum(c1))+10]);
end 
legend(num2str(m'))
xlabel ('A_T_S [m^2]');           
ylabel('Frequency')  
clear parameter temp p1 c1 d1 

% % % % % % % % 
% Third row - normalized cumulative (Cum. norm.) distributions of RMSE 
% lines represent the best 10% of model runs binned into 1% increments; 
% each line represents 1% of all model simulations
% % % % % % % % 
clear h of dat tmx tmy j

of=Order_Working_matrix_10(:,8);       % Select objective function RMSE
                               % of is automatically sorted from the lower
                               % RMSE to the higher RMSE
of=of./max(of);                % normalise of
of=1-of;                       % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods
[y,i]=sort(of);
dat=Order_Working_matrix_10(i,:);
cls=floor(length(dat)/10);
tmx=zeros(cls,10);
tmy=tmx;

subplot (3,5,11)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
    tm=dat(cls*(j-1)+1:cls*j,1);     % Take first column -> velocity
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls;
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['v [m/s]'])
  ylabel(['cum. norm. RMSE'])
  
colormap(jet(10));
h=axes('position',[.07 .07 .02 .33]);
h=colorbar(h);
set(h,'ytick',[2 10]);
set(h,'yticklabel',['L';'H']);
ylabel(['Likelihood RMSE'])

clear h 

subplot (3,5,12)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
    tm=dat(cls*(j-1)+1:cls*j,2);    % Take second column now -> Area
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls;
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['A [m^2]'])
  ylabel(['cum. norm. RMSE'])
  
subplot (3,5,13)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
    tm=dat(cls*(j-1)+1:cls*j,3);    % Take third column now -> Disp
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls;
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['D [m^2/s]'])
  ylabel(['cum. norm. RMSE'])
  
subplot (3,5,14)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
    tm=dat(cls*(j-1)+1:cls*j,4);    % Take fourth column now -> Alpha
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls;
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['Alpha [1/s]'])
  ylabel(['cum. norm. RMSE'])
  
  subplot (3,5,15)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
    tm=dat(cls*(j-1)+1:cls*j,5);    % Take fifth column now -> Area transient storage
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls;
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['A_T_S [m^2]'])
  ylabel(['cum. norm. RMSE'])
  
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'top 10% latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({'OTIS latin hypercube';Description1;Description2},'FontSize',14);
  
clear h clear of dat tmx tmy j of cls 
set(gcf, 'WindowState', 'maximized');

saveas(stat1,[pwd '/Output_files_OTIS/OTIS_stat1.fig']);
saveas(stat1,[pwd '/Output_files_OTIS/OTIS_stat1.tif']);
close (stat1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %                    _____ _        _     ___  
% %                   / ____| |      | |   |__ \ 
% %                  | (___ | |_ __ _| |_     ) |
% %                   \___ \| __/ _` | __|   / / 
% %                   ____) | || (_| | |_   / /_ 
% %                  |_____/ \__\__,_|\__| |____|
% %                             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Figure 2 = v, A and D data vs RMSE; A posteriori parameter distribution;
% parameter identifiability

stat2=figure;

subplot(3,5,1)
plot(Order_Working_matrix_10(:,1),Order_Working_matrix_10(:,8),'.k')
hold on
% [P,PP] = min(Order_Working_matrix_10(:,8)); % Ir's already sorted so it's
% of course the first line
plot(Order_Working_matrix_10(1,1),Order_Working_matrix_10(1,8),'ro')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_10(:,1)) max(Order_Working_matrix_10(:,1))])
ylim([0 max(Order_Working_matrix_10(:,8))])

subplot(3,5,2)
plot(Order_Working_matrix_10(:,2),Order_Working_matrix_10(:,8),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
hold on
plot(Order_Working_matrix_10(1,2),Order_Working_matrix_10(1,8),'ro')
xlim([min(Order_Working_matrix_10(:,2)) max(Order_Working_matrix_10(:,2))])
ylim([0 max(Order_Working_matrix_10(:,8))])

subplot(3,5,3)
plot(Order_Working_matrix_10(:,3),Order_Working_matrix_10(:,8),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_10(:,3)) max(Order_Working_matrix_10(:,3))])
ylim([0 max(Order_Working_matrix_10(:,8))])
hold on
plot(Order_Working_matrix_10(1,3),Order_Working_matrix_10(1,8),'ro')

subplot(3,5,4)
plot(Order_Working_matrix_10(:,4),Order_Working_matrix_10(:,8),'.k')
xlabel ('Alpha [1/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_10(:,4)) max(Order_Working_matrix_10(:,4))])
ylim([0 max(Order_Working_matrix_10(:,8))])
hold on
plot(Order_Working_matrix_10(1,4),Order_Working_matrix_10(1,8),'ro')

subplot(3,5,5)
plot(Order_Working_matrix_10(:,5),Order_Working_matrix_10(:,8),'.k')
xlabel ('A_T_S [m^2]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix_10(:,5)) max(Order_Working_matrix_10(:,5))])
ylim([0 max(Order_Working_matrix_10(:,8))])
hold on
plot(Order_Working_matrix_10(1,5),Order_Working_matrix_10(1,8),'ro')

% % % % % % % % 
% ROW 2 -> A POSTERIORI PARAMETER DISTRIBUTION
% % % % % % % % 
% Note that a posteriori parameter distribution plots are not always
% descriptive. Especially in not-identifiable problems is not rare to have
% the top 10% of the results displaying a cloud or a plateau of values 
% where the best results approach a certain RMSE (or general likelyhood) 
% value, and having also few tens of values, isolated from the cloud,
% having the same likelyhood. This is a classical result for undersampling.
% The results of this plot will give the bar including that single/few
% value(s) a very high height, which is the mean likelyhood of the
% container... but if the container has just one value with a relatively
% high likelyhood then we will have an artificial large height for this
% parameter range, which do not really mean anything for our distribution. 
% In this case the result has to be interpreted as the model is
% not-identifiable and we need larger amount of simulation and, probably,
% some restriction in the parameter range before the montecarlo/latin
% hypercube sampling

ncontainers=20;
slider_value=100;   % We take the 100% of the values because we already have seleceted Hyperspace_best as top 10%
% In case we have no initial cut in the Hyperspace matrix, we should have
% set slider_value=10, to select top 10% of the results

% calculate likelihood
of=Order_Working_matrix_10(:,8);  % select the criterion -> We've already sorted from the lower to the higher RMSE
                                  % RMSE -> low values indicate better models
                                
LL=1-of; 						% likelihood (high values indicate more likely [probable] models)
if min(LL)<0|min(LL)==0, LL=LL-min(LL)+1000*eps;end; % transform negative lhoods

LL=LL./sum(LL);     % sum(likelihoods)=1      problem if NaN in vector

%Eliminate data that is below the slider threshold
LL=sortrows(LL);
dat=sortrows(Order_Working_matrix_10,8);     % sort Hyperspace based on the RMSE values
LL = flipud(LL);
numdat = floor((slider_value / 100) * size(LL));
LL(numdat+1:size(LL),:) = [];
dat(numdat+1:size(dat),:) = [];

subplot(3,5,6) % velocity

width=(max(dat(:,1))-min(dat(:,1)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,1)>min(dat(:,1))+(n-1)*width&dat(:,1)<min(dat(:,1))+n*width); % find all values within the container
      temp1=dat(K,1);
      temp2=LL(K);
      
      bx(1,n)=min(dat(:,1))+.5*width+(n-1)*width;   % middle of container
      by(1,n)=sum(temp2)/length(temp2);             % mean likelihood in container
      
      clear temp 1 temp2 K
end 
bar(bx(1,:),by(1,:),'k');hold on;
ymax=(max(by(1,:))+0.1*(max(by(1,:))));
ymin=(min(by(1,:))-0.1*(min(by(1,:))));
axis([min(dat(:,1)) max(dat(:,1)) ymin ymax]);
ylabel(['D (RMSE)'])
xlabel('v [m/s]');
   
subplot(3,5,7) % Area

width=(max(dat(:,2))-min(dat(:,2)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,2)>min(dat(:,2))+(n-1)*width&dat(:,2)<min(dat(:,2))+n*width); % find all values within the container
      temp1=dat(K,2);
      temp2=LL(K);
      
      bx(2,n)=min(dat(:,2))+.5*width+(n-1)*width;   % middle of container
      by(2,n)=sum(temp2)/length(temp2);             % mean likelihood in container
      
      clear temp 1 temp2 K
end 
bar(bx(2,:),by(2,:),'k');hold on;
ymax=(max(by(2,:))+0.1*(max(by(2,:))));
ymin=(min(by(2,:))-0.1*(min(by(2,:))));
axis([min(dat(:,2)) max(dat(:,2)) ymin ymax]);
ylabel(['D (RMSE)'])
xlabel('A [m^2]');

subplot(3,5,8) % Disp

width=(max(dat(:,3))-min(dat(:,3)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,3)>min(dat(:,3))+(n-1)*width&dat(:,3)<min(dat(:,3))+n*width); % find all values within the container
      temp1=dat(K,3);
      temp2=LL(K);
      
      bx(3,n)=min(dat(:,3))+.5*width+(n-1)*width;   % middle of container
      by(3,n)=sum(temp2)/length(temp2);             % mean likelihood in container
      
      clear temp 1 temp2 K
end 
bar(bx(3,:),by(3,:),'k');hold on;
ymax=(max(by(3,:))+0.1*(max(by(3,:))));
ymin=(min(by(3,:))-0.1*(min(by(3,:))));
axis([min(dat(:,3)) max(dat(:,3)) ymin ymax]);
ylabel(['D (RMSE)'])
xlabel('D [m^2/s]');

subplot(3,5,9) % Alpha

width=(max(dat(:,4))-min(dat(:,4)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,4)>min(dat(:,4))+(n-1)*width&dat(:,4)<min(dat(:,4))+n*width); % find all values within the container
      temp1=dat(K,4);
      temp2=LL(K);
      
      bx(4,n)=min(dat(:,4))+.5*width+(n-1)*width;   % middle of container
      by(4,n)=sum(temp2)/length(temp2);             % mean likelihood in container
      
      clear temp 1 temp2 K
end 
bar(bx(4,:),by(4,:),'k');hold on;
ymax=(max(by(4,:))+0.1*(max(by(4,:))));
ymin=(min(by(4,:))-0.1*(min(by(4,:))));
axis([min(dat(:,4)) max(dat(:,4)) ymin ymax]);
ylabel(['D (RMSE)'])
xlabel('Alpha [1/s]');

subplot(3,5,10) % A_trans storage

width=(max(dat(:,5))-min(dat(:,5)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,5)>min(dat(:,5))+(n-1)*width&dat(:,5)<min(dat(:,5))+n*width); % find all values within the container
      temp1=dat(K,5);
      temp2=LL(K);
      
      bx(5,n)=min(dat(:,5))+.5*width+(n-1)*width;   % middle of container
      by(5,n)=sum(temp2)/length(temp2);             % mean likelihood in container
      
      clear temp 1 temp2 K
end 
bar(bx(5,:),by(5,:),'k');hold on;
ymax=(max(by(5,:))+0.1*(max(by(5,:))));
ymin=(min(by(5,:))-0.1*(min(by(5,:))));
axis([min(dat(:,5)) max(dat(:,5)) ymin ymax]);
ylabel(['D (RMSE)'])
xlabel('A_T_S [m^2]');

% % % % % % % % 
% 3rd ROW - IDENTIFIABILITY PLOTS
% % % % % % % % 
clear of I J dat numdat cls tmx i
dat=Order_Working_matrix_10;   % dat values sorted from the lower RMSE to the higher RMSE
containers=10; % division of each parameter range
grouping=10;  % number of groups

% calculate likelihood
of=Order_Working_matrix_10(:,8);  % criteria (low values indicate better models)
of=1-of; % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods

% sort data according to selected perf
[I,J]=sort(of); 
dat=dat(J,:);
% Eliminate data that is below the slider threshold
numdat = floor((slider_value / 100) * size(dat));
dat(numdat+1:size(dat),:) = [];

cls=floor(length(dat)/grouping);
tmx=zeros(cls,grouping);tmy=tmx;

subplot(3,5,11) % velocity
i=1;
for j=1:grouping
    tm=dat(cls*(j-1)+1:cls*j,i);
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; 
end
% calculate and plot gradients ***************************
  
   step=(max(dat(:,i))-min(dat(:,i)))/containers; % it is necessary that the boundaries are the same
   XI=[min(dat(:,i)):step:max(dat(:,i))]; % if a parameter id should be compared between two models
   
   % test for monotony
   for ii=1:5
      [II]=find(diff(tmx(1:length(tmx(:,grouping))-1,grouping))==0);
      tmx(II+1,grouping)=tmx(II+1,grouping)+[(tmx(II+2,grouping)-tmx(II+1,grouping))/2];
   end               
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping); max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
   
   [FX] = gradient(YI); % calculate gradient within containers
   
   ID_max(i)=max(FX); % keep the maximum gradient for information
   III=find(FX==max(FX));
   pos_max(i)=XI(III);
   
   ID(i,:)=FX; % keep the whole vector for each parameter
   pos(i,:)=XI;
   
   hpatches=bar('v6',XI,2*FX,'g');hold on;
   xd = get(hpatches,'xdata');
   yd = get(hpatches,'ydata');
   
   for n=1:size(yd,2)
     temp=max(yd(:,n));     % color value needs to be between 0 and 1
     tcolor=[abs(1-temp) abs(1-temp) abs(1-temp)];
     patch(xd(:,n),yd(:,n),tcolor);hold on;
   end

%    xlswrite(['wid_FX',num2str(i),'.xls'],FX);
   
   clear FX
   
   % ********************************************************
   
   plot([min(dat(:,i)); tmx(:,grouping); max(dat(:,i))],[0 tmy 1],'color','k','linewidth',2);hold on;
   grid off;
   axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
   ylabel(['Cum. Dist. RMSE'])
   xlabel('v [m/s]');
    
subplot(3,5,12) % Area
i=2;   
 for j=1:grouping
    tm=dat(cls*(j-1)+1:cls*j,i);
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; 
end
% calculate and plot gradients ***************************
  
   step=(max(dat(:,i))-min(dat(:,i)))/containers; % it is necessary that the boundaries are the same
   XI=[min(dat(:,i)):step:max(dat(:,i))]; % if a parameter id should be compared between two models
   
   % test for monotony
   for ii=1:5
      [II]=find(diff(tmx(1:length(tmx(:,grouping))-1,grouping))==0);
      tmx(II+1,grouping)=tmx(II+1,grouping)+[(tmx(II+2,grouping)-tmx(II+1,grouping))/2];
   end               
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping); max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
   
   [FX] = gradient(YI); % calculate gradient within containers
   
   ID_max(i)=max(FX); % keep the maximum gradient for information
   III=find(FX==max(FX));
   pos_max(i)=XI(III);
   
   ID(i,:)=FX; % keep the whole vector for each parameter
   pos(i,:)=XI;
   
   hpatches=bar('v6',XI,2*FX,'g');hold on;
   xd = get(hpatches,'xdata');
   yd = get(hpatches,'ydata');
   
   for n=1:size(yd,2)
     temp=max(yd(:,n));     % color value needs to be between 0 and 1
     tcolor=[abs(1-temp) abs(1-temp) abs(1-temp)];
     patch(xd(:,n),yd(:,n),tcolor);hold on;
   end

%    xlswrite(['wid_FX',num2str(i),'.xls'],FX);
   
   clear FX
   
   % ********************************************************
   
   plot([min(dat(:,i)); tmx(:,grouping); max(dat(:,i))],[0 tmy 1],'color','k','linewidth',2);hold on;
   grid off;
   axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);  
   ylabel(['Cum. Dist. RMSE'])
   xlabel('A [m^2]');
   
subplot(3,5,13) % Disp
i=3;   
 for j=1:grouping
    tm=dat(cls*(j-1)+1:cls*j,i);
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; 
end
% calculate and plot gradients ***************************
  
   step=(max(dat(:,i))-min(dat(:,i)))/containers; % it is necessary that the boundaries are the same
   XI=[min(dat(:,i)):step:max(dat(:,i))]; % if a parameter id should be compared between two models
   
   % test for monotony
   for ii=1:5
      [II]=find(diff(tmx(1:length(tmx(:,grouping))-1,grouping))==0);
      tmx(II+1,grouping)=tmx(II+1,grouping)+[(tmx(II+2,grouping)-tmx(II+1,grouping))/2];
   end               
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping); max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
   
   [FX] = gradient(YI); % calculate gradient within containers
   
   ID_max(i)=max(FX); % keep the maximum gradient for information
   III=find(FX==max(FX));
   pos_max(i)=XI(III);
   
   ID(i,:)=FX; % keep the whole vector for each parameter
   pos(i,:)=XI;
   
   hpatches=bar('v6',XI,2*FX,'g');hold on;
   xd = get(hpatches,'xdata');
   yd = get(hpatches,'ydata');
   
   for n=1:size(yd,2)
     temp=max(yd(:,n));     % color value needs to be between 0 and 1
     tcolor=[abs(1-temp) abs(1-temp) abs(1-temp)];
     patch(xd(:,n),yd(:,n),tcolor);hold on;
   end

%    xlswrite(['wid_FX',num2str(i),'.xls'],FX);
   
   clear FX
   
   % ********************************************************
   
   plot([min(dat(:,i)); tmx(:,grouping); max(dat(:,i))],[0 tmy 1],'color','k','linewidth',2);hold on;
   grid off;
   axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);     
   ylabel(['Cum. Dist. RMSE'])
   xlabel('D [m^2/s]');
   
   subplot(3,5,14) % Alpha
i=4;   
 for j=1:grouping
    tm=dat(cls*(j-1)+1:cls*j,i);
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; 
end
% calculate and plot gradients ***************************
  
   step=(max(dat(:,i))-min(dat(:,i)))/containers; % it is necessary that the boundaries are the same
   XI=[min(dat(:,i)):step:max(dat(:,i))]; % if a parameter id should be compared between two models
   
   % test for monotony
   for ii=1:5
      [II]=find(diff(tmx(1:length(tmx(:,grouping))-1,grouping))==0);
      tmx(II+1,grouping)=tmx(II+1,grouping)+[(tmx(II+2,grouping)-tmx(II+1,grouping))/2];
   end               
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping); max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
   
   [FX] = gradient(YI); % calculate gradient within containers
   
   ID_max(i)=max(FX); % keep the maximum gradient for information
   III=find(FX==max(FX));
   pos_max(i)=XI(III);
   
   ID(i,:)=FX; % keep the whole vector for each parameter
   pos(i,:)=XI;
   
   hpatches=bar('v6',XI,2*FX,'g');hold on;
   xd = get(hpatches,'xdata');
   yd = get(hpatches,'ydata');
   
   for n=1:size(yd,2)
     temp=max(yd(:,n));     % color value needs to be between 0 and 1
     tcolor=[abs(1-temp) abs(1-temp) abs(1-temp)];
     patch(xd(:,n),yd(:,n),tcolor);hold on;
   end

%    xlswrite(['wid_FX',num2str(i),'.xls'],FX);
   
   clear FX
   
   % ********************************************************
   
   plot([min(dat(:,i)); tmx(:,grouping); max(dat(:,i))],[0 tmy 1],'color','k','linewidth',2);hold on;
   grid off;
   axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);     
   ylabel(['Cum. Dist. RMSE'])
   xlabel('Alpha [1/s]');
   
   subplot(3,5,15) % A_TS
i=5;   
 for j=1:grouping
    tm=dat(cls*(j-1)+1:cls*j,i);
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; 
end
% calculate and plot gradients ***************************
  
   step=(max(dat(:,i))-min(dat(:,i)))/containers; % it is necessary that the boundaries are the same
   XI=[min(dat(:,i)):step:max(dat(:,i))]; % if a parameter id should be compared between two models
   
   % test for monotony
   for ii=1:5
      [II]=find(diff(tmx(1:length(tmx(:,grouping))-1,grouping))==0);
      tmx(II+1,grouping)=tmx(II+1,grouping)+[(tmx(II+2,grouping)-tmx(II+1,grouping))/2];
   end               
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping); max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
   
   [FX] = gradient(YI); % calculate gradient within containers
   
   ID_max(i)=max(FX); % keep the maximum gradient for information
   III=find(FX==max(FX));
   pos_max(i)=XI(III);
   
   ID(i,:)=FX; % keep the whole vector for each parameter
   pos(i,:)=XI;
   
   hpatches=bar('v6',XI,2*FX,'g');hold on;
   xd = get(hpatches,'xdata');
   yd = get(hpatches,'ydata');
   
   for n=1:size(yd,2)
     temp=max(yd(:,n));     % color value needs to be between 0 and 1
     tcolor=[abs(1-temp) abs(1-temp) abs(1-temp)];
     patch(xd(:,n),yd(:,n),tcolor);hold on;
   end

%    xlswrite(['wid_FX',num2str(i),'.xls'],FX);
   
   clear FX
   
   % ********************************************************
   
   plot([min(dat(:,i)); tmx(:,grouping); max(dat(:,i))],[0 tmy 1],'color','k','linewidth',2);hold on;
   grid off;
   axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);     
   ylabel(['Cum. Dist. RMSE'])
   xlabel('A_T_S [m^2]');
   
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'top 10% latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({Description1;Description2},'FontSize',14);   

set(gcf, 'WindowState', 'maximized');
saveas(stat2,[pwd '/Output_files_OTIS/OTIS_stat2.fig']);
saveas(stat2,[pwd '/Output_files_OTIS/OTIS_stat2.tif']);
close (stat2)

% Let's build the Stat1 and Stat2 images for ALL the simulations 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %                 _____ _        _   ____  
% %                / ____| |      | | |___ \ 
% %               | (___ | |_ __ _| |_  __) |
% %                \___ \| __/ _` | __||__ < 
% %                ____) | || (_| | |_ ___) |
% %               |_____/ \__\__,_|\__|____/                        
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%%%
[a,b] = sort(Order_Working_matrix(:,8)); % Sort depending on RMSE

% vector "a" puts the RMSE in crescent order while vector "b" returns a 
% the location of every RMSE value in the "a" vector

n = length(a);
nbins = 25;      %number of intervals to create the frequency plots

% thresholds corresponding to the top 20%, 10%, 5%, 1%, 0.5% and 0.01% 
% of the RMSE values --> if you want, you can change them
m = [0.2 0.1 0.05 0.01 0.005 0.001];   

% y-axis ranges
minhist = 0;
maxhist = 80;     
MAP = jet(256);

%%%%%%%%%%%%%%%%%
% Figure 1 = v, A, D, As and Alpha data vs RMSE; 
% Frequency plots and normalized cumulative distribution

stat3=figure;

subplot(3,5,1)
plot(Order_Working_matrix(:,1),Order_Working_matrix(:,8),'.k')
hold on
% [P,PP] = min(Order_Working_matrix_10(:,8)); % Ir's already sorted so it's
% of course the first line
plot(Order_Working_matrix(1,1),Order_Working_matrix(1,8),'ro')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,1)) max(Order_Working_matrix(:,1))])
ylim([0 max(Order_Working_matrix(:,8))])

subplot(3,5,2)
plot(Order_Working_matrix(:,2),Order_Working_matrix(:,8),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
hold on
plot(Order_Working_matrix(1,2),Order_Working_matrix(1,8),'ro')
xlim([min(Order_Working_matrix(:,2)) max(Order_Working_matrix(:,2))])
ylim([0 max(Order_Working_matrix(:,8))])

subplot(3,5,3)
plot(Order_Working_matrix(:,3),Order_Working_matrix(:,8),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,3)) max(Order_Working_matrix(:,3))])
ylim([0 max(Order_Working_matrix(:,8))])
hold on
plot(Order_Working_matrix(1,3),Order_Working_matrix(1,8),'ro')

subplot(3,5,4)
plot(Order_Working_matrix(:,4),Order_Working_matrix(:,8),'.k')
xlabel ('Alpha [1/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,4)) max(Order_Working_matrix(:,4))])
ylim([0 max(Order_Working_matrix(:,8))])
hold on
plot(Order_Working_matrix(1,4),Order_Working_matrix(1,8),'ro')

subplot(3,5,5)
plot(Order_Working_matrix(:,5),Order_Working_matrix(:,8),'.k')
xlabel ('A_T_S [m^2]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,5)) max(Order_Working_matrix(:,5))])
ylim([0 max(Order_Working_matrix(:,8))])
hold on
plot(Order_Working_matrix(1,5),Order_Working_matrix(1,8),'ro')

% Second row - frequency curves 

clear parameter temp p1 c1 d1 

subplot (3,5,6)
parameter=Order_Working_matrix(:,1);     % Select velocity 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end
plot([parameter(b(1)) parameter(b(1))],[0 100*max(c1)./(sum(c1))+10],'Color',[0.5 0.5 0.5],'LineWidth',1);
if sum(c1)==0
        axis([min(parameter) max(parameter) 0 100]);
    else
        axis([min(parameter) max(parameter) 0 100*max(c1)./(sum(c1))+10]);
end 
legend(num2str(m'))
xlabel ('v [m/s]');           
ylabel('Frequency')   
clear parameter temp p1 c1 d1 

subplot (3,5,7)
parameter=Order_Working_matrix(:,2);     % Select Area 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end
plot([parameter(b(1)) parameter(b(1))],[0 100*max(c1)./(sum(c1))+10],'Color',[0.5 0.5 0.5],'LineWidth',1);
if sum(c1)==0
        axis([min(parameter) max(parameter) 0 100]);
    else
        axis([min(parameter) max(parameter) 0 100*max(c1)./(sum(c1))+10]);
end 
legend(num2str(m'))
xlabel ('A [m^2]');           
ylabel('Frequency')  
clear parameter temp p1 c1 d1 

subplot (3,5,8)
parameter=Order_Working_matrix(:,3);     % Select Dispersion 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end
plot([parameter(b(1)) parameter(b(1))],[0 100*max(c1)./(sum(c1))+10],'Color',[0.5 0.5 0.5],'LineWidth',1);
if sum(c1)==0
        axis([min(parameter) max(parameter) 0 100]);
    else
        axis([min(parameter) max(parameter) 0 100*max(c1)./(sum(c1))+10]);
end 
legend(num2str(m'))
xlabel ('D [m^2/s]');           
ylabel('Frequency')  
clear parameter temp p1 c1 d1 

subplot (3,5,9)
parameter=Order_Working_matrix(:,4);     % Select Alpha 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end
plot([parameter(b(1)) parameter(b(1))],[0 100*max(c1)./(sum(c1))+10],'Color',[0.5 0.5 0.5],'LineWidth',1);
if sum(c1)==0
        axis([min(parameter) max(parameter) 0 100]);
    else
        axis([min(parameter) max(parameter) 0 100*max(c1)./(sum(c1))+10]);
end 
legend(num2str(m'))
xlabel ('Alpha [1/s]');           
ylabel('Frequency')  
clear parameter temp p1 c1 d1 

subplot (3,5,10)
parameter=Order_Working_matrix(:,5);     % Select Area transient storage 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end
plot([parameter(b(1)) parameter(b(1))],[0 100*max(c1)./(sum(c1))+10],'Color',[0.5 0.5 0.5],'LineWidth',1);
if sum(c1)==0
        axis([min(parameter) max(parameter) 0 100]);
    else
        axis([min(parameter) max(parameter) 0 100*max(c1)./(sum(c1))+10]);
end 
legend(num2str(m'))
xlabel ('A_T_S [m^2]');           
ylabel('Frequency')  
clear parameter temp p1 c1 d1 

% % % % % % % % 
% Third row - normalized cumulative (Cum. norm.) distributions of RMSE 
% lines represent the best 10% of model runs binned into 1% increments; 
% each line represents 1% of all model simulations
% % % % % % % % 
clear h of dat tmx tmy j

of=Order_Working_matrix(:,8);       % Select objective function RMSE
                               % of is automatically sorted from the lower
                               % RMSE to the higher RMSE
of=of./max(of);                % normalise of
of=1-of;                       % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods
[y,i]=sort(of);
dat=Order_Working_matrix(i,:);
cls=floor(length(dat)/10);
tmx=zeros(cls,10);
tmy=tmx;

subplot (3,5,11)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
    tm=dat(cls*(j-1)+1:cls*j,1);     % Take first column -> velocity
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls;
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['v [m/s]'])
  ylabel(['cum. norm. RMSE'])
  
colormap(jet(10));
h=axes('position',[.07 .07 .02 .33]);
h=colorbar(h);
set(h,'ytick',[2 10]);
set(h,'yticklabel',['L';'H']);
ylabel(['Likelihood RMSE'])

clear h 

subplot (3,5,12)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
    tm=dat(cls*(j-1)+1:cls*j,2);    % Take second column now -> Area
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls;
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['A [m^2]'])
  ylabel(['cum. norm. RMSE'])
  
subplot (3,5,13)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
    tm=dat(cls*(j-1)+1:cls*j,3);    % Take third column now -> Disp
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls;
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['D [m^2/s]'])
  ylabel(['cum. norm. RMSE'])
  
subplot (3,5,14)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
    tm=dat(cls*(j-1)+1:cls*j,4);    % Take fourth column now -> Alpha
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls;
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['Alpha [1/s]'])
  ylabel(['cum. norm. RMSE'])
  
  subplot (3,5,15)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
    tm=dat(cls*(j-1)+1:cls*j,5);    % Take fifth column now -> Area transient storage
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls;
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['A_T_S [m^2]'])
  ylabel(['cum. norm. RMSE'])
  
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'All latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({'OTIS latin hypercube';Description1;Description2},'FontSize',14);
  
clear h clear of dat tmx tmy j of cls 
set(gcf, 'WindowState', 'maximized');

saveas(stat3,[pwd '/Output_files_OTIS/OTIS_stat3.fig']);
saveas(stat3,[pwd '/Output_files_OTIS/OTIS_stat3.tif']);
close (stat3);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %                   _____ _        _   _  _   
% %                  / ____| |      | | | || |  
% %                 | (___ | |_ __ _| |_| || |_ 
% %                  \___ \| __/ _` | __|__   _|
% %                  ____) | || (_| | |_   | |  
% %                 |_____/ \__\__,_|\__|  |_|                     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Figure 4 = v, A and D data vs RMSE; A posteriori parameter distribution;
% parameter identifiability for ALL the results

stat4=figure;

subplot(3,5,1)
plot(Order_Working_matrix(:,1),Order_Working_matrix(:,8),'.k')
hold on
% [P,PP] = min(Order_Working_matrix_10(:,8)); % Ir's already sorted so it's
% of course the first line
plot(Order_Working_matrix(1,1),Order_Working_matrix(1,8),'ro')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,1)) max(Order_Working_matrix(:,1))])
ylim([0 max(Order_Working_matrix(:,8))])

subplot(3,5,2)
plot(Order_Working_matrix(:,2),Order_Working_matrix(:,8),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
hold on
plot(Order_Working_matrix(1,2),Order_Working_matrix(1,8),'ro')
xlim([min(Order_Working_matrix(:,2)) max(Order_Working_matrix(:,2))])
ylim([0 max(Order_Working_matrix(:,8))])

subplot(3,5,3)
plot(Order_Working_matrix(:,3),Order_Working_matrix(:,8),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,3)) max(Order_Working_matrix(:,3))])
ylim([0 max(Order_Working_matrix(:,8))])
hold on
plot(Order_Working_matrix(1,3),Order_Working_matrix(1,8),'ro')

subplot(3,5,4)
plot(Order_Working_matrix(:,4),Order_Working_matrix(:,8),'.k')
xlabel ('Alpha [1/s]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,4)) max(Order_Working_matrix(:,4))])
ylim([0 max(Order_Working_matrix(:,8))])
hold on
plot(Order_Working_matrix(1,4),Order_Working_matrix(1,8),'ro')

subplot(3,5,5)
plot(Order_Working_matrix(:,5),Order_Working_matrix(:,8),'.k')
xlabel ('A_T_S [m^2]');
ylabel ('RMSE');
xlim([min(Order_Working_matrix(:,5)) max(Order_Working_matrix(:,5))])
ylim([0 max(Order_Working_matrix(:,8))])
hold on
plot(Order_Working_matrix(1,5),Order_Working_matrix(1,8),'ro')

% % % % % % % % 
% ROW 2 -> A POSTERIORI PARAMETER DISTRIBUTION
% % % % % % % % 

ncontainers=20;
slider_value=100;   

% calculate likelihood
of=Order_Working_matrix(:,8);  % select the criterion -> We've already sorted from the lower to the higher RMSE
                                  % RMSE -> low values indicate better models
                                
LL=1-of; 						% likelihood (high values indicate more likely [probable] models)
if min(LL)<0|min(LL)==0, LL=LL-min(LL)+1000*eps;end; % transform negative lhoods

LL=LL./sum(LL);     % sum(likelihoods)=1      problem if NaN in vector

%Eliminate data that is below the slider threshold
LL=sortrows(LL);
dat=sortrows(Order_Working_matrix,8);     % sort Hyperspace based on the RMSE values
LL = flipud(LL);
numdat = floor((slider_value / 100) * size(LL));
LL(numdat+1:size(LL),:) = [];
dat(numdat+1:size(dat),:) = [];

subplot(3,5,6) % velocity

width=(max(dat(:,1))-min(dat(:,1)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,1)>min(dat(:,1))+(n-1)*width&dat(:,1)<min(dat(:,1))+n*width); % find all values within the container
      temp1=dat(K,1);
      temp2=LL(K);
      
      bx(1,n)=min(dat(:,1))+.5*width+(n-1)*width;   % middle of container
      by(1,n)=sum(temp2)/length(temp2);             % mean likelihood in container
      
      clear temp 1 temp2 K
end 
bar(bx(1,:),by(1,:),'k');hold on;
ymax=(max(by(1,:))+0.1*(max(by(1,:))));
ymin=(min(by(1,:))-0.1*(min(by(1,:))));
axis([min(dat(:,1)) max(dat(:,1)) ymin ymax]);
ylabel(['D (RMSE)'])
xlabel('v [m/s]');
   
subplot(3,5,7) % Area

width=(max(dat(:,2))-min(dat(:,2)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,2)>min(dat(:,2))+(n-1)*width&dat(:,2)<min(dat(:,2))+n*width); % find all values within the container
      temp1=dat(K,2);
      temp2=LL(K);
      
      bx(2,n)=min(dat(:,2))+.5*width+(n-1)*width;   % middle of container
      by(2,n)=sum(temp2)/length(temp2);             % mean likelihood in container
      
      clear temp 1 temp2 K
end 
bar(bx(2,:),by(2,:),'k');hold on;
ymax=(max(by(2,:))+0.1*(max(by(2,:))));
ymin=(min(by(2,:))-0.1*(min(by(2,:))));
axis([min(dat(:,2)) max(dat(:,2)) ymin ymax]);
ylabel(['D (RMSE)'])
xlabel('A [m^2]');

subplot(3,5,8) % Disp

width=(max(dat(:,3))-min(dat(:,3)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,3)>min(dat(:,3))+(n-1)*width&dat(:,3)<min(dat(:,3))+n*width); % find all values within the container
      temp1=dat(K,3);
      temp2=LL(K);
      
      bx(3,n)=min(dat(:,3))+.5*width+(n-1)*width;   % middle of container
      by(3,n)=sum(temp2)/length(temp2);             % mean likelihood in container
      
      clear temp 1 temp2 K
end 
bar(bx(3,:),by(3,:),'k');hold on;
ymax=(max(by(3,:))+0.1*(max(by(3,:))));
ymin=(min(by(3,:))-0.1*(min(by(3,:))));
axis([min(dat(:,3)) max(dat(:,3)) ymin ymax]);
ylabel(['D (RMSE)'])
xlabel('D [m^2/s]');

subplot(3,5,9) % Alpha

width=(max(dat(:,4))-min(dat(:,4)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,4)>min(dat(:,4))+(n-1)*width&dat(:,4)<min(dat(:,4))+n*width); % find all values within the container
      temp1=dat(K,4);
      temp2=LL(K);
      
      bx(4,n)=min(dat(:,4))+.5*width+(n-1)*width;   % middle of container
      by(4,n)=sum(temp2)/length(temp2);             % mean likelihood in container
      
      clear temp 1 temp2 K
end 
bar(bx(4,:),by(4,:),'k');hold on;
ymax=(max(by(4,:))+0.1*(max(by(4,:))));
ymin=(min(by(4,:))-0.1*(min(by(4,:))));
axis([min(dat(:,4)) max(dat(:,4)) ymin ymax]);
ylabel(['D (RMSE)'])
xlabel('Alpha [1/s]');

subplot(3,5,10) % A_trans storage

width=(max(dat(:,5))-min(dat(:,5)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,5)>min(dat(:,5))+(n-1)*width&dat(:,5)<min(dat(:,5))+n*width); % find all values within the container
      temp1=dat(K,5);
      temp2=LL(K);
      
      bx(5,n)=min(dat(:,5))+.5*width+(n-1)*width;   % middle of container
      by(5,n)=sum(temp2)/length(temp2);             % mean likelihood in container
      
      clear temp 1 temp2 K
end 
bar(bx(5,:),by(5,:),'k');hold on;
ymax=(max(by(5,:))+0.1*(max(by(5,:))));
ymin=(min(by(5,:))-0.1*(min(by(5,:))));
axis([min(dat(:,5)) max(dat(:,5)) ymin ymax]);
ylabel(['D (RMSE)'])
xlabel('A_T_S [m^2]');

% % % % % % % % 
% 3rd ROW - IDENTIFIABILITY PLOTS
% % % % % % % % 
clear of I J dat numdat cls tmx i
dat=Order_Working_matrix;   % dat values sorted from the lower RMSE to the higher RMSE
containers=10; % division of each parameter range
grouping=10;  % number of groups

% calculate likelihood
of=Order_Working_matrix(:,8);  % criteria (low values indicate better models)
of=1-of; % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods

% sort data according to selected perf
[I,J]=sort(of); 
dat=dat(J,:);
% Eliminate data that is below the slider threshold
numdat = floor((slider_value / 100) * size(dat));
dat(numdat+1:size(dat),:) = [];

cls=floor(length(dat)/grouping);
tmx=zeros(cls,grouping);tmy=tmx;

subplot(3,5,11) % velocity
i=1;
for j=1:grouping
    tm=dat(cls*(j-1)+1:cls*j,i);
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; 
end
% calculate and plot gradients ***************************
  
   step=(max(dat(:,i))-min(dat(:,i)))/containers; % it is necessary that the boundaries are the same
   XI=[min(dat(:,i)):step:max(dat(:,i))]; % if a parameter id should be compared between two models
   
   % test for monotony
   for ii=1:5
      [II]=find(diff(tmx(1:length(tmx(:,grouping))-1,grouping))==0);
      tmx(II+1,grouping)=tmx(II+1,grouping)+[(tmx(II+2,grouping)-tmx(II+1,grouping))/2];
   end               
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping); max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
   
   [FX] = gradient(YI); % calculate gradient within containers
   
   ID_max(i)=max(FX); % keep the maximum gradient for information
   III=find(FX==max(FX));
   pos_max(i)=XI(III);
   
   ID(i,:)=FX; % keep the whole vector for each parameter
   pos(i,:)=XI;
   
   hpatches=bar('v6',XI,2*FX,'g');hold on;
   xd = get(hpatches,'xdata');
   yd = get(hpatches,'ydata');
   
   for n=1:size(yd,2)
     temp=max(yd(:,n));     % color value needs to be between 0 and 1
     tcolor=[abs(1-temp) abs(1-temp) abs(1-temp)];
     patch(xd(:,n),yd(:,n),tcolor);hold on;
   end

%    xlswrite(['wid_FX',num2str(i),'.xls'],FX);
   
   clear FX
   
   % ********************************************************
   
   plot([min(dat(:,i)); tmx(:,grouping); max(dat(:,i))],[0 tmy 1],'color','k','linewidth',2);hold on;
   grid off;
   axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
   ylabel(['Cum. Dist. RMSE'])
   xlabel('v [m/s]');
    
subplot(3,5,12) % Area
i=2;   
 for j=1:grouping
    tm=dat(cls*(j-1)+1:cls*j,i);
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; 
end
% calculate and plot gradients ***************************
  
   step=(max(dat(:,i))-min(dat(:,i)))/containers; % it is necessary that the boundaries are the same
   XI=[min(dat(:,i)):step:max(dat(:,i))]; % if a parameter id should be compared between two models
   
   % test for monotony
   for ii=1:5
      [II]=find(diff(tmx(1:length(tmx(:,grouping))-1,grouping))==0);
      tmx(II+1,grouping)=tmx(II+1,grouping)+[(tmx(II+2,grouping)-tmx(II+1,grouping))/2];
   end               
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping); max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
   
   [FX] = gradient(YI); % calculate gradient within containers
   
   ID_max(i)=max(FX); % keep the maximum gradient for information
   III=find(FX==max(FX));
   pos_max(i)=XI(III);
   
   ID(i,:)=FX; % keep the whole vector for each parameter
   pos(i,:)=XI;
   
   hpatches=bar('v6',XI,2*FX,'g');hold on;
   xd = get(hpatches,'xdata');
   yd = get(hpatches,'ydata');
   
   for n=1:size(yd,2)
     temp=max(yd(:,n));     % color value needs to be between 0 and 1
     tcolor=[abs(1-temp) abs(1-temp) abs(1-temp)];
     patch(xd(:,n),yd(:,n),tcolor);hold on;
   end

%    xlswrite(['wid_FX',num2str(i),'.xls'],FX);
   
   clear FX
   
   % ********************************************************
   
   plot([min(dat(:,i)); tmx(:,grouping); max(dat(:,i))],[0 tmy 1],'color','k','linewidth',2);hold on;
   grid off;
   axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);  
   ylabel(['Cum. Dist. RMSE'])
   xlabel('A [m^2]');
   
subplot(3,5,13) % Disp
i=3;   
 for j=1:grouping
    tm=dat(cls*(j-1)+1:cls*j,i);
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; 
end
% calculate and plot gradients ***************************
  
   step=(max(dat(:,i))-min(dat(:,i)))/containers; % it is necessary that the boundaries are the same
   XI=[min(dat(:,i)):step:max(dat(:,i))]; % if a parameter id should be compared between two models
   
   % test for monotony
   for ii=1:5
      [II]=find(diff(tmx(1:length(tmx(:,grouping))-1,grouping))==0);
      tmx(II+1,grouping)=tmx(II+1,grouping)+[(tmx(II+2,grouping)-tmx(II+1,grouping))/2];
   end               
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping); max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
   
   [FX] = gradient(YI); % calculate gradient within containers
   
   ID_max(i)=max(FX); % keep the maximum gradient for information
   III=find(FX==max(FX));
   pos_max(i)=XI(III);
   
   ID(i,:)=FX; % keep the whole vector for each parameter
   pos(i,:)=XI;
   
   hpatches=bar('v6',XI,2*FX,'g');hold on;
   xd = get(hpatches,'xdata');
   yd = get(hpatches,'ydata');
   
   for n=1:size(yd,2)
     temp=max(yd(:,n));     % color value needs to be between 0 and 1
     tcolor=[abs(1-temp) abs(1-temp) abs(1-temp)];
     patch(xd(:,n),yd(:,n),tcolor);hold on;
   end

%    xlswrite(['wid_FX',num2str(i),'.xls'],FX);
   
   clear FX
   
   % ********************************************************
   
   plot([min(dat(:,i)); tmx(:,grouping); max(dat(:,i))],[0 tmy 1],'color','k','linewidth',2);hold on;
   grid off;
   axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);     
   ylabel(['Cum. Dist. RMSE'])
   xlabel('D [m^2/s]');
   
   subplot(3,5,14) % Alpha
i=4;   
 for j=1:grouping
    tm=dat(cls*(j-1)+1:cls*j,i);
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; 
end
% calculate and plot gradients ***************************
  
   step=(max(dat(:,i))-min(dat(:,i)))/containers; % it is necessary that the boundaries are the same
   XI=[min(dat(:,i)):step:max(dat(:,i))]; % if a parameter id should be compared between two models
   
   % test for monotony
   for ii=1:5
      [II]=find(diff(tmx(1:length(tmx(:,grouping))-1,grouping))==0);
      tmx(II+1,grouping)=tmx(II+1,grouping)+[(tmx(II+2,grouping)-tmx(II+1,grouping))/2];
   end               
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping); max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
   
   [FX] = gradient(YI); % calculate gradient within containers
   
   ID_max(i)=max(FX); % keep the maximum gradient for information
   III=find(FX==max(FX));
   pos_max(i)=XI(III);
   
   ID(i,:)=FX; % keep the whole vector for each parameter
   pos(i,:)=XI;
   
   hpatches=bar('v6',XI,2*FX,'g');hold on;
   xd = get(hpatches,'xdata');
   yd = get(hpatches,'ydata');
   
   for n=1:size(yd,2)
     temp=max(yd(:,n));     % color value needs to be between 0 and 1
     tcolor=[abs(1-temp) abs(1-temp) abs(1-temp)];
     patch(xd(:,n),yd(:,n),tcolor);hold on;
   end

%    xlswrite(['wid_FX',num2str(i),'.xls'],FX);
   
   clear FX
   
   % ********************************************************
   
   plot([min(dat(:,i)); tmx(:,grouping); max(dat(:,i))],[0 tmy 1],'color','k','linewidth',2);hold on;
   grid off;
   axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);     
   ylabel(['Cum. Dist. RMSE'])
   xlabel('Alpha [1/s]');
   
   subplot(3,5,15) % A_TS
i=5;   
 for j=1:grouping
    tm=dat(cls*(j-1)+1:cls*j,i);
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; 
end
% calculate and plot gradients ***************************
  
   step=(max(dat(:,i))-min(dat(:,i)))/containers; % it is necessary that the boundaries are the same
   XI=[min(dat(:,i)):step:max(dat(:,i))]; % if a parameter id should be compared between two models
   
   % test for monotony
   for ii=1:5
      [II]=find(diff(tmx(1:length(tmx(:,grouping))-1,grouping))==0);
      tmx(II+1,grouping)=tmx(II+1,grouping)+[(tmx(II+2,grouping)-tmx(II+1,grouping))/2];
   end               
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping); max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
   
   [FX] = gradient(YI); % calculate gradient within containers
   
   ID_max(i)=max(FX); % keep the maximum gradient for information
   III=find(FX==max(FX));
   pos_max(i)=XI(III);
   
   ID(i,:)=FX; % keep the whole vector for each parameter
   pos(i,:)=XI;
   
   hpatches=bar('v6',XI,2*FX,'g');hold on;
   xd = get(hpatches,'xdata');
   yd = get(hpatches,'ydata');
   
   for n=1:size(yd,2)
     temp=max(yd(:,n));     % color value needs to be between 0 and 1
     tcolor=[abs(1-temp) abs(1-temp) abs(1-temp)];
     patch(xd(:,n),yd(:,n),tcolor);hold on;
   end

%    xlswrite(['wid_FX',num2str(i),'.xls'],FX);
   
   clear FX
   
   % ********************************************************
   
   plot([min(dat(:,i)); tmx(:,grouping); max(dat(:,i))],[0 tmy 1],'color','k','linewidth',2);hold on;
   grid off;
   axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);     
   ylabel(['Cum. Dist. RMSE'])
   xlabel('A_T_S [m^2]');
   
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'All latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({Description1;Description2},'FontSize',14);   

set(gcf, 'WindowState', 'maximized');
saveas(stat4,[pwd '/Output_files_OTIS/OTIS_stat4.fig']);
saveas(stat4,[pwd '/Output_files_OTIS/OTIS_stat4.tif']);
close (stat4)









