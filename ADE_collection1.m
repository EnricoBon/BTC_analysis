function [ADE_1] = ADE_collection1(ADE1,ADE2,ADE3,Interp_curve_01,time,BTC_input,Description1,Description2);

% In this figure we have 3 plots, in the first one we plot the results of 
% the 1st ADE with the best-fitting (according to RMSE) curve. Same for the
% second plot (ADE2). In the third plot the best results (RMSE) after the
% latin hypercube sampling with the the top 1% of the results 

ADE_1=figure;

formatSpec1="v_m_i_n=%0.4f m/s";
formatSpec2="v_m_a_x=%0.4f m/s";
formatSpec3="A_m_i_n=%0.4f m^2";
formatSpec4="A_m_a_x=%0.4f m^2";
formatSpec5="D_m_i_n=%0.4f m^2/s"; 
formatSpec6="D_m_a_x=%0.4f m^2/s";
formatSpec7="RMSE_m_i_n=%0.4f"; 
formatSpec8="RMSE_m_a_x=%0.4f";
formatSpec9="r^2_m_i_n=%0.4f"; 
formatSpec10="r^2_m_a_x=%0.4f";
str(1,1)=sprintf(formatSpec1,cell2mat(ADE3.Hyperspace_01_Limits(2,1)));
str(2,1)=sprintf(formatSpec2,cell2mat(ADE3.Hyperspace_01_Limits(2,2)));
str(3,1)=sprintf(formatSpec3,cell2mat(ADE3.Hyperspace_01_Limits(2,3)));
str(4,1)=sprintf(formatSpec4,cell2mat(ADE3.Hyperspace_01_Limits(2,4)));
str(5,1)=sprintf(formatSpec5,cell2mat(ADE3.Hyperspace_01_Limits(2,5)));
str(6,1)=sprintf(formatSpec6,cell2mat(ADE3.Hyperspace_01_Limits(2,6)));
str(7,1)=sprintf(formatSpec7,cell2mat(ADE3.Hyperspace_01_Limits(2,7)));
str(8,1)=sprintf(formatSpec8,cell2mat(ADE3.Hyperspace_01_Limits(2,8)));
str(9,1)=sprintf(formatSpec9,cell2mat(ADE3.Hyperspace_01_Limits(2,9)));
str(10,1)=sprintf(formatSpec10,cell2mat(ADE3.Hyperspace_01_Limits(2,10)));

subplot (1,3,3)
plot(Interp_curve_01(:,1), Interp_curve_01(:,2), 'k', 'LineWidth', 1,'HandleVisibility','off');
hold on;
plot(Interp_curve_01(:,1), Interp_curve_01(:,3), 'k', 'LineWidth', 1,'HandleVisibility','off');
x2 = [Interp_curve_01(:,1)', fliplr(Interp_curve_01(:,1)')];
inBetween = [Interp_curve_01(:,2)', fliplr(Interp_curve_01(:,3)')];
fill(x2, inBetween, [0.86,0.86,0.86]);
hold on
plot(time(:,1),BTC_input(:,2),'-r','LineWidth',2)
legend('Top 0.1% ADE results','Observed BTC')
xlabel ('time [s]');
ylabel ('Cl Conc [mg/l]');
annotation('textbox', [0.789, 0.67, 0.1, 0.1], 'String',str,'FontSize',11,'LineStyle','-')
% annotation('textbox', [0.75, 0.9, 0.1, 0.1], 'String', {"Top 0.1% ADE results from Latin","Hypercube sampling (100'000)"},...
%     'FontSize',12,'LineStyle','none','FitBoxToText','on','HorizontalAlignment','center')
title({"Top 0.1% ADE results from Latin","Hypercube sampling"},'FontSize',12,'LineStyle','none')

clear formatSpec1 formatSpec2 formatSpec3 formatSpec4 formatSpec5 formatSpec6
clear formatSpec7 formatSpec8 formatSpec9 formatSpec10 str

formatSpec1="v=%0.4f m/s";
formatSpec2="A=%0.4f m^2";
formatSpec3="D=%0.4f m^2/s";
formatSpec4="M=%0.4f g";
formatSpec5="RMSE=%0.4f";
formatSpec6="r^2=%0.4f";
str(1,1)=sprintf(formatSpec1,(ADE1.v(1,1)));
str(2,1)=sprintf(formatSpec2,(ADE1.A(1,1)));
str(3,1)=sprintf(formatSpec3,cell2mat(ADE1.Best(2,1)));
str(4,1)=sprintf(formatSpec4,ADE1.InjectedMass(1,1));
str(5,1)=sprintf(formatSpec5,cell2mat(ADE1.Best(2,2)));
str(6,1)=sprintf(formatSpec6,cell2mat(ADE1.Best(2,4)));

subplot (1,3,1)
plot(time(:,1),BTC_input(:,2),'-r','LineWidth',2)   % Observed curve
hold on
plot(time(:,1),ADE1.best_BTC(:,2),'-k','LineWidth',2)    % Best-fitting RSS ADE
legend('Observed BTC','Calibrated ADE')
annotation('textbox',[.25 .67 .1 .1],'String',str,'FitBoxToText','on');
title({"Best-fitting ADE from","dilution gauging"},'FontSize',12,'LineStyle','none')
xlabel ('Time [s]');

clear formatSpec1 formatSpec2 formatSpec3 formatSpec4 formatSpec5 formatSpec6 str
formatSpec1="v=%0.4f m/s";
formatSpec2="A=%0.4f m^2";
formatSpec3="D=%0.4f m^2/s";
formatSpec4="M=%0.4f g";
formatSpec5="RMSE=%0.4f";
formatSpec6="r^2=%0.4f";
str(1,1)=sprintf(formatSpec1,(ADE2.v(1,1)));
str(2,1)=sprintf(formatSpec2,(ADE2.A(1,1)));
str(3,1)=sprintf(formatSpec3,cell2mat(ADE2.Best(2,1)));
str(4,1)=sprintf(formatSpec4,ADE2.RecoveredMass(1,1));
str(5,1)=sprintf(formatSpec5,cell2mat(ADE2.Best(2,2)));
str(6,1)=sprintf(formatSpec6,cell2mat(ADE2.Best(2,4)));

subplot (1,3,2)
plot(time(:,1),BTC_input(:,2),'-r','LineWidth',2)   % Observed curve
hold on
plot(time(:,1),ADE2.best_BTC(:,2),'-k','LineWidth',2)    % Best-fitting RSS ADE
legend('Observed BTC','Calibrated ADE')
annotation('textbox',[.53 .67 .1 .1],'String',str,'FitBoxToText','on');
title({"Best-fitting ADE by","using measured discharge"},'FontSize',12,'LineStyle','none')
xlabel ('Time [s]');

sgtitle({Description1;Description2},'FontSize',14);


annotation('textbox', [0.12, 0.87, 0.1, 0.1],'String', "RMSE and r^2 (NSE)",...
    'FontSize',12,'LineStyle','none','FitBoxToText','on')

set(gcf, 'WindowState', 'maximized');
end

