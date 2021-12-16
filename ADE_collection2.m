function [ADE_2,ADE3_limits] = ADE_collection2(Hyperspace, time, Length, BTC_input,M_g,L)

top01=length(Hyperspace(:,1))*0.001;      % Limit for top 0.1% of the results

Hyperspace_temp1=sortrows(Hyperspace,4);   % Hyperspace sorted depending on RMSE (or r^2 (NSE) values
Hyperspace_temp2=sortrows(Hyperspace,7);   % Hyperspace sorted depending on logRMSE values
Hyperspace_temp3=sortrows(Hyperspace,11,'descend');   % Hyperspace sorted depending on KGE values

% Select only the top 1% for each one
Hyperspace_temp1=Hyperspace_temp1(1:top01,:);
Hyperspace_temp2=Hyperspace_temp2(1:top01,:);
Hyperspace_temp3=Hyperspace_temp3(1:top01,:);

% ADE collection for top 0.1 of the RMSE values
for k=1:1:length(Hyperspace_temp1(:,1))
for i=1:1:Length
        ADE_temp1(i,k)=(M_g/((Hyperspace_temp1(k,2)*(4*3.14*Hyperspace_temp1(k,3)*time(i,1))^(1/2))))*...
            exp(-((L-Hyperspace_temp1(k,1)*time(i,1))^2/(4*Hyperspace_temp1(k,3)*time(i,1))));
end
end

% ADE collection for top 0.1 of the logRMSE values
for k=1:1:length(Hyperspace_temp2(:,1))
for i=1:1:Length
        ADE_temp2(i,k)=(M_g/((Hyperspace_temp2(k,2)*(4*3.14*Hyperspace_temp2(k,3)*time(i,1))^(1/2))))*...
            exp(-((L-Hyperspace_temp2(k,1)*time(i,1))^2/(4*Hyperspace_temp2(k,3)*time(i,1))));
end
end

% ADE collection for top 0.1 of the KGE values
for k=1:1:length(Hyperspace_temp3(:,1))
for i=1:1:Length
        ADE_temp3(i,k)=(M_g/((Hyperspace_temp3(k,2)*(4*3.14*Hyperspace_temp3(k,3)*time(i,1))^(1/2))))*...
            exp(-((L-Hyperspace_temp3(k,1)*time(i,1))^2/(4*Hyperspace_temp3(k,3)*time(i,1))));
end
end

for i=1:1:Length
    Interp_curve_temp1(i,1)=time(i,1);
    Interp_curve_temp1(i,2)=min(ADE_temp1(i,:));
    Interp_curve_temp1(i,3)=max(ADE_temp1(i,:));
    
    Interp_curve_temp2(i,1)=time(i,1);
    Interp_curve_temp2(i,2)=min(ADE_temp2(i,:));
    Interp_curve_temp2(i,3)=max(ADE_temp2(i,:));
    
    Interp_curve_temp3(i,1)=time(i,1);
    Interp_curve_temp3(i,2)=min(ADE_temp3(i,:));
    Interp_curve_temp3(i,3)=max(ADE_temp3(i,:));
end

clear ADE_temp1 ADE_temp2 ADE_temp3

ADE3_limits.RMSE(1,1)={'v min'};          
ADE3_limits.RMSE(2,1)=num2cell(min(Hyperspace_temp1(:,1)));
ADE3_limits.RMSE(1,2)={'v max'};          
ADE3_limits.RMSE(2,2)=num2cell(max(Hyperspace_temp1(:,1)));
ADE3_limits.RMSE(1,3)={'A min'};          
ADE3_limits.RMSE(2,3)=num2cell(min(Hyperspace_temp1(:,2)));
ADE3_limits.RMSE(1,4)={'A max'};          
ADE3_limits.RMSE(2,4)=num2cell(max(Hyperspace_temp1(:,2)));
ADE3_limits.RMSE(1,5)={'D min'};          
ADE3_limits.RMSE(2,5)=num2cell(min(Hyperspace_temp1(:,3)));
ADE3_limits.RMSE(1,6)={'D max'};          
ADE3_limits.RMSE(2,6)=num2cell(max(Hyperspace_temp1(:,3)));
ADE3_limits.RMSE(1,7)={'RMSE min'};          
ADE3_limits.RMSE(2,7)=num2cell(min(Hyperspace_temp1(:,4)));
ADE3_limits.RMSE(1,8)={'RMSE max'};          
ADE3_limits.RMSE(2,8)=num2cell(max(Hyperspace_temp1(:,4)));
ADE3_limits.RMSE(1,9)={'r^2 (NSE) min'};          
ADE3_limits.RMSE(2,9)=num2cell(min(Hyperspace_temp1(:,5)));
ADE3_limits.RMSE(1,10)={'r^2 (NSE) max'};          
ADE3_limits.RMSE(2,10)=num2cell(max(Hyperspace_temp1(:,5)));

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
str(1,1)=sprintf(formatSpec1,cell2mat(ADE3_limits.RMSE(2,1)));
str(2,1)=sprintf(formatSpec2,cell2mat(ADE3_limits.RMSE(2,2)));
str(3,1)=sprintf(formatSpec3,cell2mat(ADE3_limits.RMSE(2,3)));
str(4,1)=sprintf(formatSpec4,cell2mat(ADE3_limits.RMSE(2,4)));
str(5,1)=sprintf(formatSpec5,cell2mat(ADE3_limits.RMSE(2,5)));
str(6,1)=sprintf(formatSpec6,cell2mat(ADE3_limits.RMSE(2,6)));
str(7,1)=sprintf(formatSpec7,cell2mat(ADE3_limits.RMSE(2,7)));
str(8,1)=sprintf(formatSpec8,cell2mat(ADE3_limits.RMSE(2,8)));
str(9,1)=sprintf(formatSpec9,cell2mat(ADE3_limits.RMSE(2,9)));
str(10,1)=sprintf(formatSpec10,cell2mat(ADE3_limits.RMSE(2,10)));

ADE_2=figure;

%%%% RMSE and r^2
subplot (1,3,1)
plot(Interp_curve_temp1(:,1), Interp_curve_temp1(:,2), 'k', 'LineWidth', 1,'HandleVisibility','off');
hold on;
plot(Interp_curve_temp1(:,1), Interp_curve_temp1(:,3), 'k', 'LineWidth', 1,'HandleVisibility','off');
x2 = [Interp_curve_temp1(:,1)', fliplr(Interp_curve_temp1(:,1)')];
inBetween = [Interp_curve_temp1(:,2)', fliplr(Interp_curve_temp1(:,3)')];
fill(x2, inBetween, [0.86,0.86,0.86]);
hold on
plot(time(:,1),BTC_input(:,2),'-r','LineWidth',2)
legend('Top 0.1% ADE results','Observed BTC')
xlabel ('time [s]');
ylabel ('Cl Conc [mg/l]');
annotation('textbox', [0.23, 0.73, 0.1, 0.1], 'String',str,'FontSize',11,'LineStyle','-')
% annotation('textbox', [0.75, 0.9, 0.1, 0.1], 'String', {"Top 0.1% ADE results from Latin","Hypercube sampling (100'000)"},...
%     'FontSize',12,'LineStyle','none','FitBoxToText','on','HorizontalAlignment','center')
title({"Top 0.1% ADE for RMSE (and NSE) results","from Latin Hypercube sampling"},'FontSize',12,'LineStyle','none')

clear formatSpec1 formatSpec2 formatSpec3 formatSpec4 formatSpec5 formatSpec6
clear formatSpec7 formatSpec8 formatSpec9 formatSpec10 str


%%%% logRMSE

ADE3_limits.logRMSE(1,1)={'v min'};          
ADE3_limits.logRMSE(2,1)=num2cell(min(Hyperspace_temp2(:,1)));
ADE3_limits.logRMSE(1,2)={'v max'};          
ADE3_limits.logRMSE(2,2)=num2cell(max(Hyperspace_temp2(:,1)));
ADE3_limits.logRMSE(1,3)={'A min'};          
ADE3_limits.logRMSE(2,3)=num2cell(min(Hyperspace_temp2(:,2)));
ADE3_limits.logRMSE(1,4)={'A max'};          
ADE3_limits.logRMSE(2,4)=num2cell(max(Hyperspace_temp2(:,2)));
ADE3_limits.logRMSE(1,5)={'D min'};          
ADE3_limits.logRMSE(2,5)=num2cell(min(Hyperspace_temp2(:,3)));
ADE3_limits.logRMSE(1,6)={'D max'};          
ADE3_limits.logRMSE(2,6)=num2cell(max(Hyperspace_temp2(:,3)));
ADE3_limits.logRMSE(1,7)={'logRMSE min'};          
ADE3_limits.logRMSE(2,7)=num2cell(min(Hyperspace_temp2(:,7)));
ADE3_limits.logRMSE(1,8)={'logRMSE max'};          
ADE3_limits.logRMSE(2,8)=num2cell(max(Hyperspace_temp2(:,7)));

formatSpec1="v_m_i_n=%0.4f m/s";
formatSpec2="v_m_a_x=%0.4f m/s";
formatSpec3="A_m_i_n=%0.4f m^2";
formatSpec4="A_m_a_x=%0.4f m^2";
formatSpec5="D_m_i_n=%0.4f m^2/s"; 
formatSpec6="D_m_a_x=%0.4f m^2/s";
formatSpec7="logRMSE_m_i_n=%0.4f"; 
formatSpec8="logRMSE_m_a_x=%0.4f";

str(1,1)=sprintf(formatSpec1,cell2mat(ADE3_limits.logRMSE(2,1)));
str(2,1)=sprintf(formatSpec2,cell2mat(ADE3_limits.logRMSE(2,2)));
str(3,1)=sprintf(formatSpec3,cell2mat(ADE3_limits.logRMSE(2,3)));
str(4,1)=sprintf(formatSpec4,cell2mat(ADE3_limits.logRMSE(2,4)));
str(5,1)=sprintf(formatSpec5,cell2mat(ADE3_limits.logRMSE(2,5)));
str(6,1)=sprintf(formatSpec6,cell2mat(ADE3_limits.logRMSE(2,6)));
str(7,1)=sprintf(formatSpec7,cell2mat(ADE3_limits.logRMSE(2,7)));
str(8,1)=sprintf(formatSpec8,cell2mat(ADE3_limits.logRMSE(2,8)));

subplot (1,3,2)
plot(Interp_curve_temp2(:,1), Interp_curve_temp2(:,2), 'k', 'LineWidth', 1,'HandleVisibility','off');
hold on;
plot(Interp_curve_temp2(:,1), Interp_curve_temp2(:,3), 'k', 'LineWidth', 1,'HandleVisibility','off');
x2 = [Interp_curve_temp2(:,1)', fliplr(Interp_curve_temp2(:,1)')];
inBetween = [Interp_curve_temp2(:,2)', fliplr(Interp_curve_temp2(:,3)')];
fill(x2, inBetween, [0.86,0.86,0.86]);
hold on
plot(time(:,1),BTC_input(:,2),'-r','LineWidth',2)
legend('Top 0.1% ADE results','Observed BTC')
xlabel ('time [s]');
ylabel ('Cl Conc [mg/l]');
annotation('textbox', [0.51, 0.73, 0.1, 0.1], 'String',str,'FontSize',11,'LineStyle','-')
% annotation('textbox', [0.75, 0.9, 0.1, 0.1], 'String', {"Top 0.1% ADE results from Latin","Hypercube sampling (100'000)"},...
%     'FontSize',12,'LineStyle','none','FitBoxToText','on','HorizontalAlignment','center')
title({"Top 0.1% ADE for logRMSE results","from Latin Hypercube sampling"},'FontSize',12,'LineStyle','none')

clear formatSpec1 formatSpec2 formatSpec3 formatSpec4 formatSpec5 formatSpec6
clear formatSpec7 formatSpec8 formatSpec9 formatSpec10 str


%%%% KGE

ADE3_limits.KGE(1,1)={'v min'};          
ADE3_limits.KGE(2,1)=num2cell(min(Hyperspace_temp3(:,1)));
ADE3_limits.KGE(1,2)={'v max'};          
ADE3_limits.KGE(2,2)=num2cell(max(Hyperspace_temp3(:,1)));
ADE3_limits.KGE(1,3)={'A min'};          
ADE3_limits.KGE(2,3)=num2cell(min(Hyperspace_temp3(:,2)));
ADE3_limits.KGE(1,4)={'A max'};          
ADE3_limits.KGE(2,4)=num2cell(max(Hyperspace_temp3(:,2)));
ADE3_limits.KGE(1,5)={'D min'};          
ADE3_limits.KGE(2,5)=num2cell(min(Hyperspace_temp3(:,3)));
ADE3_limits.KGE(1,6)={'D max'};          
ADE3_limits.KGE(2,6)=num2cell(max(Hyperspace_temp3(:,3)));
ADE3_limits.KGE(1,7)={'KGE min'};          
ADE3_limits.KGE(2,7)=num2cell(min(Hyperspace_temp3(:,11)));
ADE3_limits.KGE(1,8)={'KGE max'};          
ADE3_limits.KGE(2,8)=num2cell(max(Hyperspace_temp3(:,11)));

formatSpec1="v_m_i_n=%0.4f m/s";
formatSpec2="v_m_a_x=%0.4f m/s";
formatSpec3="A_m_i_n=%0.4f m^2";
formatSpec4="A_m_a_x=%0.4f m^2";
formatSpec5="D_m_i_n=%0.4f m^2/s"; 
formatSpec6="D_m_a_x=%0.4f m^2/s";
formatSpec7="KGE_m_i_n=%0.4f"; 
formatSpec8="KGE_m_a_x=%0.4f";

str(1,1)=sprintf(formatSpec1,cell2mat(ADE3_limits.KGE(2,1)));
str(2,1)=sprintf(formatSpec2,cell2mat(ADE3_limits.KGE(2,2)));
str(3,1)=sprintf(formatSpec3,cell2mat(ADE3_limits.KGE(2,3)));
str(4,1)=sprintf(formatSpec4,cell2mat(ADE3_limits.KGE(2,4)));
str(5,1)=sprintf(formatSpec5,cell2mat(ADE3_limits.KGE(2,5)));
str(6,1)=sprintf(formatSpec6,cell2mat(ADE3_limits.KGE(2,6)));
str(7,1)=sprintf(formatSpec7,cell2mat(ADE3_limits.KGE(2,7)));
str(8,1)=sprintf(formatSpec8,cell2mat(ADE3_limits.KGE(2,8)));

subplot (1,3,3)
plot(Interp_curve_temp3(:,1), Interp_curve_temp3(:,2), 'k', 'LineWidth', 1,'HandleVisibility','off');
hold on;
plot(Interp_curve_temp3(:,1), Interp_curve_temp3(:,3), 'k', 'LineWidth', 1,'HandleVisibility','off');
x2 = [Interp_curve_temp3(:,1)', fliplr(Interp_curve_temp3(:,1)')];
inBetween = [Interp_curve_temp3(:,2)', fliplr(Interp_curve_temp3(:,3)')];
fill(x2, inBetween, [0.86,0.86,0.86]);
hold on
plot(time(:,1),BTC_input(:,2),'-r','LineWidth',2)
legend('Top 0.1% ADE results','Observed BTC')
xlabel ('time [s]');
ylabel ('Cl Conc [mg/l]');
annotation('textbox', [0.79, 0.73, 0.1, 0.1], 'String',str,'FontSize',11,'LineStyle','-')
% annotation('textbox', [0.75, 0.9, 0.1, 0.1], 'String', {"Top 0.1% ADE results from Latin","Hypercube sampling (100'000)"},...
%     'FontSize',12,'LineStyle','none','FitBoxToText','on','HorizontalAlignment','center')
title({"Top 0.1% ADE for KGE results","from Latin Hypercube sampling"},'FontSize',12,'LineStyle','none')
set(gcf, 'WindowState', 'maximized');

clear formatSpec1 formatSpec2 formatSpec3 formatSpec4 formatSpec5 formatSpec6
clear formatSpec7 formatSpec8 formatSpec9 formatSpec10 str

% let's collect in the ADE3_limits all the other obj function as well, for
% completeness

Hyperspace_temp4=sortrows(Hyperspace,6);   % Hyperspace sorted depending on nRMSE values
Hyperspace_temp5=sortrows(Hyperspace,8,'descend');   % Hyperspace sorted depending on log r^2 values
Hyperspace_temp6=sortrows(Hyperspace,9,'descend');   % Hyperspace sorted depending on Pearson values
Hyperspace_temp7=sortrows(Hyperspace,10,'descend');  % Hyperspace sorted depending on logPearson values

Hyperspace_temp4=Hyperspace_temp4(1:top01,:);
Hyperspace_temp5=Hyperspace_temp5(1:top01,:);
Hyperspace_temp6=Hyperspace_temp6(1:top01,:);
Hyperspace_temp7=Hyperspace_temp7(1:top01,:);

%%%% nRMSE

ADE3_limits.nRMSE(1,1)={'v min'};          
ADE3_limits.nRMSE(2,1)=num2cell(min(Hyperspace_temp4(:,1)));
ADE3_limits.nRMSE(1,2)={'v max'};          
ADE3_limits.nRMSE(2,2)=num2cell(max(Hyperspace_temp4(:,1)));
ADE3_limits.nRMSE(1,3)={'A min'};          
ADE3_limits.nRMSE(2,3)=num2cell(min(Hyperspace_temp4(:,2)));
ADE3_limits.nRMSE(1,4)={'A max'};          
ADE3_limits.nRMSE(2,4)=num2cell(max(Hyperspace_temp4(:,2)));
ADE3_limits.nRMSE(1,5)={'D min'};          
ADE3_limits.nRMSE(2,5)=num2cell(min(Hyperspace_temp4(:,3)));
ADE3_limits.nRMSE(1,6)={'D max'};          
ADE3_limits.nRMSE(2,6)=num2cell(max(Hyperspace_temp4(:,3)));
ADE3_limits.nRMSE(1,7)={'nRMSE min'};          
ADE3_limits.nRMSE(2,7)=num2cell(min(Hyperspace_temp4(:,6)));
ADE3_limits.nRMSE(1,8)={'nRMSE max'};          
ADE3_limits.nRMSE(2,8)=num2cell(max(Hyperspace_temp4(:,6)));

%%%% logr^2

ADE3_limits.logr2(1,1)={'v min'};          
ADE3_limits.logr2(2,1)=num2cell(min(Hyperspace_temp5(:,1)));
ADE3_limits.logr2(1,2)={'v max'};          
ADE3_limits.logr2(2,2)=num2cell(max(Hyperspace_temp5(:,1)));
ADE3_limits.logr2(1,3)={'A min'};          
ADE3_limits.logr2(2,3)=num2cell(min(Hyperspace_temp5(:,2)));
ADE3_limits.logr2(1,4)={'A max'};          
ADE3_limits.logr2(2,4)=num2cell(max(Hyperspace_temp5(:,2)));
ADE3_limits.logr2(1,5)={'D min'};          
ADE3_limits.logr2(2,5)=num2cell(min(Hyperspace_temp5(:,3)));
ADE3_limits.logr2(1,6)={'D max'};          
ADE3_limits.logr2(2,6)=num2cell(max(Hyperspace_temp5(:,3)));
ADE3_limits.logr2(1,7)={'nRMSE min'};          
ADE3_limits.logr2(2,7)=num2cell(min(Hyperspace_temp5(:,8)));
ADE3_limits.logr2(1,8)={'nRMSE max'};          
ADE3_limits.logr2(2,8)=num2cell(max(Hyperspace_temp5(:,8)));

%%%% Pearson r2

ADE3_limits.Pearson_r2(1,1)={'v min'};          
ADE3_limits.Pearson_r2(2,1)=num2cell(min(Hyperspace_temp6(:,1)));
ADE3_limits.Pearson_r2(1,2)={'v max'};          
ADE3_limits.Pearson_r2(2,2)=num2cell(max(Hyperspace_temp6(:,1)));
ADE3_limits.Pearson_r2(1,3)={'A min'};          
ADE3_limits.Pearson_r2(2,3)=num2cell(min(Hyperspace_temp6(:,2)));
ADE3_limits.Pearson_r2(1,4)={'A max'};          
ADE3_limits.Pearson_r2(2,4)=num2cell(max(Hyperspace_temp6(:,2)));
ADE3_limits.Pearson_r2(1,5)={'D min'};          
ADE3_limits.Pearson_r2(2,5)=num2cell(min(Hyperspace_temp6(:,3)));
ADE3_limits.Pearson_r2(1,6)={'D max'};          
ADE3_limits.Pearson_r2(2,6)=num2cell(max(Hyperspace_temp6(:,3)));
ADE3_limits.Pearson_r2(1,7)={'nRMSE min'};          
ADE3_limits.Pearson_r2(2,7)=num2cell(min(Hyperspace_temp6(:,9)));
ADE3_limits.Pearson_r2(1,8)={'nRMSE max'};          
ADE3_limits.Pearson_r2(2,8)=num2cell(max(Hyperspace_temp6(:,9)));

%%%% log Pearson r2

ADE3_limits.Pearson_logr2(1,1)={'v min'};          
ADE3_limits.Pearson_logr2(2,1)=num2cell(min(Hyperspace_temp7(:,1)));
ADE3_limits.Pearson_logr2(1,2)={'v max'};          
ADE3_limits.Pearson_logr2(2,2)=num2cell(max(Hyperspace_temp7(:,1)));
ADE3_limits.Pearson_logr2(1,3)={'A min'};          
ADE3_limits.Pearson_logr2(2,3)=num2cell(min(Hyperspace_temp7(:,2)));
ADE3_limits.Pearson_logr2(1,4)={'A max'};          
ADE3_limits.Pearson_logr2(2,4)=num2cell(max(Hyperspace_temp7(:,2)));
ADE3_limits.Pearson_logr2(1,5)={'D min'};          
ADE3_limits.Pearson_logr2(2,5)=num2cell(min(Hyperspace_temp7(:,3)));
ADE3_limits.Pearson_logr2(1,6)={'D max'};          
ADE3_limits.Pearson_logr2(2,6)=num2cell(max(Hyperspace_temp7(:,3)));
ADE3_limits.Pearson_logr2(1,7)={'nRMSE min'};          
ADE3_limits.Pearson_logr2(2,7)=num2cell(min(Hyperspace_temp7(:,10)));
ADE3_limits.Pearson_logr2(1,8)={'nRMSE max'};          
ADE3_limits.Pearson_logr2(2,8)=num2cell(max(Hyperspace_temp7(:,10)));

end

