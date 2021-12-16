function [stat1_r2, stat2_r2, stat1_logRMSE, stat2_logRMSE,stat1_KGE,stat2_KGE,...
    Range_param_r2,Range_param_logRMSE,Range_param_KGE,Range_param_nRMSE]...
    = StatFigures_var(Hyperspace,Description1,Description2)
clear Hyperspace_best a b MAP parameter temp of dat cls y i tmx tmy tm P PP L numdat width
% Let's compute the statistics for only the best 10% of the simulations
% based on r^2 values
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %                           ___  
% %                          |__ \ 
% %                        _ __ ) |
% %                       | '__/ / 
% %                       | | / /_ 
% %                       |_||____|
% %               
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
Hyperspace_best=sortrows(Hyperspace,5,'descend'); % let's act on r^2
% Take the 10%
rem=length(Hyperspace_best)*0.1;
Hyperspace_best=Hyperspace_best(1:rem,:);
clear rem
%%%
[a,b] = sort(Hyperspace_best(:,5),'descend');

% vector "a" puts the RMSE in crescent order while vector "b" returns a 
% the location of every RMSE value in the "a" vector

n = length(a);
nbins = 25;      %number of intervals to create the frequency plots

% thresholds corresponding to the top 20%, 10%, 5%, 1%, 0.5% and 0.01% 
% of the r^2 (NSE) values --> if you want, you can change them
m = [0.2 0.1 0.05 0.01 0.005 0.001];   

% y-axis ranges
minhist = 0;
maxhist = 80;     
MAP = jet(256);

%%%%%%%%%%%%%%%%%
% Figure 1 = v, A and D data vs r2; Frequency plots and normalized
% cumulative distribution

stat1_r2=figure;

subplot(3,3,1)
plot(Hyperspace_best(:,1),Hyperspace_best(:,5),'.k')
hold on
[P,PP] = max(Hyperspace_best(:,5));
plot(Hyperspace_best(PP,1),Hyperspace_best(PP,5),'ro')
xlabel ('v [m/s]');
ylabel ('r^2 - NSE');
xlim([min(Hyperspace_best(:,1)) max(Hyperspace_best(:,1))])
ylim([min(Hyperspace_best(:,5)) 1])

subplot(3,3,2)
plot(Hyperspace_best(:,2),Hyperspace_best(:,5),'.k')
xlabel ('A [m^2]');
ylabel ('r^2 - NSE');
hold on
[P,PP] = max(Hyperspace_best(:,5));
plot(Hyperspace_best(PP,2),Hyperspace_best(PP,5),'ro')
xlim([min(Hyperspace_best(:,2)) max(Hyperspace_best(:,2))])
ylim([min(Hyperspace_best(:,5)) 1])

subplot(3,3,3)
plot(Hyperspace_best(:,3),Hyperspace_best(:,5),'.k')
xlabel ('D [m^2/s]');
ylabel ('r^2 - NSE');
xlim([min(Hyperspace_best(:,3)) max(Hyperspace_best(:,3))])
ylim([min(Hyperspace_best(:,5)) 1])
hold on
[P,PP] = max(Hyperspace_best(:,5));
plot(Hyperspace_best(PP,3),Hyperspace_best(PP,5),'ro')

% Second row - frequency curves 
% this is parameter value probability distribution functions that display
% the distribution of the parameter values CORRESPONDING to a certin
% threshold range. If the narrow around a certain range/value -> the
% parameter is certain and the peak is around the best value
clear parameter temp p1 c1 d1 

subplot (3,3,4)
parameter=Hyperspace_best(:,1);     % Select velocity 
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

subplot (3,3,5)
parameter=Hyperspace_best(:,2);     % Select Area 
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

subplot (3,3,6)
parameter=Hyperspace_best(:,3);     % Select Dispersion 
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
% % % % % % % % 
% Third row - normalized cumulative (Cum. norm.) distributions of RMSE 
% lines represent the best 10% of model runs binned into 1% increments; 
% each line represents 1% of all model simulations
% % % % % % % % 
clear h of dat tmx tmy j

of=Hyperspace_best(:,5);       % Select objective function RMSE
of=of./max(of);                % normalise of
of=1-of;                       % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods
[y,i]=sort(of,'descend');
dat=Hyperspace_best(i,:);
cls=floor(length(dat)/10);
tmx=zeros(cls,10);tmy=tmx;

subplot (3,3,7)
set(gcf,'DefaultAxesColorOrder',(jet(10)));
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
% h=axes('position',[.07 .07 .02 .33]);
% h=colorbar(h);
% set(h,'ytick',[2 10]);
% set(h,'yticklabel',['L';'H']);
h=colorbar;
set(h,'position',[.08 .11 .005 .22]);
set(h,'Ticks',[0 1]);
set(h,'TickLabels',{'L', 'H'});
set(h,'FontSize',11)
h.Label.String = 'Likelihood RMSE';

clear h 

subplot (3,3,8)
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
  ylabel(['cum. norm. r^2'])
  
subplot (3,3,9)
set(gcf,'DefaultAxesColorOrder',(jet(10)));
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
  ylabel(['cum. norm. r^2'])
  
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'top 10% latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({Description1;Description2},'FontSize',14);
  
clear h clear of dat tmx tmy j of cls 
set(gcf, 'WindowState', 'maximized');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Figure 2 = v, A and D data vs RMSE; A posteriori parameter distribution;
% parameter identifiability
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

stat2_r2=figure;

subplot(3,3,1)
plot(Hyperspace_best(:,1),Hyperspace_best(:,5),'.k')
hold on
[P,PP] = max(Hyperspace_best(:,5));
plot(Hyperspace_best(PP,1),Hyperspace_best(PP,5),'ro')
xlabel ('v [m/s]');
ylabel ('r^2 - NSE');
xlim([min(Hyperspace_best(:,1)) max(Hyperspace_best(:,1))])
ylim([min(Hyperspace_best(:,5)) 1])

subplot(3,3,2)
plot(Hyperspace_best(:,2),Hyperspace_best(:,5),'.k')
xlabel ('A [m^2]');
ylabel ('r^2 - NSE');
hold on
[P,PP] = max(Hyperspace_best(:,5));
plot(Hyperspace_best(PP,2),Hyperspace_best(PP,5),'ro')
xlim([min(Hyperspace_best(:,2)) max(Hyperspace_best(:,2))])
ylim([min(Hyperspace_best(:,5)) 1])

subplot(3,3,3)
plot(Hyperspace_best(:,3),Hyperspace_best(:,5),'.k')
xlabel ('D [m^2/s]');
ylabel ('r^2 - NSE');
xlim([min(Hyperspace_best(:,3)) max(Hyperspace_best(:,3))])
ylim([min(Hyperspace_best(:,5)) 1])
hold on
[P,PP] = max(Hyperspace_best(:,5));
plot(Hyperspace_best(PP,3),Hyperspace_best(PP,5),'ro')

% % % % % % % % 
% ROW 2 -> A POSTERIORI PARAMETER DISTRIBUTION
% % % % % % % % 

ncontainers=20;
slider_value=100;   % We take the 100% of the values because we already have seleceted Hyperspace_best as top 10%
% In case we have no initial cut in the Hyperspace matrix, we should have
% set slider_value=10, to select top 10% of the results

% calculate likelihood
of=Hyperspace_best(:,5);             % select the criterion 
                                     %   2^2 -> high values indicate better models
                                
LL=1-of; 						% likelihood (high values indicate more likely [probable] models)
if min(LL)<0|min(LL)==0, LL=LL-min(LL)+1000*eps;end; % transform negative lhoods

LL=LL./sum(LL);     % sum(likelihoods)=1      problem if NaN in vector

%Eliminate data that is below the slider threshold
LL=sortrows(LL);
dat=sortrows(Hyperspace_best,5,'descend');     % sort Hyperspace based on the r^2 values
LL = flipud(LL);
numdat = floor((slider_value / 100) * size(LL));
LL(numdat+1:size(LL),:) = [];
dat(numdat+1:size(dat),:) = [];

subplot(3,3,4) % velocity

width=(max(dat(:,1))-min(dat(:,1)))/ncontainers; % container width -> select velocity

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
ylabel(['D (r2-NSE)'])
xlabel('v [m/s]');
   
subplot(3,3,5) % Area

width=(max(dat(:,2))-min(dat(:,2)))/ncontainers; % container width -> select area

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
ylabel(['D (r2-NSE)'])
xlabel('A [m^2]');

subplot(3,3,6) % Disp

width=(max(dat(:,3))-min(dat(:,3)))/ncontainers; % container width -> select disp

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
ylabel(['D (r2-NSE)'])
xlabel('D [m^2/s]');

% % % % % % % % 
% Third row - IDENTIFIABILITY PLOTS
% % % % % % % % 
clear of I J dat numdat cls tmx

dat=Hyperspace_best;%(i,:);
containers=10; % division of each parameter range
grouping=10;  % number of groups

% calculate likelihood
of=Hyperspace_best(:,5);  % criteria (high values indicate better models)
of=1-of; % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods

% sort data according to selected perf
[I,J]=sort(of,'descend');
dat=dat(J,:);
% Eliminate data that is below the slider threshold
numdat = floor((slider_value / 100) * size(dat));
dat(numdat+1:size(dat),:) = [];

cls=floor(length(dat)/grouping);
tmx=zeros(cls,grouping);tmy=tmx;

subplot(3,3,7) % velocity
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
   ylabel(['Cum. Dist. r^2'])
   xlabel('v [m/s]');
    
subplot(3,3,8) % Area
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
   ylabel(['Cum. Dist. r^2'])
   xlabel('A [m^2]');
   
subplot(3,3,9) % Disp
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
   ylabel(['Cum. Dist. r^2'])
   xlabel('D [m^2/s]');
   
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'top 10% latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({Description1;Description2},'FontSize',14);   

set(gcf, 'WindowState', 'maximized');


% Let's compute the statistics for only the best 10% of the simulations
% based on logRMSE
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %              _             _____  __  __  _____ ______ 
% %             | |           |  __ \|  \/  |/ ____|  ____|
% %             | | ___   __ _| |__) | \  / | (___ | |__   
% %             | |/ _ \ / _` |  _  /| |\/| |\___ \|  __|  
% %             | | (_) | (_| | | \ \| |  | |____) | |____ 
% %             |_|\___/ \__, |_|  \_\_|  |_|_____/|______|
% %                       __/ |                            
% %                      |___/                             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear Hyperspace_best a b MAP parameter temp of dat cls y i tmx tmy tm P PP L numdat width

% Let's compute the statistics for only the best 10% of the simulations
% based on logRMSE values

Hyperspace_best=sortrows(Hyperspace,7);
% Take the 10%
rem=length(Hyperspace_best)*0.1;
Hyperspace_best=Hyperspace_best(1:rem,:);
clear rem
%%%
[a,b] = sort(Hyperspace_best(:,7));

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
% Figure 1 = v, A and D data vs RMSE; Frequency plots and normalized
% cumulative distribution

stat1_logRMSE=figure;

subplot(3,3,1)
plot(Hyperspace_best(:,1),Hyperspace_best(:,7),'.k')
hold on
[P,PP] = min(Hyperspace_best(:,7));
plot(Hyperspace_best(PP,1),Hyperspace_best(PP,7),'ro')
xlabel ('v [m/s]');
ylabel ('logRMSE');
xlim([min(Hyperspace_best(:,1)) max(Hyperspace_best(:,1))])
ylim([0 max(Hyperspace_best(:,7))])

subplot(3,3,2)
plot(Hyperspace_best(:,2),Hyperspace_best(:,7),'.k')
xlabel ('A [m^2]');
ylabel ('logRMSE');
hold on
[P,PP] = min(Hyperspace_best(:,7));
plot(Hyperspace_best(PP,2),Hyperspace_best(PP,7),'ro')
xlim([min(Hyperspace_best(:,2)) max(Hyperspace_best(:,2))])
ylim([0 max(Hyperspace_best(:,7))])

subplot(3,3,3)
plot(Hyperspace_best(:,3),Hyperspace_best(:,7),'.k')
xlabel ('D [m^2/s]');
ylabel ('logRMSE');
xlim([min(Hyperspace_best(:,3)) max(Hyperspace_best(:,3))])
ylim([0 max(Hyperspace_best(:,7))])
hold on
[P,PP] = min(Hyperspace_best(:,7));
plot(Hyperspace_best(PP,3),Hyperspace_best(PP,7),'ro')

% Second row - frequency curves 
% this is parameter value probability distribution functions that display
% the distribution of the parameter values CORRESPONDING to a certin
% threshold range. If the narrow around a certain range/value -> the
% parameter is certain and the peak is around the best value
clear parameter temp p1 c1 d1 

subplot (3,3,4)
parameter=Hyperspace_best(:,1);     % Select velocity 
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

subplot (3,3,5)
parameter=Hyperspace_best(:,2);     % Select Area 
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

subplot (3,3,6)
parameter=Hyperspace_best(:,3);     % Select Dispersion 
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
% % % % % % % % 
% Third row - normalized cumulative (Cum. norm.) distributions of RMSE 
% lines represent the best 10% of model runs binned into 1% increments; 
% each line represents 1% of all model simulations
% % % % % % % % 
clear h of dat tmx tmy j

of=Hyperspace_best(:,7);       % Select objective function logRMSE
of=of./max(of);                % normalise of
of=1-of;                       % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods
[y,i]=sort(of);
dat=Hyperspace_best(i,:);
cls=floor(length(dat)/10);
tmx=zeros(cls,10);tmy=tmx;

subplot (3,3,7)
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
  ylabel(['cum. norm. logRMSE'])
  
colormap(jet(10));
% h=axes('position',[.07 .07 .02 .33]);
% h=colorbar(h);
% set(h,'ytick',[2 10]);
% set(h,'yticklabel',['L';'H']);
h=colorbar;
set(h,'position',[.08 .11 .005 .22]);
set(h,'Ticks',[0 1]);
set(h,'TickLabels',{'L', 'H'});
set(h,'FontSize',11)
h.Label.String = 'Likelihood RMSE';

clear h 

subplot (3,3,8)
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
  ylabel(['cum. norm. logRMSE'])
  
subplot (3,3,9)
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
  ylabel(['cum. norm. logRMSE'])
  
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'top 10% latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({Description1;Description2},'FontSize',14);
  
clear h clear of dat tmx tmy j of cls 
set(gcf, 'WindowState', 'maximized');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Figure 2 = v, A and D data vs RMSE; A posteriori parameter distribution;
% parameter identifiability
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

stat2_logRMSE=figure;

subplot(3,3,1)
plot(Hyperspace_best(:,1),Hyperspace_best(:,7),'.k')
hold on
[P,PP] = min(Hyperspace_best(:,7));
plot(Hyperspace_best(PP,1),Hyperspace_best(PP,7),'ro')
xlabel ('v [m/s]');
ylabel ('logRMSE');
xlim([min(Hyperspace_best(:,1)) max(Hyperspace_best(:,1))])
ylim([0 max(Hyperspace_best(:,7))])

subplot(3,3,2)
plot(Hyperspace_best(:,2),Hyperspace_best(:,7),'.k')
xlabel ('A [m^2]');
ylabel ('logRMSE');
hold on
[P,PP] = min(Hyperspace_best(:,7));
plot(Hyperspace_best(PP,2),Hyperspace_best(PP,7),'ro')
xlim([min(Hyperspace_best(:,2)) max(Hyperspace_best(:,2))])
ylim([0 max(Hyperspace_best(:,7))])

subplot(3,3,3)
plot(Hyperspace_best(:,3),Hyperspace_best(:,7),'.k')
xlabel ('D [m^2/s]');
ylabel ('logRMSE');
xlim([min(Hyperspace_best(:,3)) max(Hyperspace_best(:,3))])
ylim([0 max(Hyperspace_best(:,7))])
hold on
[P,PP] = min(Hyperspace_best(:,7));
plot(Hyperspace_best(PP,3),Hyperspace_best(PP,7),'ro')

% % % % % % % % 
% ROW 2 -> A POSTERIORI PARAMETER DISTRIBUTION
% % % % % % % % 

ncontainers=20;
slider_value=100;   % We take the 100% of the values because we already have seleceted Hyperspace_best as top 10%
% In case we have no initial cut in the Hyperspace matrix, we should have
% set slider_value=10, to select top 10% of the results

% calculate likelihood
of=Hyperspace_best(:,7);             % select the criterion 
                                % logRMSE -> low values indicate better models
                                
LL=1-of; 						% likelihood (high values indicate more likely [probable] models)
if min(LL)<0|min(LL)==0, LL=LL-min(LL)+1000*eps;end; % transform negative lhoods

LL=LL./sum(LL);     % sum(likelihoods)=1      problem if NaN in vector

%Eliminate data that is below the slider threshold
LL=sortrows(LL);
dat=sortrows(Hyperspace_best,7);     % sort Hyperspace based on the RMSE values
LL = flipud(LL);
numdat = floor((slider_value / 100) * size(LL));
LL(numdat+1:size(LL),:) = [];
dat(numdat+1:size(dat),:) = [];

subplot(3,3,4) % velocity

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
ylabel(['D (logRMSE)'])
xlabel('v [m/s]');
   
subplot(3,3,5) % Area

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
ylabel(['D (logRMSE)'])
xlabel('A [m^2]');

subplot(3,3,6) % Disp

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
ylabel(['D (logRMSE)'])
xlabel('D [m^2/s]');

% % % % % % % % 
% 3rd ROW - IDENTIFIABILITY PLOTS
% % % % % % % % 
clear of I J dat numdat cls tmx
dat=Hyperspace_best; %(i,:); % it is already sorted
containers=10; % division of each parameter range
grouping=10;  % number of groups

% calculate likelihood
of=Hyperspace_best(:,7);  % criteria (low values indicate better models)
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

subplot(3,3,7) % velocity
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
   ylabel(['Cum. Dist. logRMSE'])
   xlabel('v [m/s]');
    
subplot(3,3,8) % Area
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
   ylabel(['Cum. Dist. logRMSE'])
   xlabel('A [m^2]');
   
subplot(3,3,9) % Disp
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
   ylabel(['Cum. Dist. logRMSE'])
   xlabel('D [m^2/s]');
   
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'top 10% latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({Description1;Description2},'FontSize',14);   

set(gcf, 'WindowState', 'maximized');



% Let's compute the statistics for only the best 10% of the simulations
% based on KGE values
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %               _  _______ ______ 
% %              | |/ / ____|  ____|
% %              | ' / |  __| |__   
% %              |  <| | |_ |  __|  
% %              | . \ |__| | |____ 
% %              |_|\_\_____|______|
% %                                           
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear Hyperspace_best a b MAP parameter temp of dat cls y i tmx tmy tm P PP L numdat width

Hyperspace_best=sortrows(Hyperspace,11,'descend'); % let's act on KGE
% Take the 10%
rem=length(Hyperspace_best)*0.1;
Hyperspace_best=Hyperspace_best(1:rem,:);
clear rem
%%%
[a,b] = sort(Hyperspace_best(:,11),'descend');

% vector "a" puts the RMSE in crescent order while vector "b" returns a 
% the location of every RMSE value in the "a" vector

n = length(a);
nbins = 25;      %number of intervals to create the frequency plots

% thresholds corresponding to the top 20%, 10%, 5%, 1%, 0.5% and 0.01% 
% of the KGE values --> if you want, you can change them
m = [0.2 0.1 0.05 0.01 0.005 0.001];   

% y-axis ranges
minhist = 0;
maxhist = 80;     
MAP = jet(256);

%%%%%%%%%%%%%%%%%
% Figure 1 = v, A and D data vs KGE; Frequency plots and normalized
% cumulative distribution

stat1_KGE=figure;

subplot(3,3,1)
plot(Hyperspace_best(:,1),Hyperspace_best(:,11),'.k')
hold on
[P,PP] = max(Hyperspace_best(:,11));
plot(Hyperspace_best(PP,1),Hyperspace_best(PP,11),'ro')
xlabel ('v [m/s]');
ylabel ('KGE');
xlim([min(Hyperspace_best(:,1)) max(Hyperspace_best(:,1))])
ylim([min(Hyperspace_best(:,11)) 1])

subplot(3,3,2)
plot(Hyperspace_best(:,2),Hyperspace_best(:,11),'.k')
xlabel ('A [m^2]');
ylabel ('KGE');
hold on
[P,PP] = max(Hyperspace_best(:,11));
plot(Hyperspace_best(PP,2),Hyperspace_best(PP,11),'ro')
xlim([min(Hyperspace_best(:,2)) max(Hyperspace_best(:,2))])
ylim([min(Hyperspace_best(:,11)) 1])

subplot(3,3,3)
plot(Hyperspace_best(:,3),Hyperspace_best(:,11),'.k')
xlabel ('D [m^2/s]');
ylabel ('KGE');
xlim([min(Hyperspace_best(:,3)) max(Hyperspace_best(:,3))])
ylim([min(Hyperspace_best(:,11)) 1])
hold on
[P,PP] = max(Hyperspace_best(:,11));
plot(Hyperspace_best(PP,3),Hyperspace_best(PP,11),'ro')

% Second row - frequency curves 
% this is parameter value probability distribution functions that display
% the distribution of the parameter values CORRESPONDING to a certin
% threshold range. If the narrow around a certain range/value -> the
% parameter is certain and the peak is around the best value
clear parameter temp p1 c1 d1 

subplot (3,3,4)
parameter=Hyperspace_best(:,1);     % Select velocity 
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

subplot (3,3,5)
parameter=Hyperspace_best(:,2);     % Select Area 
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

subplot (3,3,6)
parameter=Hyperspace_best(:,3);     % Select Dispersion 
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
% % % % % % % % 
% Third row - normalized cumulative (Cum. norm.) distributions of RMSE 
% lines represent the best 10% of model runs binned into 1% increments; 
% each line represents 1% of all model simulations
% % % % % % % % 
clear h of dat tmx tmy j

of=Hyperspace_best(:,11);      % Select objective function KGE
of=of./max(of);                % normalise of
of=1-of;                       % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods
[y,i]=sort(of,'descend');
dat=Hyperspace_best(i,:);
cls=floor(length(dat)/10);
tmx=zeros(cls,10);tmy=tmx;

subplot (3,3,7)
set(gcf,'DefaultAxesColorOrder',(jet(10)));
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
  ylabel(['cum. norm. KGE'])
  
colormap(jet(10));
% h=axes('position',[.07 .07 .02 .33]);
% h=colorbar(h);
% set(h,'ytick',[2 10]);
% set(h,'yticklabel',['L';'H']);
h=colorbar;
set(h,'position',[.08 .11 .005 .22]);
set(h,'Ticks',[0 1]);
set(h,'TickLabels',{'L', 'H'});
set(h,'FontSize',11)
h.Label.String = 'Likelihood RMSE';

clear h 

subplot (3,3,8)
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
  ylabel(['cum. norm. KGE'])
  
subplot (3,3,9)
set(gcf,'DefaultAxesColorOrder',(jet(10)));
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
  ylabel(['cum. norm. KGE'])
  
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'top 10% latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({Description1;Description2},'FontSize',14);
  
clear h clear of dat tmx tmy j of cls 
set(gcf, 'WindowState', 'maximized');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Figure 2 = v, A and D data vs RMSE; A posteriori parameter distribution;
% parameter identifiability
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

stat2_KGE=figure;

subplot(3,3,1)
plot(Hyperspace_best(:,1),Hyperspace_best(:,11),'.k')
hold on
[P,PP] = max(Hyperspace_best(:,11));
plot(Hyperspace_best(PP,1),Hyperspace_best(PP,11),'ro')
xlabel ('v [m/s]');
ylabel ('KGE');
xlim([min(Hyperspace_best(:,1)) max(Hyperspace_best(:,1))])
ylim([min(Hyperspace_best(:,11)) 1])

subplot(3,3,2)
plot(Hyperspace_best(:,2),Hyperspace_best(:,11),'.k')
xlabel ('A [m^2]');
ylabel ('KGE');
hold on
[P,PP] = max(Hyperspace_best(:,11));
plot(Hyperspace_best(PP,2),Hyperspace_best(PP,11),'ro')
xlim([min(Hyperspace_best(:,2)) max(Hyperspace_best(:,2))])
ylim([min(Hyperspace_best(:,11)) 1])

subplot(3,3,3)
plot(Hyperspace_best(:,3),Hyperspace_best(:,11),'.k')
xlabel ('D [m^2/s]');
ylabel ('KGE');
xlim([min(Hyperspace_best(:,3)) max(Hyperspace_best(:,3))])
ylim([min(Hyperspace_best(:,11)) 1])
hold on
[P,PP] = max(Hyperspace_best(:,11));
plot(Hyperspace_best(PP,3),Hyperspace_best(PP,11),'ro')

% % % % % % % % 
% ROW 2 -> A POSTERIORI PARAMETER DISTRIBUTION
% % % % % % % % 

ncontainers=20;
slider_value=100;   % We take the 100% of the values because we already have seleceted Hyperspace_best as top 10%
% In case we have no initial cut in the Hyperspace matrix, we should have
% set slider_value=10, to select top 10% of the results

% calculate likelihood
of=Hyperspace_best(:,11);             % select the criterion 
                                      % KGE -> high values indicate better models
                                
LL=1-of; 						% likelihood (high values indicate more likely [probable] models)
if min(LL)<0|min(LL)==0, LL=LL-min(LL)+1000*eps;end; % transform negative lhoods

LL=LL./sum(LL);     % sum(likelihoods)=1      problem if NaN in vector

%Eliminate data that is below the slider threshold
LL=sortrows(LL);
dat=sortrows(Hyperspace_best,11,'descend');     % sort Hyperspace based on the KGE values
LL = flipud(LL);
numdat = floor((slider_value / 100) * size(LL));
LL(numdat+1:size(LL),:) = [];
dat(numdat+1:size(dat),:) = [];

subplot(3,3,4) % velocity

width=(max(dat(:,1))-min(dat(:,1)))/ncontainers; % container width -> select velocity

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
ylabel(['D (KGE)'])
xlabel('v [m/s]');
   
subplot(3,3,5) % Area

width=(max(dat(:,2))-min(dat(:,2)))/ncontainers; % container width -> select area

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
ylabel(['D (KGE)'])
xlabel('A [m^2]');

subplot(3,3,6) % Disp

width=(max(dat(:,3))-min(dat(:,3)))/ncontainers; % container width -> select disp

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
ylabel(['D (KGE)'])
xlabel('D [m^2/s]');

% % % % % % % % 
% Third row - IDENTIFIABILITY PLOTS
% % % % % % % % 
clear of I J dat numdat cls tmx

dat=Hyperspace_best;
containers=10; % division of each parameter range
grouping=10;  % number of groups

% calculate likelihood
of=Hyperspace_best(:,11);  % criteria (high values indicate better models)
of=1-of; % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods

% sort data according to selected perf
[I,J]=sort(of,'descend');
dat=dat(J,:);
% Eliminate data that is below the slider threshold
numdat = floor((slider_value / 100) * size(dat));
dat(numdat+1:size(dat),:) = [];

cls=floor(length(dat)/grouping);
tmx=zeros(cls,grouping);tmy=tmx;

subplot(3,3,7) % velocity
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
   ylabel(['Cum. Dist. KGE'])
   xlabel('v [m/s]');
    
subplot(3,3,8) % Area
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
   ylabel(['Cum. Dist. KGE'])
   xlabel('A [m^2]');
   
subplot(3,3,9) % Disp
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
   ylabel(['Cum. Dist. KGE'])
   xlabel('D [m^2/s]');
   
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'top 10% latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({Description1;Description2},'FontSize',14);   

set(gcf, 'WindowState', 'maximized');



%%%%%% NEW SAVE ALL THE PARAMETER RANGE AND LIKELIHOOD FOR ALL THE
%%%%%% OBJECTIVE FUNCTIONS

%%%% r2 %%%%
Hyperspace_best=sortrows(Hyperspace,5,'descend'); % let's act on r^2
% Take the 10%
rem=length(Hyperspace_best)*0.1;
Hyperspace_best=Hyperspace_best(1:rem,:);
clear rem
[a,b] = sort(Hyperspace_best(:,5),'descend');
%%%
parameter=Hyperspace_best(:,1);     % Select velocity 
Range_param_r2.v.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_r2.v.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_r2.v.five=parameter(b(1:round((n*m(3)))),1);
Range_param_r2.v.one=parameter(b(1:round((n*m(4)))),1);
Range_param_r2.v.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_r2.v.zeroone=parameter(b(1:round((n*m(6)))),1);

parameter=Hyperspace_best(:,2);     % Select Area 
Range_param_r2.A.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_r2.A.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_r2.A.five=parameter(b(1:round((n*m(3)))),1);
Range_param_r2.A.one=parameter(b(1:round((n*m(4)))),1);
Range_param_r2.A.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_r2.A.zeroone=parameter(b(1:round((n*m(6)))),1);

parameter=Hyperspace_best(:,2);     % Select D 
Range_param_r2.D.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_r2.D.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_r2.D.five=parameter(b(1:round((n*m(3)))),1);
Range_param_r2.D.one=parameter(b(1:round((n*m(4)))),1);
Range_param_r2.D.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_r2.D.zeroone=parameter(b(1:round((n*m(6)))),1);

%%%% logRMSE %%%%
Hyperspace_best=sortrows(Hyperspace,7);
% Take the 10%
rem=length(Hyperspace_best)*0.1;
Hyperspace_best=Hyperspace_best(1:rem,:);
clear rem
[a,b] = sort(Hyperspace_best(:,7));
%%%
parameter=Hyperspace_best(:,1);     % Select velocity 
Range_param_logRMSE.v.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_logRMSE.v.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_logRMSE.v.five=parameter(b(1:round((n*m(3)))),1);
Range_param_logRMSE.v.one=parameter(b(1:round((n*m(4)))),1);
Range_param_logRMSE.v.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_logRMSE.v.zeroone=parameter(b(1:round((n*m(6)))),1);

parameter=Hyperspace_best(:,2);     % Select Area 
Range_param_logRMSE.A.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_logRMSE.A.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_logRMSE.A.five=parameter(b(1:round((n*m(3)))),1);
Range_param_logRMSE.A.one=parameter(b(1:round((n*m(4)))),1);
Range_param_logRMSE.A.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_logRMSE.A.zeroone=parameter(b(1:round((n*m(6)))),1);

parameter=Hyperspace_best(:,2);     % Select D 
Range_param_logRMSE.D.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_logRMSE.D.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_logRMSE.D.five=parameter(b(1:round((n*m(3)))),1);
Range_param_logRMSE.D.one=parameter(b(1:round((n*m(4)))),1);
Range_param_logRMSE.D.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_logRMSE.D.zeroone=parameter(b(1:round((n*m(6)))),1);


%%%%% KGE %%%%
Hyperspace_best=sortrows(Hyperspace,11,'descend'); % let's act on KGE
% Take the 10%
rem=length(Hyperspace_best)*0.1;
Hyperspace_best=Hyperspace_best(1:rem,:);
clear rem
%%%
[a,b] = sort(Hyperspace_best(:,11),'descend');

parameter=Hyperspace_best(:,1);     % Select velocity 
Range_param_KGE.v.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_KGE.v.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_KGE.v.five=parameter(b(1:round((n*m(3)))),1);
Range_param_KGE.v.one=parameter(b(1:round((n*m(4)))),1);
Range_param_KGE.v.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_KGE.v.zeroone=parameter(b(1:round((n*m(6)))),1);

parameter=Hyperspace_best(:,2);     % Select Area 
Range_param_KGE.A.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_KGE.A.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_KGE.A.five=parameter(b(1:round((n*m(3)))),1);
Range_param_KGE.A.one=parameter(b(1:round((n*m(4)))),1);
Range_param_KGE.A.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_KGE.A.zeroone=parameter(b(1:round((n*m(6)))),1);

parameter=Hyperspace_best(:,2);     % Select D 
Range_param_KGE.D.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_KGE.D.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_KGE.D.five=parameter(b(1:round((n*m(3)))),1);
Range_param_KGE.D.one=parameter(b(1:round((n*m(4)))),1);
Range_param_KGE.D.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_KGE.D.zeroone=parameter(b(1:round((n*m(6)))),1);

 
%%%% nRMSE %%%%
Hyperspace_best=sortrows(Hyperspace,6);
% Take the 10%
rem=length(Hyperspace_best)*0.1;
Hyperspace_best=Hyperspace_best(1:rem,:);
clear rem
[a,b] = sort(Hyperspace_best(:,6));
%%%
parameter=Hyperspace_best(:,1);     % Select velocity 
Range_param_nRMSE.v.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_nRMSE.v.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_nRMSE.v.five=parameter(b(1:round((n*m(3)))),1);
Range_param_nRMSE.v.one=parameter(b(1:round((n*m(4)))),1);
Range_param_nRMSE.v.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_nRMSE.v.zeroone=parameter(b(1:round((n*m(6)))),1);

parameter=Hyperspace_best(:,2);     % Select Area 
Range_param_nRMSE.A.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_nRMSE.A.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_nRMSE.A.five=parameter(b(1:round((n*m(3)))),1);
Range_param_nRMSE.A.one=parameter(b(1:round((n*m(4)))),1);
Range_param_nRMSE.A.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_nRMSE.A.zeroone=parameter(b(1:round((n*m(6)))),1);

parameter=Hyperspace_best(:,2);     % Select D 
Range_param_nRMSE.D.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_nRMSE.D.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_nRMSE.D.five=parameter(b(1:round((n*m(3)))),1);
Range_param_nRMSE.D.one=parameter(b(1:round((n*m(4)))),1);
Range_param_nRMSE.D.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_nRMSE.D.zeroone=parameter(b(1:round((n*m(6)))),1);
end
