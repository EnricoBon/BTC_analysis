function [stat1, stat2,Range_param_RMSE,Likelihood_RMSE] = StatFigures(Hyperspace,Description1,Description2)
clear Hyperspace_best a b MAP parameter temp of dat cls y i tmx tmy tm P PP L numdat width

% Let's compute the statistics for only the best 10% of the simulations
% based on RMSE values
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
Hyperspace_best=sortrows(Hyperspace,4);
% Take the 10%
rem=length(Hyperspace_best)*0.1;

% All the stats are computed on the top 10% pf the results
Hyperspace_best=Hyperspace_best(1:rem,:);  
clear rem
%%%
[a,b] = sort(Hyperspace_best(:,4));

% vector "a" puts the RMSE in crescent order while vector "b" returns  
% the location of every RMSE value in the "a" vector

n = length(a);
nbins = 25;      % number of intervals to create the frequency plots (25 -> following Ward et al., 2017)

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

stat1=figure;

subplot(3,3,1)
plot(Hyperspace_best(:,1),Hyperspace_best(:,4),'.k')
hold on
[P,PP] = min(Hyperspace_best(:,4));
plot(Hyperspace_best(PP,1),Hyperspace_best(PP,4),'ro')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Hyperspace_best(:,1)) max(Hyperspace_best(:,1))])
ylim([0 max(Hyperspace_best(:,4))])

subplot(3,3,2)
plot(Hyperspace_best(:,2),Hyperspace_best(:,4),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
hold on
[P,PP] = min(Hyperspace_best(:,4));
plot(Hyperspace_best(PP,2),Hyperspace_best(PP,4),'ro')
xlim([min(Hyperspace_best(:,2)) max(Hyperspace_best(:,2))])
ylim([0 max(Hyperspace_best(:,4))])

subplot(3,3,3)
plot(Hyperspace_best(:,3),Hyperspace_best(:,4),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Hyperspace_best(:,3)) max(Hyperspace_best(:,3))])
ylim([0 max(Hyperspace_best(:,4))])
hold on
[P,PP] = min(Hyperspace_best(:,4));
plot(Hyperspace_best(PP,3),Hyperspace_best(PP,4),'ro')

% Second row - frequency curves 
% We plot here PDFs. So we plot the the distribution of the parameter 
% values CORRESPONDING to a certin threshold range. If the pdf is narrow around 
% a certain range/value -> the parameter is certain and the peak is around the best value
clear parameter temp p1 c1 d1 

subplot (3,3,4)
parameter=Hyperspace_best(:,1);     % Select velocity 
temp = linspace(min(parameter),max(parameter),nbins); % Divide the parameter vector in "nbins" intervals
for i = 1:length(m)
        p1 = parameter(b(1:round((n*m(i)))),1);
        [c1,d1] = hist(p1,temp);
        plot(d1,100*c1./(sum(c1)),'Color',MAP(i*41,:),'LineWidth',2);hold on;
end

% NEW save parameter list to perform a distribution later
Range_param_RMSE.v.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_RMSE.v.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_RMSE.v.five=parameter(b(1:round((n*m(3)))),1);
Range_param_RMSE.v.one=parameter(b(1:round((n*m(4)))),1);
Range_param_RMSE.v.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_RMSE.v.zeroone=parameter(b(1:round((n*m(6)))),1);
% END NEW

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

% NEW save parameter list to perform a distribution later
Range_param_RMSE.A.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_RMSE.A.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_RMSE.A.five=parameter(b(1:round((n*m(3)))),1);
Range_param_RMSE.A.one=parameter(b(1:round((n*m(4)))),1);
Range_param_RMSE.A.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_RMSE.A.zeroone=parameter(b(1:round((n*m(6)))),1);
% END NEW

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

% NEW save parameter list to perform a distribution later
Range_param_RMSE.D.twenty=parameter(b(1:round((n*m(1)))),1);
Range_param_RMSE.D.ten=parameter(b(1:round((n*m(2)))),1);
Range_param_RMSE.D.five=parameter(b(1:round((n*m(3)))),1);
Range_param_RMSE.D.one=parameter(b(1:round((n*m(4)))),1);
Range_param_RMSE.D.zerofive=parameter(b(1:round((n*m(5)))),1);
Range_param_RMSE.D.zeroone=parameter(b(1:round((n*m(6)))),1);
% END NEW

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
% each line represents 1% of all model simulations.
% % % % % % % % 

clear h of dat tmx tmy j

of=Hyperspace_best(:,4);       % Select objective function RMSE
                               % of is automatically sorted from the lower
                               % RMSE to the higher RMSE
of=of./max(of);                % normalise of
of=1-of;                       % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods
[y,i]=sort(of);
dat=Hyperspace_best(i,:);   % take Hyperspace_best in the likely order (from the smaller to the higher)
cls=floor(length(dat)/10);
tmx=zeros(cls,10);
tmy=tmx;

subplot (3,3,7)
set(gcf,'DefaultAxesColorOrder',jet(10));
for j=1:10
% tm takes the dat in the likelyhood order (from the lowest CLS-interval, to the 
% highest CLS-interval). Then we do the cumulative distribution of every
% likelyhood interval to see where the parameters are.
    tm=dat(cls*(j-1)+1:cls*j,1);     % Take first column -> velocity
    tm=sort(tm);
    tmx(:,j)=tm;
    tmy=(1:length(tmx))/cls; % it's a cumulative distribution, every tmx is a step to go from 0 to 1
end
  plot(tmx,tmy,'linewidth',1);hold on;
  plot(tmx(:,10),tmy,'r','linewidth',3);hold on;
  plot(tmx(:,1),tmy,'b','linewidth',3);hold off;
  axis([min(min(tmx)) max(max(tmx)) min(min(tmy)) max(max(tmy))]);
  xlabel(['v [m/s]'])
  ylabel(['cum. norm. RMSE'])
  
% NEW Save likelyhood and values associated
offf=sort(of);
j=1;
Likelihood_RMSE.v.ten(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.v.ten(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=2;
Likelihood_RMSE.v.nine(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.v.nine(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=3;
Likelihood_RMSE.v.eight(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.v.eight(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=4;
Likelihood_RMSE.v.seven(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.v.seven(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=5;
Likelihood_RMSE.v.six(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.v.six(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=6;
Likelihood_RMSE.v.five(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.v.five(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=7;
Likelihood_RMSE.v.four(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.v.four(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=8;
Likelihood_RMSE.v.three(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.v.three(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=9;
Likelihood_RMSE.v.two(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.v.two(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=10;
Likelihood_RMSE.v.one(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.v.one(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

% END NEW
  
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
  ylabel(['cum. norm. RMSE'])
  
% NEW Save likelyhood and values associated
offf=sort(of);
j=1;
Likelihood_RMSE.A.ten(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.A.ten(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=2;
Likelihood_RMSE.A.nine(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.A.nine(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=3;
Likelihood_RMSE.A.eight(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.A.eight(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=4;
Likelihood_RMSE.A.seven(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.A.seven(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=5;
Likelihood_RMSE.A.six(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.A.six(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=6;
Likelihood_RMSE.A.five(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.A.five(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=7;
Likelihood_RMSE.A.four(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.A.four(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=8;
Likelihood_RMSE.A.three(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.A.three(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=9;
Likelihood_RMSE.A.two(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.A.two(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=10;
Likelihood_RMSE.A.one(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.A.one(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

% END NEW  
  
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
  ylabel(['cum. norm. RMSE'])

% NEW Save likelyhood and values associated
offf=sort(of);
j=1;
Likelihood_RMSE.D.ten(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.D.ten(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=2;
Likelihood_RMSE.D.nine(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.D.nine(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=3;
Likelihood_RMSE.D.eight(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.D.eight(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=4;
Likelihood_RMSE.D.seven(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.D.seven(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=5;
Likelihood_RMSE.D.six(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.D.six(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=6;
Likelihood_RMSE.D.five(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.D.five(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=7;
Likelihood_RMSE.D.four(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.D.four(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=8;
Likelihood_RMSE.D.three(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.D.three(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=9;
Likelihood_RMSE.D.two(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.D.two(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

j=10;
Likelihood_RMSE.D.one(:,1)=dat(cls*(j-1)+1:cls*j,1);    % velocity values
Likelihood_RMSE.D.one(:,2)=offf(cls*(j-1)+1:cls*j,1);   % likelihood

% END NEW  
  
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'top 10% latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({Description1;Description2},'FontSize',14);
  
clear h clear of dat tmx tmy j of cls 
set(gcf, 'WindowState', 'maximized');

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

subplot(3,3,1)
plot(Hyperspace_best(:,1),Hyperspace_best(:,4),'.k')
hold on
[P,PP] = min(Hyperspace_best(:,4));
plot(Hyperspace_best(PP,1),Hyperspace_best(PP,4),'ro')
xlabel ('v [m/s]');
ylabel ('RMSE');
xlim([min(Hyperspace_best(:,1)) max(Hyperspace_best(:,1))])
ylim([0 max(Hyperspace_best(:,4))])

subplot(3,3,2)
plot(Hyperspace_best(:,2),Hyperspace_best(:,4),'.k')
xlabel ('A [m^2]');
ylabel ('RMSE');
hold on
[P,PP] = min(Hyperspace_best(:,4));
plot(Hyperspace_best(PP,2),Hyperspace_best(PP,4),'ro')
xlim([min(Hyperspace_best(:,2)) max(Hyperspace_best(:,2))])
ylim([0 max(Hyperspace_best(:,4))])

subplot(3,3,3)
plot(Hyperspace_best(:,3),Hyperspace_best(:,4),'.k')
xlabel ('D [m^2/s]');
ylabel ('RMSE');
xlim([min(Hyperspace_best(:,3)) max(Hyperspace_best(:,3))])
ylim([0 max(Hyperspace_best(:,4))])
hold on
[P,PP] = min(Hyperspace_best(:,4));
plot(Hyperspace_best(PP,3),Hyperspace_best(PP,4),'ro')

% % % % % % % % 
% ROW 2 -> A POSTERIORI PARAMETER DISTRIBUTION
% % % % % % % % 

ncontainers=20;     % 20 containers for the boxes (from OTIS-MCAT, Ward et al., 2017)
slider_value=100;   % We take the 100% of the values because we already 
% have seleceted Hyperspace_best as top 10%
% In case we have no initial cut in the Hyperspace matrix, we should have
% set slider_value=10, to select top 10% of the results

% calculate likelihood
of=Hyperspace_best(:,4);        % select the criterion -> already sorted from the lower to the higher RMSE
                                % RMSE -> low values indicate better models
% NEW
of=of./max(of);
% END                                
LL=1-of; 						% likelihood (high values indicate more likely [probable] models)
if min(LL)<0|min(LL)==0, LL=LL-min(LL)+1000*eps;end; % transform negative lhoods

LL=LL./sum(LL);     % sum(likelihoods)=1      problem if NaN in vector

%Eliminate data that is below the slider threshold
LL=sortrows(LL);
dat=sortrows(Hyperspace_best,4);     % sort Hyperspace based on the RMSE values
LL = flipud(LL);
numdat = floor((slider_value / 100) * size(LL));
LL(numdat+1:size(LL),:) = [];
dat(numdat+1:size(dat),:) = [];

subplot(3,3,4) % velocity

width=(max(dat(:,1))-min(dat(:,1)))/ncontainers; % container width

for n=1:ncontainers
      [K,J]=find(dat(:,1)>min(dat(:,1))+(n-1)*width & dat(:,1)<min(dat(:,1))+n*width); % find all values within the container
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
ylabel(['D (RMSE)'])
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
ylabel(['D (RMSE)'])
xlabel('D [m^2/s]');

% % % % % % % % 
% 3rd ROW - IDENTIFIABILITY PLOTS
% % % % % % % % 
clear of I J dat numdat cls tmx i
dat=Hyperspace_best;%(i,:);   % dat values sorted from the lower RMSE to the higher RMSE
containers=10; % division of each parameter range
grouping=10;   % number of groups

% calculate likelihood
of=Hyperspace_best(:,4);  % criteria (low values indicate better models)
of=1-of; % likelihood (high values indicate more likely [probable] models)
if min(of)<0|min(of)==0, of=of-min(of)+1000*eps;end; % transform negative lhoods

% sort data according to selected perf
[I,J]=sort(of); 
dat=dat(J,:);
% Eliminate data that is below the slider threshold
numdat = floor((slider_value / 100) * size(dat));
dat(numdat+1:size(dat),:) = [];

cls=floor(length(dat)/grouping);
tmx=zeros(cls,grouping); % preallocating tmx and tmy
tmy=tmx;

subplot(3,3,7) % velocity
i=1; % Select velocity
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
   
   [YI]=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping);...
       max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);

% What are we doing with this? We are deducing for every XI interval, where XI
% divides the range of the parameters (from min to max) with the chosen grouping interval,
% the value of YI, which is the number of step from 0 to 1 for our
% comulative curve. You will understand better once you plot this:
% figure
% YI=interp1([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping);...
%        max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],XI);
% plot([min(dat(:,i))-.0001; tmx(1:length(tmx(:,grouping))-1,grouping);...
%        max(tmx(:,grouping))+0.000001; max(dat(:,i))+.0001],[0 tmy 1],'o')
% hold on
% plot(XI,YI,':.') 
% Steeper the line and densier will be the concentration of parameter
% values that give us the top10% of results
   
   [FX] = gradient(YI); % calculate gradient within containers
% Note that sum(FX)~1. So each bin, summed with the other has to give us 
% the total distriution (cdf=1). However, for aesthetic purposes we will
% increase the bins (multiply them by 2 in the "hpatches") to better represent 
% which interval is linked with the higher "slope". So with the higher 
% increase of the cdf. Higher is our bar plot, and more "concentrated" are the
% parameter in that interval that give us better performances

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
   ylabel(['Cum. Dist. RMSE'])
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
   ylabel(['Cum. Dist. RMSE'])
   xlabel('D [m^2/s]');
   
annotation('textbox', [0.12, 0.87, 0.1, 0.1], 'String',{'top 10% latin hypercube results'},'FontSize',12,'LineStyle','none');
sgtitle({Description1;Description2},'FontSize',14);   

set(gcf, 'WindowState', 'maximized');

end





