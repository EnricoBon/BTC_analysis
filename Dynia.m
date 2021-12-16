function [f,Dynia_param_range] = Dynia(TEMP_BTC,Order_Working_matrix,BTC_input,Description1, Description2)
% 
% Code based on MCAT (Thorsten Wagener, Penn State, October 2004)

t=size(TEMP_BTC,1); % length of time series;

ff(1,1) = 5;    % = Number of parameters; v, A, D, Alpha, A_TS
ff(1,2) = 8;    % = Number of Objective functions; RMSE, r2, nRMSE, logRMSE, logr2, Pearson_r2, Pearson_logr2, KGE
ff(1,3) = 1;    % nvars
ff(1,4) = length(Order_Working_matrix(1:100,1));  % Number of simulations (take the top 100 BTC)       

pars(:,1)=Order_Working_matrix(1:100,1);  % Velocity
pars(:,2)=Order_Working_matrix(1:100,2);  % Area
pars(:,3)=Order_Working_matrix(1:100,3);  % Dispersion
pars(:,4)=Order_Working_matrix(1:100,4);  % Alpha
pars(:,5)=Order_Working_matrix(1:100,5);  % A_TS

crit=[Order_Working_matrix(1:100,8:15)];
dat=[pars crit]; 	      % data matrix, in the original formulation dat=[pars crit vars]; but we don't have vars 

mct=TEMP_BTC;%';      % Every i-th line has the BTC simulation for the i-th parameter in the pars matrix. Number of column
                      % is the number of timestep of the BTC
                    
obs=BTC_input(:,2); % observed BTC (just the concentration);

t=size(mct,1);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % %            _            _ _         
% % % % %           | |          (_) |        
% % % % % __   _____| | ___   ___ _| |_ _   _ 
% % % % % \ \ / / _ \ |/ _ \ / __| | __| | | |
% % % % %  \ V /  __/ | (_) | (__| | |_| |_| |
% % % % %   \_/ \___|_|\___/ \___|_|\__|\__, |
% % % % %                                __/ |
% % % % %                               |___/ 


p1      = 1;    % selected parameter
window  = 3;    % * 2 + 1 = size of moving window
% Window size equal to 3 because "the concentration history is sensitive to 
% the parameters over relatively narrow time spans that are associated with 
% specific segments of the concentration history" (Wagner & Harvey 1997).
% And used also by by Wagener et al., 2002.

% 'fixed' algorithm parameters ********************************************
containers  = 20; % split of parameter range
aa          =  1; % [1] medium window (running mean) [2] right boundary window (regressive)
grouping    = 10; % number of 'horizontal' groups

% 'fixed' algorithm parameters in the code
ff(1)       =  p1; 	% last selected parameter equals first one
timestep    = ' ';  % time-series time step
timestep='samples';

% DYNIA *******************************************************************
h_wait=waitbar(0,'Running DYNIA Algorithm - velocity');  
for dt=1:t

   % calculate "likelihood" at every time step dt
   if aa==1 % medium window (running mean)
   	for no=1:size(mct,2)
   	   if dt<window+1
   	   	    residual(no,dt)=mean(abs(mct(1:dt+window,no)-obs(1:dt+window))); % mean absolute error in moving window
   	   elseif t-dt<window+1
   	        residual(no,dt)=mean(abs(mct(dt-window:t,no)-obs(dt-window:t))); % mean absolute error in moving window
   	   else
			residual(no,dt)=mean(abs(mct(dt-window:dt+window,no)-obs(dt-window:dt+window))); % mean absolute error in moving window         
   	   end      
   	end
	elseif aa==2 % right boundary window (recursive)
   	for no=1:size(mct,2)
   	    if dt<2*window+1
   	   	    residual(no,dt)=mean(abs(mct(1:dt+2*window,no)-obs(1:dt+2*window))); % mean absolute error in moving window
        else
			residual(no,dt)=mean(abs(mct(dt-2*window:dt,no)-obs(dt-2*window:dt))); % mean absolute error in moving window         
   	    end      
   	end
	end            
   
   if max(residual)~=0
      L=residual(:,dt)./max(residual(:,dt)); % normalise criterion
   else
      L=residual(:,dt);
   end
   L=1-L; % likelihood (high values indicate more likely [probable] models)
   if min(L)<0|min(L)==0, L=L-min(L)+1000*eps;end; % transform negative lhoods
   L=L';
   
   % sort data according to selected perf measure
	[I,J]=sort(L);
	pop=dat(J,:);
	cls=floor(length(pop)/grouping);
	tmx=zeros(cls,grouping);tmy=tmx;
      
	for i=p1:ff(1)
      
      % calculate cumulative distribution of top 10% of the model population ******************************************************************
      
      tm=pop(cls*(grouping-1)+1:cls*grouping,i);
      tmx=sort(tm);
	  tmy=(1:length(tmx))/cls; 
      
      % calculate the 90% confidence limits *************************************************************************************************
      
      tucl=find(tmy>0.95); ucl(dt)=tmx(tucl(1));
      tlcl=find(tmy>0.05); lcl(dt)=tmx(tlcl(1));
      
	  % calculate gradients *******************************************************************************************************************
	   
	  step=(max(dat(:,i))-min(dat(:,i)))/containers;
	  XI(:,i)=[min(dat(:,i)):step:max(dat(:,i))]';
      
      temp1=min(dat(:,i))-.001;
      if temp1<0
         temp1=0;
      end
      temp2=max(dat(:,i))+.001;    

      for kk=1:length(tmx)-2
          if tmx(kk)==tmx(kk+1)
              tmx(kk+1)=tmx(kk+1)+[tmx(kk+2)-tmx(kk)]/1000;
          end
      end
      
      [YI]=interp_special(sort([temp1;tmx;temp2]),[0 tmy 1],XI(:,i));
      
      [FX(:,dt,i)] = gradient(YI); % calculate gradient within containers
      
      % select best performing model as a function of the highest gradient
      % try select best model as point estimate!!!!
      
      location=find(FX(:,dt,i)==max(FX(:,dt,i)));
      best(dt)=XI(location,i);
      
      % find number of models within best segment (i.e. pixel)
      
      if location<size(XI,1)
          [YY,ZZ]=find(tmx>XI(location,i)&tmx<XI(location+1,i));
      else
          [YY,ZZ]=find(tmx>XI(location,i));
      end
                 
      clusterline_pixel(dt)=size(YY,1);
      clear YY ZZ
      
      % find number of models within 90% cfls
      
      [YY,ZZ]=find(tmx>lcl(dt)&tmx<ucl(dt));
      clusterline_cfls(dt)=size(YY,1);      
      clear YY
      
	  % **********************************************************************************************************************************

      dtcolor(:,dt,i)=FX(:,dt,i); % let color values run from min (bottom) to max (top) parameter value
       
      clear tm YI
       
	end
   
   clear L 
   waitbar(dt/t);
   
end
close(h_wait);
G = find(isnan(dtcolor)); dtcolor(G)=zeros(size(G));

for i=p1:ff(1) % every parameter

	infocontent=1-[(ucl-lcl)/(max(dat(:,i))-min(dat(:,i)))];
   
	nucl=ucl/(max(dat(:,i))-min(dat(:,i))); % normalise ucl
	nlcl=lcl/(max(dat(:,i))-min(dat(:,i))); % normalise lcl

end

% PLOTTING ****************************************************************

% create x matrix for patch function
% one x vector for every time step
% matrix is the same for all parameters
for dt=1:t
   xmatrix(1:2,dt)=dt-0.5;
   xmatrix(3:4,dt)=dt+0.5;
end
% create y matrix for patch function
for C=1:containers
   ymatrix(1,C)=(C-1)*(1/containers);
   ymatrix(2:3,C)=C*(1/containers);
   ymatrix(4,C)=(C-1)*(1/containers);
end

f=figure; %('visible','off'); % Let's help our GPU a lil' bit...
f.WindowState='fullscreen';
for i=p1:ff(1) % every parameter
  
  subplot(5,2,1); % DYNIA PLOT ********************************************
  yyaxis left
  for dt=1:t % every time step
     for C=1:containers % every container
		tcolor=[abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i)))) abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i)))) abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i))))];
        tcolor=[1 tcolor(1,2) tcolor(1,3)]; % let's plot them in RED shades
        patch(xmatrix(:,dt),ymatrix(:,C).*(max(XI(:,i))-min(XI(:,i)))+min(XI(:,i)),tcolor,'edgecolor',tcolor);hold on;
     end
  end
  
  plot(ucl,':','color','k','linewidth',1);hold on; % upper confidence limit
  plot(lcl,':','color','k','linewidth',1);hold on; % lower confidence limit
  axis([1 t min(XI(:,i)) max(XI(:,i))]);
  ylabel(['v [m/s]']);
  yyaxis right
  plot((1:1:t),BTC_input(:,2),'k','linewidth',2);
  hold on;
  axis('auto'); 
  ylabel(['Cl [mg/l]']);
  xlim([0 t]);

title({'Likelihood distribution as function of';'parameter values at each time step'}); 

  set(gca,'layer','top');
  set(get(gca,'xlabel'),'fontsize',16); % Before it was 14
  set(get(gca,'ylabel'),'fontsize',16); % Before it was 14
  set(get(gca,'yaxis'),'Color','black');
  
  set(gca,'layer','top');
  set(gca,'fontsize',16,'linewidth',2);
  set(get(gca,'title'),'fontsize',18); % Before it was 14

  subplot(5,2,2); % INFORMATION CONTENT PLOT *************************************************************************************************************
  
  infocontent=1-[(ucl-lcl)/(max(dat(:,i))-min(dat(:,i)))];
  
  bar(infocontent,'r');hold on;

  ylim([0 1]);
  ylabel({'Inform.';'cont. v []'});
  yyaxis right
  plot((1:1:t),BTC_input(:,2),'color','black','linewidth',2); % 
  hold on;
  axis('auto'); 
  ylabel(['Cl [mg/l]']);
%   save infocontent infocontent;
  xlim([0 t]);
  
title({'Parameter information content';'at each time step'}); 
  
  set(gca,'layer','top');
  set(get(gca,'xlabel'),'fontsize',16); % Before it was 14
  set(get(gca,'ylabel'),'fontsize',16); % Before it was 14
  
  set(gca,'layer','top');
  set(gca,'fontsize',16,'linewidth',2);
  set(get(gca,'yaxis'),'Color','black');
  set(get(gca,'title'),'fontsize',18); % Before it was 14
  
end
% 
Dynia_param_range.v(:,1)=BTC_input(:,1);
Dynia_param_range.v(:,2)=lcl';
Dynia_param_range.v(:,3)=ucl';
Dynia_param_range.v(:,4)=infocontent;
% save Dynia_v.mat

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % %                          
% % % % %     /\                   
% % % % %    /  \   _ __ ___  __ _ 
% % % % %   / /\ \ | '__/ _ \/ _` |
% % % % %  / ____ \| | |  __/ (_| |
% % % % % /_/    \_\_|  \___|\__,_|

p1      = 2;    % selected parameter
window  = 3;    % * 2 + 1 = size of moving window
% Window size equal to 3 because "the concentration history is sensitive to 
% the parameters over relatively narrow time spans that are associated with 
% specific segments of the concentration history" (Wagner & Harvey 1997).
% And used also by by Wagener et al., 2002.

% 'fixed' algorithm parameters ********************************************************************************************************************
containers  = 20; % split of parameter range
aa          =  1; % [1] medium window (running mean) [2] right boundary window (regressive)
grouping    = 10; % number of 'horizontal' groups

% 'fixed' algorithm parameters in the code
ff(1)       =  p1; 	% last selected parameter equals first one
timestep    = ' ';  % time-series time step

   timestep='samples';
% end

% DYNIA *******************************************************************
h_wait=waitbar(0,'Running DYNIA Algorithm - Area');  
for dt=1:t

   % calculate "likelihood" at every time step dt
   if aa==1 % medium window (running mean)
   	for no=1:size(mct,2)
   	   if dt<window+1
   	   	    residual(no,dt)=mean(abs(mct(1:dt+window,no)-obs(1:dt+window))); % mean absolute error in moving window
   	   elseif t-dt<window+1
   	        residual(no,dt)=mean(abs(mct(dt-window:t,no)-obs(dt-window:t))); % mean absolute error in moving window
   	   else
			residual(no,dt)=mean(abs(mct(dt-window:dt+window,no)-obs(dt-window:dt+window))); % mean absolute error in moving window         
   	   end      
   	end
	elseif aa==2 % right boundary window (recursive)
   	for no=1:size(mct,2)
   	    if dt<2*window+1
   	   	    residual(no,dt)=mean(abs(mct(1:dt+2*window,no)-obs(1:dt+2*window))); % mean absolute error in moving window
        else
			residual(no,dt)=mean(abs(mct(dt-2*window:dt,no)-obs(dt-2*window:dt))); % mean absolute error in moving window         
   	    end      
   	end
	end            
   
   if max(residual)~=0
      L=residual(:,dt)./max(residual(:,dt)); % normalise criterion
   else
      L=residual(:,dt);
   end
   L=1-L; % likelihood (high values indicate more likely [probable] models)
   if min(L)<0|min(L)==0, L=L-min(L)+1000*eps;end; % transform negative lhoods
   L=L';
   
   % sort data according to selected perf measure
	[I,J]=sort(L);
	pop=dat(J,:);
	cls=floor(length(pop)/grouping);
	tmx=zeros(cls,grouping);tmy=tmx;
      
	for i=p1:ff(1)
      
      % calculate cumulative distribution of top 10% of the model population ******************************************************************
      
      tm=pop(cls*(grouping-1)+1:cls*grouping,i);
      tmx=sort(tm);
	  tmy=(1:length(tmx))/cls; 
      
      % calculate the 90% confidence limits *************************************************************************************************
      
      tucl=find(tmy>0.95); ucl(dt)=tmx(tucl(1));
      tlcl=find(tmy>0.05); lcl(dt)=tmx(tlcl(1));
      
	  % calculate gradients *******************************************************************************************************************
	   
	  step=(max(dat(:,i))-min(dat(:,i)))/containers;
	  XI(:,i)=[min(dat(:,i)):step:max(dat(:,i))]';
      
      temp1=min(dat(:,i))-.001;
      if temp1<0
         temp1=0;
      end
      temp2=max(dat(:,i))+.001;    

      for kk=1:length(tmx)-2
          if tmx(kk)==tmx(kk+1)
              tmx(kk+1)=tmx(kk+1)+[tmx(kk+2)-tmx(kk)]/1000;
          end
      end
      
      [YI]=interp_special(sort([temp1;tmx;temp2]),[0 tmy 1],XI(:,i));
      
      [FX(:,dt,i)] = gradient(YI); % calculate gradient within containers
      
      % select best performing model as a function of the highest gradient
      % try select best model as point estimate!!!!
      
      location=find(FX(:,dt,i)==max(FX(:,dt,i)));
      best(dt)=XI(location,i);
      
      % find number of models within best segment (i.e. pixel)
      
      if location<size(XI,1)
          [YY,ZZ]=find(tmx>XI(location,i)&tmx<XI(location+1,i));
      else
          [YY,ZZ]=find(tmx>XI(location,i));
      end
                 
      clusterline_pixel(dt)=size(YY,1);
      clear YY ZZ
      
      % find number of models within 90% cfls
      
      [YY,ZZ]=find(tmx>lcl(dt)&tmx<ucl(dt));
      clusterline_cfls(dt)=size(YY,1);      
      clear YY
      
	  % **********************************************************************************************************************************

      dtcolor(:,dt,i)=FX(:,dt,i); % let color values run from min (bottom) to max (top) parameter value
       
      clear tm YI
       
	end
   
   clear L 
   waitbar(dt/t);
   
end
close(h_wait);
G = find(isnan(dtcolor)); dtcolor(G)=zeros(size(G));

for i=p1:ff(1) % every parameter

	infocontent=1-[(ucl-lcl)/(max(dat(:,i))-min(dat(:,i)))];
   
	nucl=ucl/(max(dat(:,i))-min(dat(:,i))); % normalise ucl
	nlcl=lcl/(max(dat(:,i))-min(dat(:,i))); % normalise lcl

end

% PLOTTING ****************************************************************
% 
% create x matrix for patch function
% one x vector for every time step
% matrix is the same for all parameters
for dt=1:t
   xmatrix(1:2,dt)=dt-0.5;
   xmatrix(3:4,dt)=dt+0.5;
end
% create y matrix for patch function
for C=1:containers
   ymatrix(1,C)=(C-1)*(1/containers);
   ymatrix(2:3,C)=C*(1/containers);
   ymatrix(4,C)=(C-1)*(1/containers);
end

% figure
for i=p1:ff(1) % every parameter
  
  subplot(5,2,3); % DYNIA PLOT ****************************************************************************************************************************
   yyaxis left
  for dt=1:t % every time step
     for C=1:containers % every container
		tcolor=[abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i)))) abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i)))) abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i))))];
        tcolor=[1 tcolor(1,2) tcolor(1,3)]; % let's plot them in RED shades
        patch(xmatrix(:,dt),ymatrix(:,C).*(max(XI(:,i))-min(XI(:,i)))+min(XI(:,i)),tcolor,'edgecolor',tcolor);hold on;
     end
  end
  plot(ucl,':','color','k','linewidth',1);hold on; % upper confidence limit
  plot(lcl,':','color','k','linewidth',1);hold on; % lower confidence limit
  axis([1 t min(XI(:,i)) max(XI(:,i))]);
  ylabel(['A [m^2]']);
  
%   plot([obs./(max(obs)+0.1*max(obs))]*[max(XI(:,i))-min(XI(:,i))]+[min(XI(:,i))+0.01*(max(XI(:,i))-min(XI(:,i)))],'k','linewidth',2);hold on; % normalized observed flow
  yyaxis right
%   plot([obs./(max(obs)+0.1*max(obs))]*[max(XI(:,i))-min(XI(:,i))]+[min(XI(:,i))+0.01*(max(XI(:,i))-min(XI(:,i)))],'k','linewidth',2);hold on; % normalized observed flow
  plot((1:1:t),BTC_input(:,2),'k','linewidth',2);
  hold on;
  axis('auto'); 
  ylabel(['Cl [mg/l]']);
  xlim([0 t]);
  
  set(gca,'layer','top');
  set(get(gca,'xlabel'),'fontsize',16);
  set(get(gca,'ylabel'),'fontsize',16);
%   set(get(gca,'title'),'fontsize',14);
  set(get(gca,'yaxis'),'Color','black');

  set(gca,'layer','top');
  set(gca,'fontsize',16,'linewidth',2);
  
  subplot(5,2,4); % INFORMATION CONTENT PLOT *************************************************************************************************************
  
  infocontent=1-[(ucl-lcl)/(max(dat(:,i))-min(dat(:,i)))];
  
  bar(infocontent,'r');hold on;
  ylim([0 1]);
  ylabel({'Inform.';'cont. A []'});
%   plot([obs./(max(obs)+0.1*max(obs))]+0.025,'linewidth',2,'color',[.7 .7 .7]);hold on;
  yyaxis right
%   plot([obs./(max(obs)+0.1*max(obs))]*[max(XI(:,i))-min(XI(:,i))]+[min(XI(:,i))+0.01*(max(XI(:,i))-min(XI(:,i)))],'k','linewidth',2);hold on; % normalized observed flow
  plot((1:1:t),BTC_input(:,2),'color','black','linewidth',2);
  hold on;
  axis('auto'); 
  ylabel(['Cl [mg/l]']);
%   save infocontent infocontent;
  xlim([0 t]);
  save infocontent infocontent;
  
  set(gca,'layer','top');
  set(get(gca,'xlabel'),'fontsize',14);
  set(get(gca,'ylabel'),'fontsize',14);
  set(get(gca,'title'),'fontsize',14);
  set(get(gca,'yaxis'),'Color','black');
  
  set(gca,'layer','top');
  set(gca,'fontsize',16,'linewidth',2);
  
%  sgtitle({'DYNIA algorithm for TSM parameters';Description1;Description2},'FontSize',14);
end
Dynia_param_range.A(:,1)=BTC_input(:,1);
Dynia_param_range.A(:,2)=lcl';
Dynia_param_range.A(:,3)=ucl';
Dynia_param_range.A(:,4)=infocontent;
% save Dynia_A.mat

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % %  _____  _                         _             
% % % % % |  __ \(_)                       (_)            
% % % % % | |  | |_ ___ _ __   ___ _ __ ___ _  ___  _ __  
% % % % % | |  | | / __| '_ \ / _ \ '__/ __| |/ _ \| '_ \ 
% % % % % | |__| | \__ \ |_) |  __/ |  \__ \ | (_) | | | |
% % % % % |_____/|_|___/ .__/ \___|_|  |___/_|\___/|_| |_|
% % % % %              | |                                
% % % % %              |_|                                

p1      = 3;    % selected parameter
window  = 3;    % * 2 + 1 = size of moving window
% Window size equal to 3 because "the concentration history is sensitive to 
% the parameters over relatively narrow time spans that are associated with 
% specific segments of the concentration history" (Wagner & Harvey 1997).
% And used also by by Wagener et al., 2002.

% 'fixed' algorithm parameters ********************************************************************************************************************
containers  = 20; % split of parameter range
aa          =  1; % [1] medium window (running mean) [2] right boundary window (regressive)
grouping    = 10; % number of 'horizontal' groups

% 'fixed' algorithm parameters in the code
ff(1)       =  p1; 	% last selected parameter equals first one
timestep    = ' ';  % time-series time step

timestep='samples';

% DYNIA *******************************************************************
h_wait=waitbar(0,'Running DYNIA Algorithm - Disp');  
for dt=1:t

   % calculate "likelihood" at every time step dt
   if aa==1 % medium window (running mean)
   	for no=1:size(mct,2)
   	   if dt<window+1
   	   	    residual(no,dt)=mean(abs(mct(1:dt+window,no)-obs(1:dt+window))); % mean absolute error in moving window
   	   elseif t-dt<window+1
   	        residual(no,dt)=mean(abs(mct(dt-window:t,no)-obs(dt-window:t))); % mean absolute error in moving window
   	   else
			residual(no,dt)=mean(abs(mct(dt-window:dt+window,no)-obs(dt-window:dt+window))); % mean absolute error in moving window         
   	   end      
   	end
	elseif aa==2 % right boundary window (recursive)
   	for no=1:size(mct,2)
   	    if dt<2*window+1
   	   	    residual(no,dt)=mean(abs(mct(1:dt+2*window,no)-obs(1:dt+2*window))); % mean absolute error in moving window
        else
			residual(no,dt)=mean(abs(mct(dt-2*window:dt,no)-obs(dt-2*window:dt))); % mean absolute error in moving window         
   	    end      
   	end
	end            
   
   if max(residual)~=0
      L=residual(:,dt)./max(residual(:,dt)); % normalise criterion
   else
      L=residual(:,dt);
   end
   L=1-L; % likelihood (high values indicate more likely [probable] models)
   if min(L)<0|min(L)==0, L=L-min(L)+1000*eps;end; % transform negative lhoods
   L=L';
   
   % sort data according to selected perf measure
	[I,J]=sort(L);
	pop=dat(J,:);
	cls=floor(length(pop)/grouping);
	tmx=zeros(cls,grouping);tmy=tmx;
      
	for i=p1:ff(1)
      
      % calculate cumulative distribution of top 10% of the model population ******************************************************************
      
      tm=pop(cls*(grouping-1)+1:cls*grouping,i);
      tmx=sort(tm);
	  tmy=(1:length(tmx))/cls; 
      
      % calculate the 90% confidence limits *************************************************************************************************
      
      tucl=find(tmy>0.95); ucl(dt)=tmx(tucl(1));
      tlcl=find(tmy>0.05); lcl(dt)=tmx(tlcl(1));
      
	  % calculate gradients *******************************************************************************************************************
	   
	  step=(max(dat(:,i))-min(dat(:,i)))/containers;
	  XI(:,i)=[min(dat(:,i)):step:max(dat(:,i))]';
      
      temp1=min(dat(:,i))-.001;
      if temp1<0
         temp1=0;
      end
      temp2=max(dat(:,i))+.001;    

      for kk=1:length(tmx)-2
          if tmx(kk)==tmx(kk+1)
              tmx(kk+1)=tmx(kk+1)+[tmx(kk+2)-tmx(kk)]/1000;
          end
      end
      
      [YI]=interp_special(sort([temp1;tmx;temp2]),[0 tmy 1],XI(:,i));
      
      [FX(:,dt,i)] = gradient(YI); % calculate gradient within containers
      
      % select best performing model as a function of the highest gradient
      % try select best model as point estimate!!!!
      
      location=find(FX(:,dt,i)==max(FX(:,dt,i)));
      best(dt)=XI(location,i);
      
      % find number of models within best segment (i.e. pixel)
      
      if location<size(XI,1)
          [YY,ZZ]=find(tmx>XI(location,i)&tmx<XI(location+1,i));
      else
          [YY,ZZ]=find(tmx>XI(location,i));
      end
                 
      clusterline_pixel(dt)=size(YY,1);
      clear YY ZZ
      
      % find number of models within 90% cfls
      
      [YY,ZZ]=find(tmx>lcl(dt)&tmx<ucl(dt));
      clusterline_cfls(dt)=size(YY,1);      
      clear YY
      
	  % **********************************************************************************************************************************

      dtcolor(:,dt,i)=FX(:,dt,i); % let color values run from min (bottom) to max (top) parameter value
       
      clear tm YI
       
	end
   
   clear L 
   waitbar(dt/t);
   
end
close(h_wait);
G = find(isnan(dtcolor)); dtcolor(G)=zeros(size(G));

for i=p1:ff(1) % every parameter

	infocontent=1-[(ucl-lcl)/(max(dat(:,i))-min(dat(:,i)))];
   
	nucl=ucl/(max(dat(:,i))-min(dat(:,i))); % normalise ucl
	nlcl=lcl/(max(dat(:,i))-min(dat(:,i))); % normalise lcl

end

% PLOTTING ****************************************************************

% create x matrix for patch function
% one x vector for every time step
% matrix is the same for all parameters
for dt=1:t
   xmatrix(1:2,dt)=dt-0.5;
   xmatrix(3:4,dt)=dt+0.5;
end
% create y matrix for patch function
for C=1:containers
   ymatrix(1,C)=(C-1)*(1/containers);
   ymatrix(2:3,C)=C*(1/containers);
   ymatrix(4,C)=(C-1)*(1/containers);
end

% figure
for i=p1:ff(1) % every parameter
  
  subplot(5,2,5); % DYNIA PLOT ****************************************************************************************************************************
  yyaxis left 
  for dt=1:t % every time step
     for C=1:containers % every container
		tcolor=[abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i)))) abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i)))) abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i))))];
        tcolor=[1 tcolor(1,2) tcolor(1,3)]; % let's plot them in RED shades
        patch(xmatrix(:,dt),ymatrix(:,C).*(max(XI(:,i))-min(XI(:,i)))+min(XI(:,i)),tcolor,'edgecolor',tcolor);hold on;
     end
  end
  plot(ucl,':','color','k','linewidth',1);hold on; % upper confidence limit
  plot(lcl,':','color','k','linewidth',1);hold on; % lower confidence limit
  axis([1 t min(XI(:,i)) max(XI(:,i))]);
  ylabel(['D [m^2/s]']);
  
%   plot([obs./(max(obs)+0.1*max(obs))]*[max(XI(:,i))-min(XI(:,i))]+[min(XI(:,i))+0.01*(max(XI(:,i))-min(XI(:,i)))],'k','linewidth',2);hold on; % normalized observed flow
  yyaxis right
%   plot([obs./(max(obs)+0.1*max(obs))]*[max(XI(:,i))-min(XI(:,i))]+[min(XI(:,i))+0.01*(max(XI(:,i))-min(XI(:,i)))],'k','linewidth',2);hold on; % normalized observed flow
  plot((1:1:t),BTC_input(:,2),'k','linewidth',2);
  hold on;
  axis('auto'); 
  ylabel(['Cl [mg/l]']);
  xlim([0 t]);
  
  set(gca,'layer','top');
  set(get(gca,'xlabel'),'fontsize',16);
  set(get(gca,'ylabel'),'fontsize',16);
  set(get(gca,'title'),'fontsize',18);
  set(get(gca,'yaxis'),'Color','black');

  set(gca,'layer','top');
  set(gca,'fontsize',16,'linewidth',2);
  
  subplot(5,2,6); % INFORMATION CONTENT PLOT *************************************************************************************************************
  
  infocontent=1-[(ucl-lcl)/(max(dat(:,i))-min(dat(:,i)))];
  
 bar(infocontent,'r');hold on;
 ylim([0 1]);
  ylabel({'Inform.';'cont. D []'});
  yyaxis right
  plot((1:1:t),BTC_input(:,2),'color','black','linewidth',2);
  hold on;
  axis('auto'); 
  ylabel(['Cl [mg/l]']);
  xlim([0 t]);
  
  save infocontent infocontent;
  
  set(gca,'layer','top');
  set(get(gca,'xlabel'),'fontsize',16);
  set(get(gca,'ylabel'),'fontsize',16);
%   set(get(gca,'title'),'fontsize',14);
  set(get(gca,'yaxis'),'Color','black');

  set(gca,'layer','top');
  set(gca,'fontsize',16,'linewidth',2);
  
end
% save Dynia_D.mat
Dynia_param_range.D(:,1)=BTC_input(:,1);
Dynia_param_range.D(:,2)=lcl';
Dynia_param_range.D(:,3)=ucl';
Dynia_param_range.D(:,4)=infocontent;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % %           _       _           
% % % % %     /\   | |     | |          
% % % % %    /  \  | |_ __ | |__   __ _ 
% % % % %   / /\ \ | | '_ \| '_ \ / _` |
% % % % %  / ____ \| | |_) | | | | (_| |
% % % % % /_/    \_\_| .__/|_| |_|\__,_|
% % % % %            | |                
% % % % %            |_|                

p1      = 4;    % selected parameter
window  = 3;    % * 2 + 1 = size of moving window
% Window size equal to 3 because "the concentration history is sensitive to 
% the parameters over relatively narrow time spans that are associated with 
% specific segments of the concentration history" (Wagner & Harvey 1997).
% And used also by by Wagener et al., 2002.

% 'fixed' algorithm parameters ********************************************************************************************************************
containers  = 20; % split of parameter range
aa          =  1; % [1] medium window (running mean) [2] right boundary window (regressive)
grouping    = 10; % number of 'horizontal' groups

% 'fixed' algorithm parameters in the code
ff(1)       =  p1; 	% last selected parameter equals first one
timestep    = ' ';  % time-series time step

   timestep='samples';

% DYNIA *******************************************************************
h_wait=waitbar(0,'Running DYNIA Algorithm - Alpha');  
for dt=1:t

   % calculate "likelihood" at every time step dt
   if aa==1 % medium window (running mean)
   	for no=1:size(mct,2)
   	   if dt<window+1
   	   	    residual(no,dt)=mean(abs(mct(1:dt+window,no)-obs(1:dt+window))); % mean absolute error in moving window
   	   elseif t-dt<window+1
   	        residual(no,dt)=mean(abs(mct(dt-window:t,no)-obs(dt-window:t))); % mean absolute error in moving window
   	   else
			residual(no,dt)=mean(abs(mct(dt-window:dt+window,no)-obs(dt-window:dt+window))); % mean absolute error in moving window         
   	   end      
   	end
	elseif aa==2 % right boundary window (recursive)
   	for no=1:size(mct,2)
   	    if dt<2*window+1
   	   	    residual(no,dt)=mean(abs(mct(1:dt+2*window,no)-obs(1:dt+2*window))); % mean absolute error in moving window
        else
			residual(no,dt)=mean(abs(mct(dt-2*window:dt,no)-obs(dt-2*window:dt))); % mean absolute error in moving window         
   	    end      
   	end
	end            
   
   if max(residual)~=0
      L=residual(:,dt)./max(residual(:,dt)); % normalise criterion
   else
      L=residual(:,dt);
   end
   L=1-L; % likelihood (high values indicate more likely [probable] models)
   if min(L)<0|min(L)==0, L=L-min(L)+1000*eps;end; % transform negative lhoods
   L=L';
   
   % sort data according to selected perf measure
	[I,J]=sort(L);
	pop=dat(J,:);
	cls=floor(length(pop)/grouping);
	tmx=zeros(cls,grouping);tmy=tmx;
      
	for i=p1:ff(1)
      
      % calculate cumulative distribution of top 10% of the model population ******************************************************************
      
      tm=pop(cls*(grouping-1)+1:cls*grouping,i);
      tmx=sort(tm);
	  tmy=(1:length(tmx))/cls; 
      
      % calculate the 90% confidence limits *************************************************************************************************
      
      tucl=find(tmy>0.95); ucl(dt)=tmx(tucl(1));
      tlcl=find(tmy>0.05); lcl(dt)=tmx(tlcl(1));
      
	  % calculate gradients *******************************************************************************************************************
	   
	  step=(max(dat(:,i))-min(dat(:,i)))/containers;
	  XI(:,i)=[min(dat(:,i)):step:max(dat(:,i))]';
      
      temp1=min(dat(:,i))-.001;
      if temp1<0
         temp1=0;
      end
      temp2=max(dat(:,i))+.001;    

      for kk=1:length(tmx)-2
          if tmx(kk)==tmx(kk+1)
              tmx(kk+1)=tmx(kk+1)+[tmx(kk+2)-tmx(kk)]/1000;
          end
      end
      
      [YI]=interp_special(sort([temp1;tmx;temp2]),[0 tmy 1],XI(:,i));
      
      [FX(:,dt,i)] = gradient(YI); % calculate gradient within containers
      
      % select best performing model as a function of the highest gradient
      % try select best model as point estimate!!!!
      
      location=find(FX(:,dt,i)==max(FX(:,dt,i)));
      best(dt)=XI(location,i);
      
      % find number of models within best segment (i.e. pixel)
      
      if location<size(XI,1)
          [YY,ZZ]=find(tmx>XI(location,i)&tmx<XI(location+1,i));
      else
          [YY,ZZ]=find(tmx>XI(location,i));
      end
                 
      clusterline_pixel(dt)=size(YY,1);
      clear YY ZZ
      
      % find number of models within 90% cfls
      
      [YY,ZZ]=find(tmx>lcl(dt)&tmx<ucl(dt));
      clusterline_cfls(dt)=size(YY,1);      
      clear YY
      
	  % **********************************************************************************************************************************

      dtcolor(:,dt,i)=FX(:,dt,i); % let color values run from min (bottom) to max (top) parameter value
       
      clear tm YI
       
	end
   
   clear L 
   waitbar(dt/t);
   
end
close(h_wait);
G = find(isnan(dtcolor)); dtcolor(G)=zeros(size(G));

for i=p1:ff(1) % every parameter

	infocontent=1-[(ucl-lcl)/(max(dat(:,i))-min(dat(:,i)))];
   
	nucl=ucl/(max(dat(:,i))-min(dat(:,i))); % normalise ucl
	nlcl=lcl/(max(dat(:,i))-min(dat(:,i))); % normalise lcl

end

% PLOTTING ****************************************************************

% create x matrix for patch function
% one x vector for every time step
% matrix is the same for all parameters
for dt=1:t
   xmatrix(1:2,dt)=dt-0.5;
   xmatrix(3:4,dt)=dt+0.5;
end
% create y matrix for patch function
for C=1:containers
   ymatrix(1,C)=(C-1)*(1/containers);
   ymatrix(2:3,C)=C*(1/containers);
   ymatrix(4,C)=(C-1)*(1/containers);
end

% figure
for i=p1:ff(1) % every parameter
  
  subplot(5,2,7); % DYNIA PLOT ****************************************************************************************************************************
  yyaxis left 
  for dt=1:t % every time step
     for C=1:containers % every container
		tcolor=[abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i)))) abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i)))) abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i))))];
        tcolor=[1 tcolor(1,2) tcolor(1,3)]; % let's plot them in RED shades
        patch(xmatrix(:,dt),ymatrix(:,C).*(max(XI(:,i))-min(XI(:,i)))+min(XI(:,i)),tcolor,'edgecolor',tcolor);hold on;
     end
  end
  plot(ucl,':','color','k','linewidth',1);hold on; % upper confidence limit
  plot(lcl,':','color','k','linewidth',1);hold on; % lower confidence limit
  axis([1 t min(XI(:,i)) max(XI(:,i))]);
  ylabel(['Alpha [1/s]']);
  
  yyaxis right
  plot((1:1:t),BTC_input(:,2),'k','linewidth',2);
  hold on;
  axis('auto'); 
  ylabel(['Cl [mg/l]']);  
  xlim([0 t]);
  
  set(gca,'layer','top');
  set(get(gca,'xlabel'),'fontsize',16);
  set(get(gca,'ylabel'),'fontsize',16);
  set(get(gca,'yaxis'),'Color','black');

  set(gca,'layer','top');
  set(gca,'fontsize',16,'linewidth',2);
  
  subplot(5,2,8); % INFORMATION CONTENT PLOT *************************************************************************************************************
  
  infocontent=1-[(ucl-lcl)/(max(dat(:,i))-min(dat(:,i)))];
  
  bar(infocontent,'r');hold on;
  ylim([0 1]);
  ylabel({'Inform.';'cont. \alpha []'});
  yyaxis right
  plot((1:1:t),BTC_input(:,2),'color','black','linewidth',2);
  hold on;
  axis('auto'); 
  ylabel(['Cl [mg/l]']);
  xlim([0 t]);
  
  save infocontent infocontent;
  
  set(gca,'layer','top');
  set(get(gca,'xlabel'),'fontsize',16);
  set(get(gca,'ylabel'),'fontsize',16);
  set(get(gca,'yaxis'),'Color','black');

  set(gca,'layer','top');
  set(gca,'fontsize',16,'linewidth',2);
   
end
Dynia_param_range.Alpha(:,1)=BTC_input(:,1);
Dynia_param_range.Alpha(:,2)=lcl';
Dynia_param_range.Alpha(:,3)=ucl';
Dynia_param_range.Alpha(:,4)=infocontent;
% save Dynia_Alpha.mat

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % %                      _______ _____ 
% % % % %     /\              |__   __/ ____|
% % % % %    /  \                | | | (___  
% % % % %   / /\ \               | |  \___ \ 
% % % % %  / ____ \              | |  ____) |
% % % % % /_/    \_\             |_| |_____/ 
% % % % %             ______                 
% % % % %            |______|                

p1      = 5;    % selected parameter
window  = 3;    % * 2 + 1 = size of moving window
% Window size equal to 3 because "the concentration history is sensitive to 
% the parameters over relatively narrow time spans that are associated with 
% specific segments of the concentration history" (Wagner & Harvey 1997).
% And used also by by Wagener et al., 2002.

% 'fixed' algorithm parameters ********************************************************************************************************************
containers  = 20; % split of parameter range
aa          =  1; % [1] medium window (running mean) [2] right boundary window (regressive)
grouping    = 10; % number of 'horizontal' groups

% 'fixed' algorithm parameters in the code
ff(1)       =  p1; 	% last selected parameter equals first one
timestep    = ' ';  % time-series time step

timestep='samples';

% DYNIA *******************************************************************
h_wait=waitbar(0,'Running DYNIA Algorithm - A_T_S');  
for dt=1:t

   % calculate "likelihood" at every time step dt
   if aa==1 % medium window (running mean)
   	for no=1:size(mct,2)
   	   if dt<window+1
   	   	    residual(no,dt)=mean(abs(mct(1:dt+window,no)-obs(1:dt+window))); % mean absolute error in moving window
   	   elseif t-dt<window+1
   	        residual(no,dt)=mean(abs(mct(dt-window:t,no)-obs(dt-window:t))); % mean absolute error in moving window
   	   else
			residual(no,dt)=mean(abs(mct(dt-window:dt+window,no)-obs(dt-window:dt+window))); % mean absolute error in moving window         
   	   end      
   	end
	elseif aa==2 % right boundary window (recursive)
   	for no=1:size(mct,2)
   	    if dt<2*window+1
   	   	    residual(no,dt)=mean(abs(mct(1:dt+2*window,no)-obs(1:dt+2*window))); % mean absolute error in moving window
        else
			residual(no,dt)=mean(abs(mct(dt-2*window:dt,no)-obs(dt-2*window:dt))); % mean absolute error in moving window         
   	    end      
   	end
	end            
   
   if max(residual)~=0
      L=residual(:,dt)./max(residual(:,dt)); % normalise criterion
   else
      L=residual(:,dt);
   end
   L=1-L; % likelihood (high values indicate more likely [probable] models)
   if min(L)<0|min(L)==0, L=L-min(L)+1000*eps;end; % transform negative lhoods
   L=L';
   
   % sort data according to selected perf measure
	[I,J]=sort(L);
	pop=dat(J,:);
	cls=floor(length(pop)/grouping);
	tmx=zeros(cls,grouping);tmy=tmx;
      
	for i=p1:ff(1)
      
      % calculate cumulative distribution of top 10% of the model population ******************************************************************
      
      tm=pop(cls*(grouping-1)+1:cls*grouping,i);
      tmx=sort(tm);
	  tmy=(1:length(tmx))/cls; 
      
      % calculate the 90% confidence limits *************************************************************************************************
      
      tucl=find(tmy>0.95); ucl(dt)=tmx(tucl(1));
      tlcl=find(tmy>0.05); lcl(dt)=tmx(tlcl(1));
      
	  % calculate gradients *******************************************************************************************************************
	   
	  step=(max(dat(:,i))-min(dat(:,i)))/containers;
	  XI(:,i)=[min(dat(:,i)):step:max(dat(:,i))]';
      
      temp1=min(dat(:,i))-.001;
      if temp1<0
         temp1=0;
      end
      temp2=max(dat(:,i))+.001;    

      for kk=1:length(tmx)-2
          if tmx(kk)==tmx(kk+1)
              tmx(kk+1)=tmx(kk+1)+[tmx(kk+2)-tmx(kk)]/1000;
          end
      end
      
      [YI]=interp_special(sort([temp1;tmx;temp2]),[0 tmy 1],XI(:,i));
      
      [FX(:,dt,i)] = gradient(YI); % calculate gradient within containers
      
      % select best performing model as a function of the highest gradient
      % try select best model as point estimate!!!!
      
      location=find(FX(:,dt,i)==max(FX(:,dt,i)));
      best(dt)=XI(location,i);
      
      % find number of models within best segment (i.e. pixel)
      
      if location<size(XI,1)
          [YY,ZZ]=find(tmx>XI(location,i)&tmx<XI(location+1,i));
      else
          [YY,ZZ]=find(tmx>XI(location,i));
      end
                 
      clusterline_pixel(dt)=size(YY,1);
      clear YY ZZ
      
      % find number of models within 90% cfls
      
      [YY,ZZ]=find(tmx>lcl(dt)&tmx<ucl(dt));
      clusterline_cfls(dt)=size(YY,1);      
      clear YY
      
	  % **********************************************************************************************************************************

      dtcolor(:,dt,i)=FX(:,dt,i); % let color values run from min (bottom) to max (top) parameter value
       
      clear tm YI
       
	end
   
   clear L 
   waitbar(dt/t);
   
end
close(h_wait);
G = find(isnan(dtcolor)); dtcolor(G)=zeros(size(G));

for i=p1:ff(1) % every parameter

	infocontent=1-[(ucl-lcl)/(max(dat(:,i))-min(dat(:,i)))];
   
	nucl=ucl/(max(dat(:,i))-min(dat(:,i))); % normalise ucl
	nlcl=lcl/(max(dat(:,i))-min(dat(:,i))); % normalise lcl

end

% PLOTTING ****************************************************************

% create x matrix for patch function
% one x vector for every time step
% matrix is the same for all parameters
for dt=1:t
   xmatrix(1:2,dt)=dt-0.5;
   xmatrix(3:4,dt)=dt+0.5;
end
% create y matrix for patch function
for C=1:containers
   ymatrix(1,C)=(C-1)*(1/containers);
   ymatrix(2:3,C)=C*(1/containers);
   ymatrix(4,C)=(C-1)*(1/containers);
end

% figure
for i=p1:ff(1) % every parameter
  
  subplot(5,2,9); % DYNIA PLOT ****************************************************************************************************************************
  yyaxis left
  for dt=1:t % every time step
     for C=1:containers % every container
		tcolor=[abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i)))) abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i)))) abs(1-dtcolor(C,dt,i)./max(max(dtcolor(:,:,i))))];
        tcolor=[1 tcolor(1,2) tcolor(1,3)]; % let's plot them in RED shades
        patch(xmatrix(:,dt),ymatrix(:,C).*(max(XI(:,i))-min(XI(:,i)))+min(XI(:,i)),tcolor,'edgecolor',tcolor);hold on;
     end
  end
  
  plot(ucl,':','color','k','linewidth',1);hold on; % upper confidence limit
  plot(lcl,':','color','k','linewidth',1);hold on; % lower confidence limit
  axis([1 t min(XI(:,i)) max(XI(:,i))]);
  ylabel(['A_T_S [m^2]']);

  yyaxis right
  plot((1:1:t),BTC_input(:,2),'k','linewidth',2);
  hold on;
  axis('auto'); 
  ylabel(['Cl [mg/l]']);

  xlim([0 t]);
  xlabel(['time step [' timestep ']']);

  set(gca,'layer','top');
  set(get(gca,'xlabel'),'fontsize',16);
  set(get(gca,'ylabel'),'fontsize',16);
  set(get(gca,'yaxis'),'Color','black');
  
  set(gca,'layer','top');
  set(gca,'fontsize',16,'linewidth',2);
  
  subplot(5,2,10); % INFORMATION CONTENT PLOT *************************************************************************************************************
  
  infocontent=1-[(ucl-lcl)/(max(dat(:,i))-min(dat(:,i)))];
  
  bar(infocontent,'r');hold on;
  ylim([0 1]);
  ylabel({'Inform.';'cont. A_T_S []'});
  yyaxis right
  plot((1:1:t),BTC_input(:,2),'color','black','linewidth',2);
  hold on;
  axis('auto'); 
  ylabel(['Cl [mg/l]']);
  xlim([0 t]);

  xlabel(['time step [' timestep ']']);

  set(gca,'layer','top');
  set(get(gca,'xlabel'),'fontsize',14);
  set(get(gca,'ylabel'),'fontsize',14);
  set(get(gca,'title'),'fontsize',14);
  
  set(gca,'layer','top');
  set(gca,'fontsize',16,'linewidth',2);
  set(get(gca,'yaxis'),'Color','black');
  
end
% save Dynia_A_TS.mat

% colormap(gray(10)); %OLD 

% NEwW -> Shades of red 
map = [1 0.0 0.0        % let't build a matrix of RGB code for shades of red
1 0.1 0.1
1 0.2 0.2
1 0.3 0.3
1 0.4 0.4
1 0.5 0.5
1 0.6 0.6
1 0.7 0.7
1 0.8 0.8
1 0.9 0.9
1 1 1];
colormap(map);
%
h=colorbar;
set(h,'position',[.06 .11 .005 .22]);
set(h,'Ticks',[0 1]);
set(h,'TickLabels',{'H', 'L'});
set(h,'FontSize',14)
h.Label.String = {'Likelihood -','Mean absolute error'};

Dynia_param_range.ATS(:,1)=BTC_input(:,1);
Dynia_param_range.ATS(:,2)=lcl';
Dynia_param_range.ATS(:,3)=ucl';
Dynia_param_range.ATS(:,4)=infocontent;
end

