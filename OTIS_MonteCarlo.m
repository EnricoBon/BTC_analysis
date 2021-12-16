function [Working_matrix,Order_Working_matrix,Order_Working_matrix_20,...
    Order_Working_matrix_10,Data] = OTIS_MonteCarlo(ADE1,ADE3,BTC_input,...
    L,dx,n,OOO,dt,M_g,BTC_result,Description1,Description2)

% Recall the best values for v, A, and D from the ADE3
% Remember the ADE3.Hyperspace_01 is already sorted followind descending
% RMSE values, so the second row already has the top results of the
% simulation
v=cell2mat(ADE3.Hyperspace_01(2,1));     
A=cell2mat(ADE3.Hyperspace_01(2,2));   
D=cell2mat(ADE3.Hyperspace_01(2,3));   

OSFLAG=OOO.OSFLAG;

Data.time=BTC_input(:,1);
Data.solute=BTC_input(:,2);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %            _____      _      ____ _______ _____  _____ 
% %           / ____|    | |    / __ \__   __|_   _|/ ____|
% %          | (___   ___| |_  | |  | | | |    | | | (___  
% %           \___ \ / _ \ __| | |  | | | |    | |  \___ \ 
% %           ____) |  __/ |_  | |__| | | |   _| |_ ____) |
% %          |_____/ \___|\__|  \____/  |_|  |_____|_____/ 
% %                                               
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%% Recall the PARAMS.INP format
% #------------------------------
% OTIS Parameter Estimation Input
% #---------------------------------------------------------
% #  PRTOPT					Format of solute output files (1=main channel only, 2=main channel and storage zone)
% #  PSTEP  [hours]			Time interval at which results are printed (check solution integration time step)
% #  TSTEP  [hours]			Integration time step (0=steady state solution, >0=time variable solution)
% #  TSTART [hour]			Simulation start time
% #  TFINAL [hour]			Simulation end time
% #  XSTART [L]				Stream distance at upstream boundary (usaully 0)
% #  DSBOUND [(L/sec)mg/L]	Downstream boundary condition describes conc. gradient at dowstream boundary (usually 0)
% #  NREACH					Number of modeled reachesPrint Option = PRTOPT

% This is another important detail in OTIS:
% If we have a slug injection we HAVE TO set the timestep equal to 1
% second. In Ward et al. 2017 there is the definition of the timestep as:
% dt/3600, where dt is the timestep fixed for the calculation.
% This might be not correct, because if dt=5s then we set the initial concentration to 
% be long a time equal 5s, which means a constant injection 5 seconds long
% this is not a detail, since it drives strong differences in terms of total
% simulated injected salt. 
% To simulate a slug injection here we always set:
dt_h=0.0003; % ALWAYS -> 1 SEC TIMESTEP 

Instate.PRTOPT=OOO.PRTOPT;  % 1=only stream; 2=stream and TS
Instate.PSTEP=dt_h;         % Print results at hours
Instate.TSTEP=dt_h;         % Calculate solution at hours
Instate.TSTART=0;           % We have slug injection -> start at 0 time
Instate.TFINAL=BTC_input((length(BTC_input(:,1))),1)*2; % Double the end time
Length=length(BTC_input(:,1));
Instate.XSTART=0;
Instate.DSBOUND=0;
Instate.NREACH=1;

% OTIS format for Print information
% #                 Print Information
% #                 for I = 1, NPRINT 
% #
% #  NPRINT		Number of locations along the stream where output is desired
% #  IOPT		Interpolation option for print locations (1=linear interpolation used, 0=nearest upstream segment value used)
% #  PRINTLOC 	Downstream distance to a given print location, repeated NPRINT times (max of 30)   
% #                            

Instate.NPRINT=1;       % We have 1 location
Instate.IOPT=0;
Instate.PRINTLOC=L;     % reach length

% OTIS format for Physical properties
% #              Physical Parameters
% #               for I = 1, NREACH
% # 
% #  NSEG				Number of x-sectional segments within the reach. Max of 5000 (RCLEN/2 seems to work well)
% #  RCHLEN(m)			Reach length   
% #  DISP(m^2/sec)	Dispersion
% #  AREA2(m^2)		Transient storage zone area (must be a non-zero value) 
% #  ALPHA(sec-1)		Transient storage exchange coefficient 

Instate.RCHLEN=round(Instate.PRINTLOC*2.2);
% Make the simulated reach twice(+) as long as your actual reach of
% interest - this keep the d/s boundary condition from having an impact
% on the reach that you are interested in -- common practice in OTIS

%set number of segments for a given reach length and spatial step
Instate.NSEG=round(Instate.RCHLEN/dx);

% The other physical properties have to be randomly generated, so let's
% skip DISP, AREA2, and ALPHA for now

% OTIS format for Solute properties
% #        Number of Solutes and flags for decay and sorption
% #
% #  NSOLUTE  	Number of solutes (max of 3 for OTIS, 1 for OTIS-P)
% #  IDECAY 	Decay option (0 = OFF)
% #  ISORB 		Sorption option (0 = OFF)

Instate.NSOLUTE=1;
% Set Decay   
if OOO.LAMMIN==0 & OOO.LAMMAX==0  & OOO.LAM2MIN==0 && OOO.LAM2MAX==0
        Instate.IDECAY=0;
 else
        Instate.IDECAY=1;
end
    
% Set Sorption 
if OOO.LAMHATMIN==0 & OOO.LAMHATMAX==0 & OOO.LAMHAT2MIN==0 & OOO.LAMHAT2MAX==0 & OOO.RHOMIN==0 & OOO.RHOMAX==0 & OOO.KDMIN==0 & OOO.KDMAX==0 & OOO.CBGRMIN==0 & OOO.CBGRMAX==0
        Instate.ISORB=0;
 else
        Instate.ISORB=1;
end

% OTIS format for the Boundary conditions
% #         Time-Variable Upstream Boundary Conditions
% #                 for I = 1, NBOUND
% #
% #  NBOUND					Define total number of upstream boundary time steps (max of 200) 
% #  IBOUND					Boundary condition option (1=concentration step profile, 2=flux step profile, 3=concentration continuous profile)
% #  USTIME(hour)				Define each time step (max of 200)    
% #  USBC(mg/l; mg/l m^3/sec) 	Solute concentration values corresponding to each USTIME (units depend on IBOUND)
% #

% We have an instantaneous concentration, therefore:

Instate.NBOUND=3;
Instate.IBOUND=1;

% dt_h this is particularly important for our definition of the slug boundary
% condition
Instate.USTIME_USBC=[0.0,0.0;         % Zero seconds -> 0 Conc
                     dt_h, 0.0;       % 1 Second -> Preallocating Slug Conc
                     2*dt_h, 0.0];    % 2 seconds -> 0 Conc

% OTIS Discharge properties (we can set just QSTEP, since Q will vary with
% our sampled v)
Instate.QSTEP=0.00;
Instate.QSTART=1;       % We pre-allocate it, but we will change the value
                        % in the OTIS_run function

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                       
% %  _           _   _         _                                     _          
% % | |         | | (_)       | |                                   | |         
% % | |     __ _| |_ _ _ __   | |__  _   _ _ __   ___ _ __ ___ _   _| |__   ___ 
% % | |    / _` | __| | '_ \  | '_ \| | | | '_ \ / _ \ '__/ __| | | | '_ \ / _ \
% % | |___| (_| | |_| | | | | | | | | |_| | |_) |  __/ | | (__| |_| | |_) |  __/
% % |______\__,_|\__|_|_| |_| |_| |_|\__, | .__/ \___|_|  \___|\__,_|_.__/ \___|
% %                                   __/ | |                                   
% %                                  |___/|_|                                                          
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                          
                        
% Latin Hypercube sampling considering A, v, D, As, Alpha
% NOTE: Since we are modifying also the velocity in terms of OTIS
% formulation we are changing the discharge as well (Q=A*v), which means
% that we will have to modify also the initial conditions in terms of
% discharge and Initial concentration at the beginning of the reach:
% Q = v*A -> in the Q.INP file
% C = Injected mass/Q -> To obtain in the second of the immission the
% initial slug concentration in mg/l --> g/m^3

% Let's build the latin hypercube
clear OTIS_hypercube_input
% 
v=cell2mat(ADE3.Hyperspace_01(2,1));     
A=cell2mat(ADE3.Hyperspace_01(2,2));   
D=cell2mat(ADE3.Hyperspace_01(2,3));  

vmin=v*0.5;     
vmax=v+v*0.5;
Amin=A*0.5;
Amax=A+A*0.5;
Dmin=0.0001;
Dmax=D*2;

Alphamin=0.00001;
Alphamax=0.1;
A2min=0.00001;           
A2max=20;

 min_ranges= [vmin; Amin; Dmin; Alphamin; A2min];     % Min values acceptable for v, A and D
 max_ranges= [vmax; Amax; Dmax; Alphamax; A2max];     % Max values acceptable for v, A and D
 p=length(min_ranges);     
 
 NN=n;    
 
slope=max_ranges-min_ranges;
offset=min_ranges;
SLOPE=ones(NN,p);
OFFSET=ones(NN,p);
for i=1:p
    SLOPE(:,i)=ones(NN,1).*slope(i);
    OFFSET(:,i)=ones(NN,1).*offset(i);
end
X_normalized = lhsdesign(NN,p);        % this one builds latin hypercube in [0 1] interval
X_scaled=SLOPE.*X_normalized+OFFSET;   % this one modify the [0 1] latin hypercube in the
                                       % latin hypercube of our parameters
OTIS_hypercube_input(:,1)=X_scaled(:,1);    % VELOCITY [m/s]
OTIS_hypercube_input(:,2)=X_scaled(:,2);    % AREA CHANNEL [m^2]
OTIS_hypercube_input(:,3)=X_scaled(:,3);    % DISP [m^2/s]
OTIS_hypercube_input(:,4)=X_scaled(:,4);    % ALPHA [1/s]
OTIS_hypercube_input(:,5)=X_scaled(:,5);    % AREA TRANS STORAGE [m^2]

OTIS_hypercube_input(:,6)=OTIS_hypercube_input(:,1).*OTIS_hypercube_input(:,2);  % DISCHARGE [m3/s] +2L used for the solution
OTIS_hypercube_input(:,7)=M_g./(OTIS_hypercube_input(:,6));    % Instantaneous Cl concentration [g/m3](or [mg/l])

% The first TSM simulation is done on a wide parameter range and with no
% discharge condition 
                                      
% build the string for the figures -> so we know the extremes of the
% sampling immediately once we see the first figure

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

str(1,1)=sprintf(formatSpec1,min_ranges(1,1));
str(2,1)=sprintf(formatSpec2,max_ranges(1,1));
str(3,1)=sprintf(formatSpec3,min_ranges(2,1));
str(4,1)=sprintf(formatSpec4,max_ranges(2,1));
str(5,1)=sprintf(formatSpec5,min_ranges(3,1));
str(6,1)=sprintf(formatSpec6,max_ranges(3,1));
str(7,1)=sprintf(formatSpec7,min_ranges(4,1));
str(8,1)=sprintf(formatSpec8,max_ranges(4,1));
str(9,1)=sprintf(formatSpec9,min_ranges(5,1));
str(10,1)=sprintf(formatSpec10,max_ranges(5,1));

% clear the workspace
clear X_normalized X_scaled SLOPE OFFSET max_ranges min_ranges p vmin vmax...
    Amin Amax Dmin Dmax Alphamin Alphamax A2min A2max A v D i NN ReferenceQ...
    Min_ReferenceQ Max_ReferenceQ

% That's why we set:
    Not_used_param.QLATIN = 0;
    Not_used_param.QLATOUT = 0;
    Not_used_param.CLATIN = 0;

    Not_used_param.LAMBDA = 0;
    Not_used_param.LAMBDA2 = 0;

    Not_used_param.LAMHAT = 0;
    Not_used_param.LAMHAT2 = 0;
    Not_used_param.RHO = 0;
    Not_used_param.KD = 0;
    Not_used_param.CSBACK = 0;

% Set the vector for the simulation time
    simtimes=[0:Instate.PSTEP:2*BTC_input((length(BTC_input(:,1))),1)+2*Instate.PSTEP];
% from 0 to 2*(last observed time)+2*dt with step equal to dt (in hours)
   
n=length(OTIS_hypercube_input(:,1));
%%%   
  
% Preallocate the info for the computed BTC
    Data.t99 = zeros(n,1);
    Data.M1 = zeros(n,1);
    Data.M1norm = zeros(n,1);
    Data.mu2 = zeros(n,1);
    Data.mu2norm = zeros(n,1);
    Data.mu3 = zeros(n,1);
    Data.mu3norm = zeros(n,1);
    Data.skewness = zeros(n,1);
    Data.skewnessnorm = zeros(n,1);
    Data.appdispersivity = zeros(n,1);
    Data.appdispersion = zeros(n,1);
    Data.Holdback = zeros(n,1);
    Data.t05 = zeros(n,1);
    Data.t10 = zeros(n,1);
    Data.t25 = zeros(n,1);
    Data.t50 = zeros(n,1);
    Data.t75 = zeros(n,1);
    Data.t90 = zeros(n,1);
    Data.t95 = zeros(n,1);
    Data.t05norm = zeros(n,1);
    Data.t10norm = zeros(n,1);
    Data.t25norm = zeros(n,1);
    Data.t50norm = zeros(n,1);
    Data.t75norm = zeros(n,1);
    Data.t90norm = zeros(n,1);
    Data.t95norm = zeros(n,1);
    Data.tpeak = zeros(n,1);
    Data.cpeak = zeros(n,1);
    Data.cpeaknorm = zeros(n,1);
% Preallocate the objective function for the computed BTC    
    Data.RMSE = zeros(n,1);
    Data.r2 = zeros(n,1);               % Remember this is mathematically equal to Nash-Sutcliff eff.
    Data.nRMSE = zeros(n,1);
    Data.logRMSE = zeros(n,1);
    Data.logr2 = zeros(n,1);
    Data.Pearson_r2 = zeros(n,1);
    Data.Pearson_logr2=zeros(n,1);
    Data.KGE=zeros(n,1);                % Kling-Gupta Efficiency

% ged rid of ANYTHING you don't need that might slower the OTIS run
clearvars -except Data Description1 Description2 Instate OTIS_hypercube_input ...
    OSFLAG Not_used_param BTC_input L n BTC_result str dt_h dx ADE1 M_g

% computational time and laptop memory are often not sufficient to run
% 100'000 OTIS simulation, that's why it's better to save the start and, in
% case on any crash, load directly the "Start" workspace and run the code
% again from this point. 
% for 100000 OTIS simulation i suggest to divide the Data and
% OTIS_hypercube_input into 3 runs of 34'000

save('C:\Program_Files\Enrico_Otis\Output_files_OTIS\Start.mat');

wbhandle=waitbar(0,'Completing LatinHypercube OTIS simulations. Please wait and cross your fingers.');

% And now the computational pain begins 
for i = 1:n             % i defines every run
    
    %Run the forward model for parameter set i 
        Model=OTIS_run(Instate,i,OTIS_hypercube_input,OSFLAG,Not_used_param);

    % Fit OBS the time of the observations with the modelled one
        Sim = interp1(Model.ttime,Model.conc_Channel,BTC_input(:,1));

    %Store model output data
%         if BTCFLAG==1
%             MOD(i,:)= Sim; %(skipped because the files get wildly large)
%         end
          
  % Evaluate BTC properties
  % create a temporary file for the function BTC_analysis
  BTC_temp(:,1)=BTC_input(:,1);
  BTC_temp(:,2)=Sim(:,1);
  
  % fix possible negative values in the computation
  for ppp=1:1:length(BTC_temp(:,2))
      if BTC_temp(ppp,2)<0
          BTC_temp(ppp,2)=0.00;
      end
  end
  
  % compute some properties of the BTC   
  BTCMOD = BTC_prop(BTC_temp,L);
                   
	%Store the absolute error for the different metrics from BTCAnalysis
        Data.t99(i) = abs(BTC_result.t99-BTCMOD.t99);
        Data.M1(i) = abs(BTC_result.M1-BTCMOD.M1);
        Data.M1norm(i) = abs(BTC_result.M1norm-BTCMOD.M1norm);
        Data.mu2(i) = abs(BTC_result.mu2 - BTCMOD.mu2);
        Data.mu2norm(i) = abs(BTC_result.mu2norm - BTCMOD.mu2norm);
        Data.mu3(i) = abs(BTC_result.mu3 - BTCMOD.mu3);
        Data.mu3norm(i) = abs(BTC_result.mu3norm - BTCMOD.mu3norm);
        Data.skewness(i) = abs(BTC_result.skewness-BTCMOD.skewness);
        Data.skewnessnorm(i) = abs(BTC_result.skewness-BTCMOD.skewnessnorm);
        Data.appdispersivity(i) = abs(BTC_result.appdispersivity-BTCMOD.appdispersivity);
        Data.appdispersion(i) = abs(BTC_result.appdispersion-BTCMOD.appdispersion);
        Data.Holdback(i) = abs(BTC_result.Holdback-BTCMOD.Holdback);
        Data.t05(i) = abs(BTC_result.t05-BTCMOD.t05);
        Data.t10(i) = abs(BTC_result.t10-BTCMOD.t10);
        Data.t25(i) = abs(BTC_result.t25-BTCMOD.t25);
        Data.t50(i) = abs(BTC_result.t50-BTCMOD.t50);
        Data.t75(i) = abs(BTC_result.t75-BTCMOD.t75);
        Data.t90(i) = abs(BTC_result.t90-BTCMOD.t90);
        Data.t95(i) = abs(BTC_result.t95-BTCMOD.t95);
        Data.t05norm(i) = abs(BTC_result.t05norm-BTCMOD.t05norm);
        Data.t10norm(i) = abs(BTC_result.t10norm-BTCMOD.t10norm);
        Data.t25norm(i) = abs(BTC_result.t25norm-BTCMOD.t25norm);
        Data.t50norm(i) = abs(BTC_result.t50norm-BTCMOD.t50norm);
        Data.t75norm(i) = abs(BTC_result.t75norm-BTCMOD.t75norm);
        Data.t90norm(i) = abs(BTC_result.t90norm-BTCMOD.t90norm);
        Data.t95norm(i) = abs(BTC_result.t95norm-BTCMOD.t95norm);        
        Data.tpeak(i) = abs(BTC_result.tpeak - BTCMOD.tpeak);
        Data.cpeak(i) = abs(BTC_result.cpeak - BTCMOD.cpeak);
        Data.cpeaknorm(i) = abs(BTC_result.cpeakNORM - BTCMOD.cpeakNORM);

        % Calcualte objective functions
        % Create a temporary CC1 vector
        CC1(:,1)=BTC_temp(:,2);
        clear BTC_temp  
        % Run the ObjFun function
        
        [RMSE,r2,nRMSE,logRMSE,logr2,Pearson_r2,Pearson_logr2,KGE] = ObjFun(BTC_input,CC1);
        clear CC1 Sim
        
        % Store the results
        Data.RMSE(i) = RMSE;
        Data.r2(i) = r2;               
        Data.nRMSE(i) = nRMSE;
        Data.logRMSE(i) = logRMSE;
        Data.logr2(i) = logr2;
        Data.Pearson_r2(i) = Pearson_r2;
        Data.Pearson_logr2(i) = Pearson_logr2;
        Data.KGE(i) = KGE;    
        
        clear RMSE r2 nRMSE logRMSE logr2 Pearson_r2 Pearson_logr2 KGE
        % get rid of the OTIS results
        clear Model BTCMOD
        
        % progress with the waitbar
        waitbar(i / n) 
end    

close(wbhandle)

save('C:\Program_Files\Enrico_Otis\Output_files_OTIS\OTIS_results.mat','Data','BTC_input','dt_h',...
    'Description1','Description2','L','n','Not_used_param','OSFLAG','OTIS_hypercube_input','dx',...
    'str','Instate');

% At this point we need to do some figures with the following function

[Working_matrix,Order_Working_matrix,Order_Working_matrix_20,Order_Working_matrix_10,...
    TEMP_BTC,Interp_curve] = figures_OTIS(OTIS_hypercube_input,Description1, Description2, ...
    Data,Instate,OSFLAG,Not_used_param,BTC_input,str);

% Dynamic identifiability Analysis

[f,Dynia_param_range] = Dynia(TEMP_BTC,Order_Working_matrix,BTC_input,Description1, Description2);

saveas(f,[pwd '/Output_files_OTIS/All_Dynia.fig']);
saveas(f,[pwd '/Output_files_OTIS/All_Dynia.tif']);

close(f)

save('C:\Program_Files\Enrico_Otis\Output_files_OTIS\OTIS_results.mat','Data','BTC_input','dt_h',...
    'Description1','Description2','L','n','Not_used_param','OSFLAG','OTIS_hypercube_input','dx',...
    'str','Working_matrix','Order_Working_matrix','Order_Working_matrix_20','Order_Working_matrix_10'...
    ,'Instate','TEMP_BTC','Interp_curve','Dynia_param_range');

end

