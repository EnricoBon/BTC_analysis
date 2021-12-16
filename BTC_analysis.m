% Code name: BTC_analysis  
% Author: Enrico Bonanno

% What does the code do?
% Given a certain breackthrough curve (BTC) it deduces:
%%%% PART 1 - General breackthrough curve properties
                % t99 % M1 % M1norm % mu2 % mu2norm % mu3 % mu3norm % skewness
                % skewnessnorm  % appdispersivity % appdispersion % Holdback % t05 
                % t10 % t25 % t50 % t75 % t90 % t95 % t05norm % t10norm % t25norm
                % t50norm % t75norm % t90norm % t95norm % tpeak % cpeak % cpeakNORM
                
%%%% PART 2 - Advection-Dispersion properties
                % Best-fitting ADE in different assumptions
                %     1) v is fixed and equal to L/t_peak; 
                %     Q is calculated via dilution gauging method;  
                %     A is calculated = Q/v;  D is calibrated
                %     2) v is fixed and equal to L/t_peak - Q is measured by
                %     v-notch upstream - A is fixed (Q/v) - D is calibrated
                %     Mass recovered from the injection =/= mass injected and it 
                %     will be equal to Q_measured*sum of concentration
                %     3) v is calibrated - A is calibrated - D is
                %     calibrated -> this option uses Monte Carlo
                %     approach/latin hypercube 

%%%% PART 3 - TSM properties - parameter sensitivity and uncertainty with 
                % dynamic identifiability analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% DATA that need to be inserted manually
  
% EC(:,1)= this is the time [s];
% EC(:,2)= this is the electrical conductivity [microS/cm];

% Change info about the experiment -> to be used in every figure to avoid
% visualization problems
Description1='XXX';
Description2='YYY';

% Convert NaCl into Cl-
M=100000;        % Injected Salt mass (mg) --> 100g=100000mg
M=M*0.6067;      % Convert Mass of NaCl into mass oc Chloride only [mg]
ts=EC(2,1)-EC(1,1);            % timestep set during the experiment [s]
Qmeasured=AAA;   % Q measured (eg: from a stream gauge upstream) [l/s]

% % % y=m*EC+c  Calibration equation to convert EC[microSiemens/cm] in Cl[mg/l] 
% % % To be decucted in the lab
m=XX.XXXX;
c=-ZZ.ZZZZ;
% % %  
L=XXX;          % Distance between EC sensor and salt injection point [m]

n=50000;        % number of parameter sets via Latin Hypercube or MonteCarlo sampling                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for OTIS simulation -> Formulation of Ward et al., 2017 - OTIS-MCAT
dt=ts;          % Time step to be used in the calculation [s]
dx=0.25;        % Spatial step to be used in the calculation [m] 
%
OOO.PRTOPT=1;   % Format of solute output files (1=main channel only, 2=main channel and storage zone)
%               % OTIS model function will only work with PRTOPT set to 1
%               
% Discharge parameters
OOO.QLATINMIN=0;    % Min lateral inflows for the reach (m3/s)
OOO.QLATINMAX=0;    % Max lateral inflows for the reach (m3/s)
OOO.QLATOUTMIN=0;   % Min lateral outflows for the reach (m3/s)
OOO.QLATOUTMAX=0;   % Max lateral outflows for the reach (m3/s)
OOO.CLATMIN=0;      % Min lateral inflow C for the reach (g/m3)
OOO.CLATMAX=0;      % Max lateral inflow C for the reach (g/m3)
% Decay parameters 
OOO.LAMMIN=0;       % Min value for lambda: in-stream first order decay (s^-1)
OOO.LAMMAX=0;       % Max value for lambda: in-stream first order decay (s^-1)
OOO.LAM2MIN=0;      % Min value for lambda2: storage zone first order decay
OOO.LAM2MAX=0;      % Max value for lambda2: storage zone first order decay
% Sorption parameters 
OOO.LAMHATMIN=0;    % Min value for lambda hat (s-1)
OOO.LAMHATMAX=0;    % Max value for lambda hat (s-1)
OOO.LAMHAT2MIN=0;   % Min value for lambda hat s (s-1)
OOO.LAMHAT2MAX=0;   % Max value for lambda hat (s-1)
OOO.RHOMIN=0;       % Min rho (g/m3)
OOO.RHOMAX=0;       % Max rho (g/m3)
OOO.KDMIN=0;        % Min Kd (m3/g)      
OOO.KDMAX=0;        % Max Kd (m3/g) 
OOO.CBGRMIN=0;      % Min C_background (g/m3)
OOO.CBGRMAX=0;      % Max C_background (g/m3)

OOO.OSFLAG=1;      % equal to 1 for Windows; Equal to 2 for UNIX/LINUX    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create our output folders
if exist('Output_files')==0
               mkdir('Output_files')
end
if exist('Output_files_OTIS')==0
               mkdir('Output_files_OTIS')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function makes a conversion from my series of EC to a time series of
% data with good unit of measurements and number of observation good for the OTIS
% BTC_input(:,1) -> Hours
% BTC_input(:,2) -> g/m^3
% Compared to OTIS-MCAT this conversion does not use interp function,
% since this function would return interpolated data, not measured one. 

[BTC_input] = conversion(m,c,EC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General breackthrough curve properties
 
[BTC_result]=BTC_prop(BTC_input,L);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ADE1, ADE2, ADE3, ADE3_limits,Hyperspace,M_g,Q1,time]=ADE_analysis(BTC_input,L, M, Qmeasured,Description1,Description2, n,BTC_result);                                           

% Save the final output obtained so far
% 
 save('Output_files\ADE123_results.mat','ADE1','ADE2','ADE3','ADE3_limits',...
    'Hyperspace','Description1','Description2','EC','M','M_g','Q1','Qmeasured','L','m','c','BTC_result',...
    'BTC_input','time','dt','dx','OOO','n');

clear all

load ('Output_files\ADE123_results.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=cell2mat(ADE3.Hyperspace_01(2,1));     

% Check courant condition
Courant=(v.*dt)./dx;
if Courant>1
        disp('CFL ERROR - CFL>1 - I change the dx to make it larger')
        dx=round((v*dt)/0.5,1);
        CourantNew=(v.*dt)./dx;
end

[Working_matrix,Order_Working_matrix,Order_Working_matrix_20,...
    Order_Working_matrix_10,Data] = OTIS_MonteCarlo(ADE1,ADE3,BTC_input,L,...
    dx,n,OOO,dt,M_g,BTC_result,Description1,Description2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OTIS_MonteCarlo function already saves in the folder "Output_files_OTIS"
% the global sensitivity analysis results and the DYNIA results.
% It's really unrealistic to achieve identifiability of all the 5 TSM
% parameters after the first TSM simulation. Let's run a 2nd (and probably
% a 3rd, 4th, 5th...) with restricted parameter intervals:
clear all

load ('Output_files\ADE123_results.mat')
load ('Output_files_OTIS\OTIS_results.mat')

% The following choices depend on the modeller
% my suggestions -> 
% IF the pdf of the top 0.1, 0.5, and 1% of the results
% are defined in a range between the range of larger performances (5%, 10%
% and 20%) then the parameter range can be defined as the top 1% of the
% results:
% EG:
top1perc=length(Order_Working_matrix)*0.01;
vmin=min(Order_Working_matrix_10(1:top1perc,1));
vmax=max(Order_Working_matrix_10(1:top1perc,1)); 
Amin=min(Order_Working_matrix_10(1:top1perc,2)); 
Amax=max(Order_Working_matrix_10(1:top1perc,2)); 
Dmin=min(Order_Working_matrix_10(1:top1perc,3)); 
Dmax=max(Order_Working_matrix_10(1:top1perc,3));
Alphamin=min(Order_Working_matrix_10(1:top1perc,4)); 
Alphamax=max(Order_Working_matrix_10(1:top1perc,4));
A2min=min(Order_Working_matrix_10(1:top1perc,5)); 
A2max=max(Order_Working_matrix_10(1:top1perc,5));

% IF the pdf of the top 0.1, 0.5, and 1% of the results COINCIDE or are very
% close to the upper/lower limit of the parameter itself... we are probably
% missing the real range of the performance peak, or underestimating the
% optimal parameter distribution. I suggest to increase the parameter range
% of the "constrained" part.
% EG:
top1perc=length(Order_Working_matrix)*0.01;
vmin=min(Order_Working_matrix_10(1:top1perc,1))*0.9;
vmax=max(Order_Working_matrix_10(1:top1perc,1))*1.1; 
Amin=min(Order_Working_matrix_10(1:top1perc,2))*0.9; 
Amax=max(Order_Working_matrix_10(1:top1perc,2))*1.1;  
Dmin=min(Order_Working_matrix_10(1:top1perc,3))*0.9; 
Dmax=max(Order_Working_matrix_10(1:top1perc,3))*1.1; 
Alphamin=min(Order_Working_matrix_10(1:top1perc,4))*0.9; 
Alphamax=max(Order_Working_matrix_10(1:top1perc,4))*1.1; 
A2min=min(Order_Working_matrix_10(1:top1perc,5))*0.9; 
A2max=max(Order_Working_matrix_10(1:top1perc,5))*1.1; 

% IF Alpha and A2 are not globally identifiable, let's check their
% identifiability on the falling limb/tail of the BTC. In this case I
% suggest:

[val, idx] = find(BTC_input==max(BTC_input(:,2)));

% Consider only the falling limb / tail for alpha and Ats parameter
% identifiability

% If needed for Alpha:
k=1;
for i=val:1:length(Dynia_param_range.Alpha(:,2))
    if Dynia_param_range.Alpha(i,4)>0.66
        limitAlpha(k,:)=Dynia_param_range.Alpha(i,:);
        k=k+1;
    end
end
Alphamin=min(limitAlpha(:,2)); % Or Alphamin=min(Order_Working_matrix_10(1:top1perc,4))
Alphamax=max(limitAlpha(:,3));

% If needed for ATS:
k=1;
for i=val:1:length(Dynia_param_range.ATS(:,2))  
    if Dynia_param_range.ATS(i,4)>0.66          % information content treshold
        limitATS(k,:)=Dynia_param_range.ATS(i,:);
        k=k+1;
    end
end
A2min=0.00001; % Or A2min=min(Order_Working_matrix_10(1:top1perc,5))
A2max=max(limitATS(:,3));

% Once the modeller chose the min/max interval for the next run for each parameter:
if exist('Output_files_OTIS')==0
               mkdir('Output_files_OTIS_nth')
end

[Working_matrix,Order_Working_matrix,Order_Working_matrix_20,...
    Order_Working_matrix_10,Data] = nthOTIS_MonteCarlo(ADE1,ADE3,BTC_input,...
    L,dx,n,OOO,dt,M_g,BTC_result,Description1,Description2, vmin, vmax, Amin, Amax,...
    Dmin, Dmax, Alphamin, Alphamax, A2min, A2max);

% This function authomatically saves the results of the n-th OTIS
% simulation. The modeller might decide to run in loop this last part,
% redefining the TSM parameter intervals using both global sensitivity
% analysis and Dynamic identifiability analysis until satisfied by the
% results (EG: when the pdf of the top 0.1, 0.5, and 1% of the results are 
% are defined in a range between the range of larger performances (5%, 10% 
% and 20%), for each TSM parameter).
% Manually save the results of the 2-nd, 3-rd, n-th iteration in another
% folder before running the "nthOTIS_MonteCarlo" function, otherwise the
% results will be overwritten

