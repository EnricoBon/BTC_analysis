# BTC_analysis
Code developed by Enrico Bonanno (2021) on the base of OTIS-MCAT software developed by Ward et al., (2017).

To properly work this code has to be saved into a "Program_files" folder directly into the local disk (C:) path.
Eg: C:\Program_Files\BTC_Analysis

Compared to the OTIS-MCAT code this script is useful for studying solute breackthrough curves without having another 
breakthrough curve available upstream as boundary condition, by including solutions for the advection-dispersion equation, and
it authomatically considers the velocity as an unknown parameters together with A (area), D (longitudinal dispersion coefficient),
alpha (exchange rate), and ATS (area of transient storage) in Transient Storage Modelling (TSM). It also has different script for defining the starting BTC to be used in OTIS 
(real measured values instead of interpolated values), it computes authomatically 8 different objective functions (RMSE, r2 (NSE), nRMSE,
logRMSE, logr2, Pearson_r2, Pearson_logr2, KGE), and it includes results and plots on Global identifiability Analysis (based on GLUE methodology,
see also Kelleher et al., 2019; Wagener et al., 2003; Wagener & Kollat, 2007; Ward et al., 2017; Wlostowski et al., 2013),
and Dynamic Identifiability Analysis (Wagener et al., 2002; Wagener & Kollat, 2007). The script indicates to the modeller the reccommended steps 
for successive TSM simulations with restricted parameter range to reduce or remove TSM parameter identifiability as indicated in Bonanno et al., (submitted).

The script has been written in MATLAB R2020a and for properly run a Latin Hypercube sampling it requires the "lhsdesign" function 
(eg. Statistics and Machine Learning Toolbox). However, exhaustive explanations are written in every script of the code, and a second 
option, using a MonteCarlo sampling, is already implemented.

Given a certain breackthrough curve (BTC) it deduces:
% PART 1 - General breackthrough curve properties
                % t99 % M1 % M1norm % mu2 % mu2norm % mu3 % mu3norm % skewness
                % skewnessnorm  % appdispersivity % appdispersion % Holdback % t05 
                % t10 % t25 % t50 % t75 % t90 % t95 % t05norm % t10norm % t25norm
                % t50norm % t75norm % t90norm % t95norm % tpeak % cpeak % cpeakNORM
                
% PART 2 - Advection-Dispersion properties
                % Best-fitting ADE in different assumptions
                %     1) v is fixed and equal to L/t_peak; Q is calculated via dilution gauging method; A is calculated = Q/v;  D is calibrated.
                %     2) v is fixed and equal to L/t_peak - Q is measured by v-notch upstream - A is fixed (Q/v) - D is calibrated. Note: Mass recovered from the injection =/= mass injected and it will be equal to Q_measured * sum of concentration 
                %     3) v is calibrated - A is calibrated - D is calibrated -> this option uses Monte Carlo or latin hypercube sampling.

% PART 3 - TSM properties - random sampling of TSM parameters (v, A, D, Alpha, ATS) via latin hypercube or MonteCarlo, 
                parameter sensitivity and uncertainty with global and dynamic identifiability analysis. After the first TSM simulation, this third part
                requires active choices by the modeller to decide, depending on the results after the first TSM simulation,
                new parameter interval for the following 2nd - 3rd - nth TSM simulations. 

The complete code list include the following scripts:

BTC_analysis -> MASTER SCRIPT

ADE_analysis -> Function to solve and deduce Advection-Dispersion properties (Part 2)

ADE_collection1 -> Needed for ADE_analysis (figures)

ADE_collection2 -> Needed for ADE_analysis (figures)

BTC_prop -> Function to deducte BTC properties (part 1)

control.ino -> File needed for OTIS

conversion -> converting Electrical conductivity into Chloride concentration / removal of 0-values / removal of background values / BTC adaptation to < 198 values

DATA.INP -> File needed for OTIS

Dynia -> Dynamic identifiability analysis (Wagener et al., 2002), Adapted script for the TSM parameters (v, A, D, Alpha, ATS) and save of 90% confidence limits and information content

echo-out -> File needed for OTIS

figure_OTIS -> Global Identifiability Analysis for the first TSM simulation

MonteCarlo -> RandomSampling and ADE solving (3rd case for PART 2 main code objectives)

nthOTIS_MonteCarlo -> RandomSampling and TSM solving for the n-th TSM simulation (new restricted parameter range defined by the modeller - discharge restrition) 

ObjFun -> "Objective function" function

OTIS.EXE -> OTIS software (Runkel, 1998)

OTIS_MonteCarlo -> RandomSampling and TSM solving for the 1st TSM simulation (large param range defined by literature and ADE solution - no discharge restrition) 

OTIS_run -> Building of the txt files to properly run OTIS software during the loop simulation for TSM

PARAMS.INP -> File needed for OTIS

params.template -> File needed for OTIS

q.template -> File needed for OTIS

Q.INP -> File needed for OTIS

soluteout.out -> File needed for OTIS

StatFigures -> Global identifiability analysis for ADE (based on RMSE)

StatFigures_var -> Global identifiability analysis for ADE (based on NSE, logRMSE, and KGE)
