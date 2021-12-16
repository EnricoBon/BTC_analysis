function [BTC_input] = conversion(m,c,EC)

%   We wanto to obtain the Chloride concentration time series starting from
%   the electrical conductivity time series.
%   This function makes a conversion from series of EC to a time series 
%   measurements to a <200 number of observation. That is because we need to
%   have number consistency between the number observation and the OTIS code. 
%   OTIS can run up to 200 observations, so everything need to be scaled at <200
%   observations and with specific units (conc= g/m3) and (time= hurs)

clear CLmin CL_net CL 
Length=length(EC(:,1));

for i=1:1:Length
    CL(i,1)=EC(i,1);            % time in the first column [s]
    CL(i,2)=EC(i,2)*m+c;        % obtain chloride conc in [mg/l]
end

CLmin=min(CL(:,2));
    
clear i
CL_net(:,1)=CL(:,1);            % time for the first column [s]
for i=1:1:Length
    CL_net(i,2)=CL(i,2)-CLmin;  % remove the background conc Cl_net=Ci-C0 % mg/l
end

% This part is pretty important, we need to have <200 values (best equal or
% less than 198 values) to run the TSM via OTIS.
% What values shall we keep and which one shall we remove?
% Moreover, even if we are below or equal 198 measurements but we rapidly
% decrease to the background concentration, there is no need to have a long
% measurement which increases the overall computational time.

% Step one -> recognize we can have some values that are still background
% concentration, but slightly above the minimum. This choices have to been
% accurately considered by the modeller. 
% Here there is a suggested solution:

clear beginning_rise
% Find the beginning of the rising limb of the BTC

% Method 1 - Visual check: open the CL_net matrix and observe the
% background concentration before the arrival of the BTC. Set background
% obscillation = 0 by hand. Values on the tail below the maximum background
% obscillation value should also be set = 0.

% Method 2 - Authomatised background removal 
% If we have 4 consecutive concentration increase -> The BTC arrived. Take
% the maximum concentration value before the arrival of the BTC.and remove
% it from the measured BTC.

[FX] = gradient(EC(:,2));
for i=1:1:length(FX(:,1))
    if FX(i,1)+FX(i+1,1)+FX(i+2,1)+FX(i+3,1)>0.1    
        iii=i;
        break
    end
end

% find the max value of the background obscillation  
 CLmin=max(CL_net(1:iii,2));

% new NET chloride concentration
for i=1:1:Length
    CL_net(i,2)=CL_net(i,2)-CLmin;  % remove again the max background 
    if CL_net(i,2)<0                % remove negative values 
       CL_net(i,2)=0; 
    end
end

% Remove excessive tailing (too many values at the end = 0)
if CL_net(Length,2)~=0
    endBTC=NaN;
else
for i=Length:-1:4
    if CL_net(i,2)~=0 % 
        endBTC=i;
    break
    end
end
end

% start: If th modeller would like to have some 0 values at the end of the BTC:
% Option 1: leave some 0-values
% eg: take still 10% more between the end of the BTC and the end of the
% measurements

% truncationEND=round(endBTC+(Length-endBTC)/10,0);

% Option 2: remove all 0-values at the end
truncationEND=endBTC+1;

% remove excessive tailing, if needed 
if truncationEND>0
CL_net(truncationEND:Length,:)=[];
end
% end removal 0- values at the end

% Remove 0-values at the beginning of the measurements
% start
for i=1:1:Length
    if CL_net(i,2)~=0 % 
        truncationBEG=i-2;
    break
    end
end

% Remove 0-values 

CL_net(1:truncationBEG,:)=[];
% end

% In this way we still have all the "maximum" amount of real values without
% the need to approximate or increase the timestep of the recorded BTC
% So far, we just removed 0-values, keeping intact the original BTC
% timestep definition... However, we might still have a number of
% measurements above the needed 198 maximum values for OTIS 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear STEP
if length(CL_net(:,1))>198
STEP=length(CL_net)/198;
STEP=floor(STEP);               % Round at the min so we got more values
else
    STEP=1;
end

clear BTC_input
k=1; 
i=1;

for i=1:STEP:length(CL_net)
    % OTIS input
    BTC_input(k,1)=CL_net(i,1)/3600;        % we need hours for OTIS
    BTC_input(k,2)=CL_net(i,2);             % Concentration in g/m^3 (1 g/m^3==1 mg/l)
    k=k+1;
end
clear k

% We might still have a number of measurements above 198

% Let's "break" the BTC
clear BTC_rise BTC_fall BTC_tail BTC_rise_fin BTC_fall_fin BTC_tail_fin
if length(BTC_input(:,2))>198
    
[val, idx] = find(BTC_input==max(BTC_input(:,2)));

% Despite this code is meant for slug injection it might happen
% we have a small plateau on the top ("val" has several values). Let's
% build a "fake" peak that do not ruin our analysis but allows an
% automatization for the BTC subdivision

if length(val(:,1))>1
    Test=length(val(:,1));
    AAAA=mod(Test,2);   % odd or even?
    if AAAA==0          % If even
        ggh=Test/2;
        VAL=val(ggh,1);
        BTC_input(VAL,2)=BTC_input(VAL,2)+0.0001;   
    else                % If odd
        ggh=Test/2+0.5;
        VAL=val(ggh,1);
        BTC_input(VAL,2)=BTC_input(VAL,2)+0.0001;   
    end
end

% % %

clear val idx
[val, idx] = find(BTC_input==max(BTC_input(:,2)));

BTC_rise=BTC_input(1:val,:);    % Rising part of the curve

clear diff DDD closest
diff=(BTC_input(val,2)-BTC_input(1,2))/5;
DDD=abs(BTC_input(val:length(BTC_input(:,1)),:)-diff);
[vall,closest]=min(DDD(:,2));

BTC_fall=BTC_input(val+1:val+closest,:);    % Falling part of the curve

Length=length(BTC_input(:,1));
BTC_tail=BTC_input(val+closest+1:Length,:); % Tailing part of the curve

clear STEP
k=1; 
i=1;
if length(BTC_rise(:,1))>50
   STEP=length(BTC_rise)/50;
   STEP=ceil(STEP);                     % Round at the max so we are sure we stay below 50
   for i=1:STEP:length(BTC_rise)
    BTC_rise_fin(k,1)=BTC_rise(i,1);       
    BTC_rise_fin(k,2)=BTC_rise(i,2);             
    k=k+1;
   end
else
    BTC_rise_fin=BTC_rise;
end
clear STEP
k=1; 
i=1;
if length(BTC_fall(:,1))>50
   STEP=length(BTC_fall)/50;
   STEP=ceil(STEP);                     % Round at the max so we are sure we stay below 50
   for i=1:STEP:length(BTC_fall)
    BTC_fall_fin(k,1)=BTC_fall(i,1);       
    BTC_fall_fin(k,2)=BTC_fall(i,2);             
    k=k+1;
   end
else
 BTC_fall_fin=BTC_fall;    
end

clear STEP
k=1; 
i=1;

final_limit1=length(BTC_fall_fin(:,1))+length(BTC_rise_fin(:,1));
final_limit=198-final_limit1;

if length(BTC_tail(:,1))>final_limit
   STEP=length(BTC_tail)/final_limit;
   STEP=ceil(STEP);                     % Round at the max so we are sure we stay below 98
   for i=1:STEP:length(BTC_tail)
    BTC_tail_fin(k,1)=BTC_tail(i,1);       
    BTC_tail_fin(k,2)=BTC_tail(i,2);             
    k=k+1;
   end
else
 BTC_tail_fin=BTC_tail;
end

% let's merge all together

clear BTC_input

BTC_input=vertcat(BTC_rise_fin,BTC_fall_fin,BTC_tail_fin);

end

% This way to arrange the BTC has the advantage to use only meaured values,
% however it might change the timestep between rising limb, falling limb
% and BTC tail. This is not a problem in terms of model performances, but
% it can probably give a visual hump between falling limb and tail of the
% BTC during the DYNIA simulation since DYNIA algorithm uses timesteps and 
% it therefore visually "merges" different timesteps. This problem is only 
% visual and does not affect DYNIA performances 

% Last condition:

if BTC_input(1,1)==0
    BTC_input(1,1)=0.001;
end
    
    
end


