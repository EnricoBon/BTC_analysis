function BTC_result=BTC_prop(BTC_input,L)

% From a rielaboration of the OTISMCAT code from Ward et al.(2017)

time=BTC_input(:,1);
conc=BTC_input(:,2);

if isnan(conc(1,1))     % This option is inserted for the OTIS 
    BTC.t99=NaN;
    BTC.M1=NaN;
    M1norm=NaN;
    BTC.M1norm=NaN;
    BTC.mu2=NaN;
    BTC.mu2norm=NaN;
    BTC.mu3=NaN;
    BTC.mu3norm=NaN;
    BTC.skewness=NaN;
    BTC.skewnessnorm=NaN;
    BTC.appdispersivity=NaN;
    BTC.appdispersion=NaN;
    BTC.Holdback=NaN;
    BTC.t05=NaN;
    BTC.t10=NaN;
    BTC.t25=NaN;
    BTC.t50=NaN;
    BTC.t75=NaN;
    BTC.t90=NaN;
    BTC.t95=NaN;
    BTC.t05norm=NaN;
    BTC.t10norm=NaN;
    BTC.t25norm=NaN;
    BTC.t50norm=NaN;
    BTC.t75norm=NaN;
    BTC.t90norm=NaN;
    BTC.t95norm=NaN;
    BTC.tpeak=NaN;
    BTC.cpeak=NaN;
    BTC.cpeakNORM=NaN; 
    
else
 
%Calculate t99
    %cumulative sum of the conc (normalized from 0 to 1)
        junk=cumtrapz(time,conc);
        junk=junk./junk(end);
            
    %pick the two points that bound the point of interest
        lim=0.99;
        indexLOW=find(junk<lim);
        indexHIGH=find(junk>lim);
        
        % linear interpolation to calc t99. 
        % interp between two junks and their corresponding time value, to 
        % find the linear-interpolated time in between at the defined value "lim"
        % EG: interp1([(0.989) (0.991)],[(1.1) (1.2)],0.99) -> 0.99junk is at 1.15h.
        BTC.t99=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[time(max(indexLOW)),time(min(indexHIGH))],lim);
        clear junk lim
        
    %clip time and conc to values of t<=t99
        cclip=conc(time<=BTC.t99);
        tclip=time(time<=BTC.t99);
        
    % Imagine our modelled BTC really sucks and the t<99 is when the tailing decreases
    % very fast. We might have just 1 value of cclip and tclip -> the
    % integration is not possible. Let's fix it with something:
        if length(cclip(:,1))<2
            cclip(2,1)=cclip(1,1)+0.0003; % +1second
            tclip(2,1)=tclip(1,1)+0.001;
        end
    % Please note that this new time series (tclip) is going to substitute time
    % So when we calculate the 750, for example, it is done NOT on
    % the total time series, as for the t99, but on the new time series
    % that goes from 0 to tclip, which is t99. So the t50 calculated using
    % tclip as maximum limit will be lower than the one calculated if we consider the
    % entire recorded time. However t99 is a meaningful truncation time so all our
    % calculations are done for times<t99 (Ward et al., 2013).
        
    %Make a normalized concentration timeseries
    %zeroth temporal moment
        M0=cumtrapz(tclip,cclip);
        
    %do the normalization
        cnorm=cclip./M0(end);
                
    %First Temporal Moment (M1)
        M1=cumtrapz(tclip,cclip .* tclip.^1);
        BTC.M1=M1(end);
        
        M1norm=cumtrapz(tclip,cnorm .* tclip.^1);
        BTC.M1norm=M1norm(end);

    %Second Central Moment (mu2)
        mu2=cumtrapz(tclip,cclip .* (tclip-BTC.M1).^2);
        BTC.mu2=mu2(end);
        
        mu2norm=cumtrapz(tclip,cnorm .* (tclip-BTC.M1norm).^2);
        BTC.mu2norm=mu2norm(end);

    %Third Central Moment (mu3)
        mu3=cumtrapz(tclip,cclip .* (tclip-BTC.M1).^3);
        BTC.mu3=mu3(end);
        
        mu3norm=cumtrapz(tclip,cnorm .* (tclip-BTC.M1norm).^3);
        BTC.mu3norm=mu3norm(end);
        
	%Skewness
        BTC.skewness=BTC.mu3./(BTC.mu2.^(3/2));
        
        BTC.skewnessnorm=BTC.mu3norm./(BTC.mu2norm.^(3/2));
        
    %Apparant Dispersivity
        BTC.appdispersivity=BTC.mu2norm.*L./2;
        
    %Apparant Dispersion
        BTC.appdispersion=(BTC.mu2norm.*L.^2)./(2.*BTC.M1norm);

    %Holdback Function
        BTC.Holdback=interp1(tclip./BTC.M1norm,cumtrapz(tclip,cnorm),1);

    % times for different percentiles
    
        %RAW
            junk=cumtrapz(tclip,cclip);
            junk=junk./junk(end);
        
            %5
            	lim=0.05;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t05=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);
        
            %10
            	lim=0.10;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t10=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);                
                
            %25
            	lim=0.25;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t25=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);
                
                
            %50
            	lim=0.50;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t50=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);                
                
            %75
            	lim=0.75;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t75=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);                
                
                
            %90
            	lim=0.90;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t90=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);                
                
            %95
            	lim=0.95;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t95
                BTC.t95=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);                
                
        %NORMALIZED
            junk=cumtrapz(tclip,cnorm);
        
            %5
            	lim=0.05;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t05norm=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);
        
            %10
            	lim=0.10;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t10norm=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);                
                
            %25
            	lim=0.25;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t25norm=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);
                
                
            %50
            	lim=0.50;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t50norm=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);                
                
            %75
            	lim=0.75;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t75norm=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);                
                
                
            %90
            	lim=0.90;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t90norm=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);                
                
            %95
            	lim=0.95;
                %pick the two points that bound the point of interest
                indexLOW=find(junk<lim);
                indexHIGH=find(junk>lim);
        
                %linear interpolation to calc t99
                BTC.t95norm=interp1([junk(max(indexLOW)),junk(min(indexHIGH))],[tclip(max(indexLOW)),tclip(min(indexHIGH))],lim);                
                 
	%Peak time
        junk=tclip(cclip==max(cclip));
        
        %if this peak value has more than one timestamp, use the mean
            if length(junk)>1
                junk=mean(junk);
            end
        
        %store it
        BTC.tpeak=junk;
            
	%Peak conc
        junk=max(cclip);
        
        %if this peak value has more than one timestamp, use the mean
            if length(junk)>1
                junk=mean(junk);
            end
        
        %store it
        BTC.cpeak=junk;
        
	%Peak conc NORM
        junk=max(cnorm);
        
        %if this peak value has more than one timestamp, use the mean
            if length(junk)>1
                junk=mean(junk);
            end
        
        %store it
        BTC.cpeakNORM=junk;       

end

BTC_result=BTC;

end

