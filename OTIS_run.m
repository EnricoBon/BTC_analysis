function Model=OTIS_run(Instate,i,OTIS_hypercube_input,OSFLAG,Not_used_param)

% Modified from OTIS-MCAT (Ward et al., 2017)

PRTOPT  = Instate.PRTOPT;                          
PSTEP   = Instate.PSTEP;                
TSTEP   = Instate.TSTEP;                 
TSTART  = Instate.TSTART;                      
TFINAL  = Instate.TFINAL;                        
XSTART  = Instate.XSTART;                    
DSBOUND = Instate.DSBOUND;
 
NPRINT  = Instate.NPRINT;
IOPT    = Instate.IOPT;
PRINTLOC= Instate.PRINTLOC;
 
NREACH  = Instate.NREACH;           
RCHLEN  = Instate.RCHLEN;
NSEG    = Instate.NSEG;
 
NSOLUTE = Instate.NSOLUTE;
IDECAY  = Instate.IDECAY;
ISORB   = Instate.ISORB;
 
NBOUND  = Instate.NBOUND;
IBOUND  = Instate.IBOUND;
USTIME_USBC = Instate.USTIME_USBC;

% Fix the USTIME_USBC 
USTIME_USBC(2,2)=OTIS_hypercube_input(i,7);
USTIME_USBC(2,2)=round(USTIME_USBC(2,2),2);
 
QSTEP   = Instate.QSTEP;
QSTART  = Instate.QSTART;

% Fix the Discharge 
QSTART  = OTIS_hypercube_input(i,6);

% Summon the i-th set of the Latin Hypercube parameters we need
AREA  = OTIS_hypercube_input(i,2);
DISP  = OTIS_hypercube_input(i,3);
ALPHA = OTIS_hypercube_input(i,4);
AREA2 = OTIS_hypercube_input(i,5);

% And the not-used one as well...
    QLATIN = Not_used_param.QLATIN(1,1);
    QLATOUT = Not_used_param.QLATOUT(1,1);
    CLATIN = Not_used_param.CLATIN(1,1);

    LAMBDA = Not_used_param.LAMBDA(1,1);
    LAMBDA2 = Not_used_param.LAMBDA2(1,1);

    LAMHAT = Not_used_param.LAMHAT(1,1);
    LAMHAT2 = Not_used_param.LAMHAT2(1,1);
    RHO = Not_used_param.RHO(1,1);
    KD = Not_used_param.KD(1,1);
    CSBACK = Not_used_param.CSBACK(1,1);

% Modify the PARAMS:INP file -> This comes directly from OTIS_Mcat

template=fopen('params.template', 'rt');
params=fopen('PARAMS.INP','wt');

for k=1:107;
    tline = fgetl(template);
        if ~ischar(tline),   break,   end  ;    
    switch k;
        case 19;
            fprintf(params, '%1i\n', PRTOPT);
        case 20;
            fprintf(params, '%9f\n', PSTEP);
        case 21;
            fprintf(params, '%9f\n', TSTEP);
        case 22;
            fprintf(params, '%9f\n', TSTART);
        case 23;
            fprintf(params, '%9f\n', TFINAL);
        case 24;
            fprintf(params, '%9f\n', XSTART);
        case 25;
            fprintf(params, '%9f\n', DSBOUND);
        case 26;
            fprintf(params, '%1i\n', NREACH);
        case 41;
            fprintf(params, '%5i%13i%13f%13f%13f\n', NSEG, RCHLEN, DISP, AREA2, ALPHA);
        case 53;
            fprintf(params, '%5i%5i%5i\n', NSOLUTE, IDECAY, ISORB);
        case 63
            if IDECAY~=0
                fprintf(params, '%13f%13f\n', LAMBDA, LAMBDA2);
            end
        case 73
            if ISORB~=0
                fprintf(params, '%13i%13f%13f%13f%13f\n', LAMHAT, LAMHAT2, RHO, KD, CSBACK);
            end            
        case 85;
            fprintf(params, '%5i%5i\n', NPRINT, IOPT);
        case 89;
            fprintf(params, '%f\n', PRINTLOC);
        case 103;
            fprintf(params, '%5i%5i\n', NBOUND, IBOUND);
        case 107;
            for h=1:length(USTIME_USBC);
                fprintf(params, '%13i%13i\n', USTIME_USBC(h,:)); % slightly modified from '%13f%13f\n' to '%13i%13i\n'
            end;
        otherwise ;   
            fprintf(params, '%s\n', tline);
    end;
end;

fclose(template);
fclose(params);    
 
% Modify the Q.INP text file
template=fopen('q.template','rt');
discharge=fopen('Q.INP','wt');

for k=1:31
    tline = fgetl(template);
        if ~ischar(tline),   break,   end ;     
    switch k;
        case 13;
            fprintf(discharge, '%5i\n', QSTEP);
        case 17;
            fprintf(discharge, '%5f\n', QSTART);
        case 31;
            fprintf(discharge, '%13f%13f%13f%13f%\n',QLATIN,QLATOUT,AREA,CLATIN);
        otherwise;    
            fprintf(discharge, '%s\n', tline);
    end
end

fclose(template);
fclose(discharge);  


% Execute the OTIS software (it has to be in the same folder)

    %FOR A PC
    if OSFLAG==1
        !OTIS.EXE
    end

    %FOR A UNIX/LINUX
    if OSFLAG==2
        !./otis
    end

% % % % % % % % % % % % % % %     
% Load the OTIS output of the soluteout.out txt
% soluteout.out(:,1) time [hrs]
% soluteout.out(:,2) conc in the stream channel [g/m3]
% soluteout.out(:,3) conc in the TS zone [g/m3]
    %a=column #1: time in hrs
    %b=column #2: stream solute concentration 
    %c=column #3: transient storage zone concentration

[a, b] = textread('soluteout.out','%s %s') ;    
% [a, b, c] = textread('soluteout.out','%s %s %s') ;     % if we want TS  
% Pre-allocation of our BTC curves    
ttime =  zeros(1,length(a));
conc_Channel = zeros(1,length(b));
% conc_TS  = zeros(1,length(c));  
    
for pp = 1:length(a);
    atest = isstrprop(a{pp,1}, 'alpha');
    btest = isstrprop(b{pp,1}, 'alpha');
%     ctest = isstrprop(c{pp,1}, 'alpha');
    
    if sum(atest) > 0.0;
        ttime(1,pp) = sscanf(a{pp,1}, '%13E');  
    else
        ttime(1,pp) = sscanf('0','%13E');
    end

    if sum(btest) > 0.0;
        conc_Channel(1,pp) = sscanf(b{pp,1},'%13E');
    else
        conc_Channel(1,pp) = sscanf('0','%13E');
    end
    
%     if sum(ctest) > 0.0;
%         conc_TS(1,pp) = sscanf(c{pp,1},'%13E');
%     else
%         conc_TS(1,pp) = sscanf('0','%13E');
%     end
end    
    
Model.ttime=ttime;
Model.conc_Channel=conc_Channel;
% Model.conc_TS=conc_TS;
    
end

