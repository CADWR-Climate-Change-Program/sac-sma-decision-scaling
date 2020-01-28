function [surf_tot, base_tot, tet_tot, uztwc_tot, uzfwc_tot, lztwc_tot, lzfpc_tot, lzfsc_tot, adimc_tot] = sma_sacsma(pet,snow_outflow,Par, IniState)

%% ----------------------------SAC_SMA PARAMETERS------------------------------

% Capacity Thresholds
uztwm  =  Par(1);    % Upper zone tension water capacity [mm]
uzfwm  =  Par(2);    % Upper zone free water capacity [mm]
lztwm  =  Par(3);    % Lower zone tension water capacity [mm]
lzfpm  =  Par(4);    % Lower zone primary free water capacity [mm]
lzfsm  =  Par(5);    % Lower zone supplementary free water capacity [mm]
% Recession Parameters
uzk    =  Par(6);    % Upper zone free water lateral depletion rate [1/day]
lzpk   =  Par(7);    % Lower zone primary free water depletion rate [1/day]
lzsk   =  Par(8);    % Lower zone supplementary free water depletion rate [1/day]
% Percolation
zperc  =  Par(9);    % Percolation demand scale parameter [-]
rexp   =  Par(10);   % Percolation demand shape parameter [-]: exponent of the percolation equation
pfree  =  Par(11);   % Percolating water split parameter (decimal fraction): Fraction of water percolating from upper zone directly to lower zone free water storage
% Impervious area
pctim  =  Par(12);   % Impervious fraction of the watershed area (decimal fraction)
adimp  =  Par(13);   % Additional impervious areas (decimal fraction)
% Others
riva   =  Par(14);   % Riparian vegetation area (decimal fraction)
side   =  Par(15);   % The ratio of deep recharge to channel base flow [-]
rserv  =  Par(16);   % Fraction of lower zone free water not transferrable to lower zone tension water (decimal fraction)



%% ------------------Execute model for every time step---------------------

% Initial Storage States
% SAC-SMA
uztwc = IniState(1);     % Upper zone tension water content
uzfwc = IniState(2);     % Upper zone free water centent
lztwc = IniState(3);     % Lower zone tension water centent
lzfsc = IniState(4);     % Lower zone supplementary free water centent
lzfpc = IniState(5);     % Upper zone primary free water centent
adimc = IniState(6);     % Additional impervious area centent
 

% RESERVOIR STATE ARRAY INITIALIZATION
% Upper zone states
uztwc_tot = NaN(size(pet));   % State of Upper zone tension water storage [mm]
uzfwc_tot = NaN(size(pet));   % State of Upper zone free water storage [mm]
% Lower zone states
lztwc_tot = NaN(size(pet));   % State of Lower zone tension water storage [mm]
lzfsc_tot = NaN(size(pet));   % State of Lower zone free water supplementary storage [mm]
lzfpc_tot = NaN(size(pet));   % State of Lower zone free water primary storage [mm]
% Additional impervious zone states
adimc_tot = NaN(size(pet));   % State of additional impervious area storages [mm]

% MODEL OUPUT ARRAY INITIALIZATION
simflow  = NaN(size(pet));  % Simulated Streamflow
tet_tot  = NaN(size(pet));  % Simulated Actual Evapotranspiration
base_tot = NaN(size(pet));  % Simulated Base Flow
surf_tot = NaN(size(pet));  % Simulated Surface&Subsurface water flow


% PERFORMS SIMULATION USING THE SAC-SMA COUPLED WITH SNOW17
thres_zero  = 0.00001;      % Threshold to be considered as zero
parea       = 1 - adimp - pctim;
for i = 1:length(pet)
    
    % Precipitation adjusted by snow process (SNOW17)
	pr = snow_outflow(i);
 
      
    
    % COMPUTE ET LOSS FOR A TIME INTERVAL
    edmnd = pet(i);  % The amount of maximum ET given by input potential evapotranspiration
    
    % ET(1), ET from Upper zone tension water storage
    et1 = edmnd * uztwc/uztwm;
    red = edmnd - et1;  % residual ET demand
    uztwc = uztwc - et1;
   
    
    % ET(2), ET from upper zone free water storage
    et2 = 0;
    if uztwc <= 0    % in case et1 > uztws, no water in the upper tension water storage
        et1 = et1 + uztwc; % et1 = uztwc
        uztwc = 0;
        red = edmnd - et1;
        if uzfwc < red    % when upper zone free water content is less than residual ET
            et2 = uzfwc;  % all content at upper zone free water zone will be gone as ET
            uzfwc = 0;
            red = red - et2;
            if uztwc < thres_zero; uztwc = 0; end
            if uzfwc < thres_zero; uzfwc = 0; end
    
        else  % when upper zone free water content is larger than residual ET
            et2 = red;  % all residual ET will be gone as ET
            uzfwc = uzfwc - et2;
            red = 0;
   
        end
        
    else % all maximum et (et1) are consumed at uztwc, so no et from uzfwc (et2=0)
        % There's possibility that upper zone free water ratio exceeds upper zone tension water ratio
        % If so, free water is transferred to tension water storage
        if (uztwc / uztwm) < (uzfwc / uzfwm)
            uzrat = (uztwc + uzfwc) / (uztwm + uzfwm);
            uztwc = uztwm * uzrat;
            uzfwc = uzfwm * uzrat;
        end
        
        if uztwc < thres_zero; uztwc = 0; end
        if uzfwc < thres_zero; uzfwc = 0; end
        
    end
    
    
    % ET(3), ET from Lower zone tension water storage when residual ET > 0
    et3 = red * lztwc / (uztwm + lztwm); % according to this calculation, residual ET is always bigger than ET(3)
    lztwc = lztwc - et3; 
    if lztwc < 0 % et3 cannot exceed lztws
        et3 = et3 + lztwc;  % et3 = lztwc  
        lztwc = 0;
    end
    % ????? why is the red not updated? not updated red is used later for ET(5) calculation
    
  
    % Water resupply from Lower free water storages to Lower tension water storage
    saved = rserv * (lzfpm + lzfsm);
    ratlzt = lztwc / lztwm;
    ratlz = (lztwc + lzfpc + lzfsc - saved) / (lztwm + lzfpm + lzfsm - saved);
    if ratlzt < ratlz   % water is first taken from supplementary water storage for resupply
        del = (ratlz - ratlzt) * lztwm;
        lztwc = lztwc + del;  % Transfer water from lzfss to lztws
        lzfsc = lzfsc - del;
        if lzfsc < 0  % if tranfer exceeds lzfsc then remainder comes from lzfps
            lzfpc = lzfpc + lzfsc;
            lzfsc = 0;
        end
    end
    if lztwc < thres_zero; lztwc = 0; end
    
    
    % ET(5), ET from additional impervious (ADIMP) area
    et5 = et1 + (red + et2) * (adimc - et1 - uztwc) / (uztwm + lztwm);  % ????? no idea where this come from, I think there's a possibility that et5 can be negative values
    adimc = adimc - et5;
    if adimc < 0   % et5 cannot exceed adims
        et5 = et5 + adimc; % et5 = adimc
        adimc = 0;
    end
    et5 = et5 * adimp;
    
    
    % Time interval available moisture in excess of uztw requirements
    twx = pr + uztwc - uztwm;
    if twx < 0     % all moisture held in uztw- no excess
        uztwc = uztwc + pr;
        twx = 0;
    else    % moisture available in excess of uztw storage
        uztwc = uztwm;
    end
    % for now twx is excess rainfall after filling the uztwc
    
    adimc = adimc + pr - twx;
    
    % Compute Impervious Area Runoff
    roimp = pr * pctim; 
    
    
    % Initialize time interval sums
    sbf   = 0;  % Sum of total baseflow(from primary and supplemental storages)
    ssur  = 0;  % Sum of surface runoff
    sif   = 0;  % Sum of interflow
    sperc = 0;  % Time interval summation of percolation
    sdro  = 0;  % Sum of direct runoff from the additional impervious area
    
        
    % Determine computational time increments for the basic time interval
    ninc = floor(1.0 + 0.2*(uzfwc+twx));  % Number of time increments that interval is divided into for further soil-moisture accountng
    dinc = 1.0 / ninc;                    % Length of each increment in days
    pinc = twx / ninc;                    % Amount of available moisture for each increment
    
    % Compute free water depletion fractions for the time increment (basic depletions are for one day)
    duz   = 1 - (1 - uzk)^dinc;
    dlzp  = 1 - (1 - lzpk)^dinc;
    dlzs  = 1 - (1 - lzsk)^dinc;
    
    
    % Start incremental for-loop for the time interval
    for n = 1:ninc
        
        adsur = 0; % Amount of surface runoff. This will be updated.
        
        % Compute direct runoff from adimp area
        ratio = (adimc - uztwc) / lztwm;
        if ratio < 0; ratio = 0; end
        addro = pinc*(ratio^2); % Amount of direct runoff from the additional impervious area     
        
        % Compute baseflow and keep track of time interval sum.
        bf_p = lzfpc * dlzp; % Baseflow from free water primary storage
        lzfpc = lzfpc - bf_p;
        if lzfpc <= 0.0001
            bf_p = bf_p + lzfpc;
            lzfpc = 0;
        end
        sbf = sbf + bf_p;
        
        bf_s = lzfsc * dlzs;  % Baseflow from free water supplemental storage
        lzfsc = lzfsc - bf_s;
        if lzfsc <= 0.0001
            bf_s = bf_s + lzfsc;
            lzfsc = 0;
        end
        sbf = sbf + bf_s; % Total Baseflow from primary and supplemental storages
        
        % Compute PERCOLATION- if no water available then skip.
        if (pinc + uzfwc) <= 0.01
            
            uzfwc = uzfwc + pinc;
 
        else
            percm = lzfpm * dlzp + lzfsm * dlzs; % Limiting drainage rate from the combined saturated lower zone storages
            perc = percm * uzfwc / uzfwm;
            
            defr = 1.0 - (lztwc + lzfpc + lzfsc)/(lztwm + lzfpm + lzfsm); % DEFR is the lower zone moisture deficiency ratio
            
            if defr < 0; defr = 0; end
   
            perc = perc * (1.0 + zperc * (defr^rexp));
            
            % Note. . . percolation occurs from uzfws before pav is added
            
            if perc >= uzfwc     % Percolation rate exceeds uzfws
                perc = uzfwc;
            end
            
            uzfwc = uzfwc - perc;    % Percolation rate is less than uzfws.
            
            check = lztwc + lzfpc + lzfsc + perc - lztwm - lzfpm - lzfsm;
            
            if check > 0    % Check to see if percolation exceeds lower zone deficiency.
                perc = perc - check;
                uzfwc = uzfwc + check;
            end
            
            sperc = sperc + perc;  % SPERC is the time interval summation of PERC
            
            
            % Compute interflow and keep track of time interval sum. Note that PINC has not yet been added.
            del = uzfwc * duz; % The amount of interflow
            sif = sif + del;
            uzfwc = uzfwc - del;
            
            % Distribute percolated water into the lower zones. Tension water
            % must be filled first except for the PFREE area. PERCT is
            % percolation to tension water and PERCF is percolation going to
            % free water.
            perct = perc * (1.0 - pfree);  % Percolation going to the tension water storage
            if (perct + lztwc) <= lztwm
                lztwc = lztwc + perct;
                percf = 0.0; % Pecolation going to th lower zone free water storages
            else
                percf = lztwc + perct - lztwm;
                lztwc = lztwm;
            end
            
            % Distribute percolation in excess of tension requirements among the free water storages.
            percf = percf + (perc * pfree);
            if percf ~= 0
               
                hpl = lzfpm / (lzfpm + lzfsm); % Relative size of the primary storage as compared with total lower zone free water storages.
                
                % Relative fullness of each storage.
                ratlp = lzfpc / lzfpm;
                ratls = lzfsc / lzfsm;
                
                fracp = hpl * 2.0 * (1.0 - ratlp) / (1.0 - ratlp + 1.0 - ratls); % The fraction going to primary
                if fracp > 1.0; fracp = 1.0; end 
       
                percp = percf * fracp; % Amount of the excess percolation going to primary
                percs = percf - percp; % Amount of the excess percolation going to supplemental
                lzfsc = lzfsc + percs;
                if lzfsc > lzfsm
                    percs = percs - lzfsc + lzfsm;
                    lzfsc = lzfsm;
                end
                lzfpc = lzfpc + percf - percs;
                
                if lzfpc >= lzfpm  % Check to make sure lzfps does not exceed lzfpm
                    excess = lzfpc - lzfpm;
                    lztwc = lztwc + excess;
                    lzfpc = lzfpm;
                end
            end
            
            % Distribute PINC between uzfws and surface runoff
            if pinc ~= 0
                if (pinc + uzfwc) <= uzfwm   % check if pinc exceeds uzfwm
                    uzfwc = uzfwc + pinc;  % no surface runoff
                else
                    sur = pinc + uzfwc - uzfwm; % Surface runoff
                    uzfwc = uzfwm;
         
                    ssur = ssur + (sur * parea);
                    % ADSUR is the amount of surface runoff which comes from
                    % that portion of adimp which is not currently generating
                    % direct runoff. ADDRO/PINC is the fraction of adimp
                    % currently generating direct runoff.
                    adsur = sur * (1.0 - addro / pinc); 
                    ssur = ssur + adsur * adimp;
                end
            end
            
         
        end
        
        adimc = adimc + pinc - addro - adsur;
        if adimc > (uztwm + lztwm)
            addro = addro + adimc - (uztwm + lztwm);
            adimc = uztwm + lztwm;
        end
        
        sdro  = sdro + (addro * adimp); % Direct runoff from the additional impervious area
        
        if adimc < thres_zero;  adimc = 0; end
        
    end  % END of incremental for loop
    
    
    % Compute sums and adjust runoff amounts by the area over which they are generated.
    
    % EUSED is the ET from PAREA which is 1.0 - adimp - pctim
    eused = et1 + et2 + et3;
    sif = sif * parea;
    
    % Separate channel component of baseflow from the non-channel component
    tbf = sbf * parea;  % TBF is the total baseflow
    bfcc = tbf * (1.0 / (1.0 + side));    % BFCC is baseflow, channel component  
        
    % Ground flow and Surface flow
    base = bfcc;                      % Baseflow and Interflow are considered as Ground inflow to the channel
    surf = roimp + sdro + ssur + sif; % Surface flow consists of Direct runoff and Surface inflow to the channel
    
    % ET(4)- ET from riparian vegetation.
    et4 = (edmnd - eused) * riva; % no effect if riva is set to zero
    
    
    % Compute total evapotransporation - TET
    eused = eused * parea;
    tet = eused + et4 + et5;

    
    % Check that adims >= uztws
    if adimc < uztwc 
        adimc = uztwc;
    end
    
    
    % Total inflow to channel for a timestep 
    tot_outflow = surf + base;
    
    tot_outflow = tot_outflow - et4;
    if tot_outflow < 0
        tot_outflow = 0;
		surf = 0;
		base = 0;
    end   

    
    % OUTPUTS FOR A TIME STEP
    % Main output: Watershed out flow
    simflow(i) = tot_outflow;
    % Optional outputs 
    tet_tot(i)   = tet;
    surf_tot(i)  = surf;
    base_tot(i)  = base;
    uztwc_tot(i) = uztwc;
    uzfwc_tot(i) = uzfwc;
    lztwc_tot(i) = lztwc;
    lzfpc_tot(i) = lzfpc;
    lzfsc_tot(i) = lzfsc;
    adimc_tot(i) = adimc;
     
end