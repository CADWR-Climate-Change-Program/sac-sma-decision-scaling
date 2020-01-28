clear; clc
tic;
%%%%%%%%%%%%%%%%%%%
% DIRECTORIES
SACSMADSdir = 'C:/Users/warnold/Box/CC-Data-Cloud/data_SAC-SMA-DS/repo/';
outputLocation = 'C:/Users/warnold/Box/CC-Data-Cloud/data_SAC-SMA-DS/Paleo_SJR/Paleo_Jitter_MonotonicTemp/';
climateData = 'C:/Users/warnold/Box/CC-Data-Cloud/data_paleoCoupling/processed/paleo_sjr_analogue_sim_ccvs/';
climatePrefix = 'data_';

% RUN CONFIG
run_type = 'watershed_specific'; % choice of 'watershed_specific' or 'callite'
outputDailyAvgHydroPlot = false;
outputDailyHydroSequence = true;
outputDailyFlowsOnly = false;
set_to_run = '11ObsInflows'; %choice of '9Unimpaired','11ObsInflows','12RimInflows'
siteStart = 5;
siteEnd = 5;

% Leap year output
leap = true; % only for daily watershed specific outputs
startYr = 894;
endYr = 2011;

% Temperature Perturbations
distributedTemp = false;
T_rangeStart = 1;
T_rangeEnd = 1;
T_change = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
T_change_label = {'0_0','0_5','1_0','1_5','2_0','2_5','3_0','3_5','4_0'};
%T_change = [0, 1.0, 2.0, 3.0, 4.0];
%T_change_label = {'0_0','1_0','2_0','3_0','4_0'};
T_distributed_sensitivity = [0.116,0.116,0.265,0.265,0.265,0.1915,0.1915,0.1915,0.0435,0.0435,0.0435,0.116];

% Precipitation Perturbations
P_rangeStart = 2;
P_rangeEnd = 2;
P_change = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3];
P_change_label = {'0_7','0_8','0_9','1_0','1_1','1_2','1_3'};

% SAC-SMA-DS Initial State
inistate = [0 0 5 5 5 0];

%%%%%%%%%%%%%%%%%%%%%
% FLOW SET SELECTION
sets = {'9Unimpaired','11ObsInflows','12RimInflows'};
watershedNames = {{'BearRiver','CacheCreek','CalaverasRiver','ChowchillaRiver','ConsumnesRiver','FresnoRiver','MokelumneRiver','PutahCreek','StonyCreek'}
    {'AMF','BLB','BND','FTO','MRC','SHA','SJF','SNS','TLG','TNL','YRS'}
    {'FOL_I','LK_MC','N_MEL','MILLE','PRD_C','N_HOG','OROVI','DPR_I','SHAST','TRINI','SMART','WKYTN'}};
watershedAreaSqMeter = {{730771760,2455821437,927537363,685393248,1384751552,717246675,1408736800,1477609964,1926664012};
    {4824193789,1926663715,25843435254,9431450110,2727572600,19346624956,4340985255,2539365534,3979829811,1861176298,2875096534}
    {4824193789,2727572601,2539365535,4340985255,1408736547,927537198,9431450110,3979829812,19346624957,1751131865,2875096534,508807202}};
watershedNameDict = containers.Map(sets,watershedNames);
watershedsToRun = watershedNameDict(set_to_run);
watershedAreaDict = containers.Map(sets,watershedAreaSqMeter);
watershedsToRunArea = watershedAreaDict(set_to_run);

% Create directories
if strcmp(run_type,'watershed_specific')
    mkdir (outputLocation,['Daily/',set_to_run])
end
if outputDailyAvgHydroPlot == true
    mkdir (outputLocation,['Daily/',set_to_run,'/Hydrologic-Plots/'])
end

fprintf(['Processing ', set_to_run,' for ',run_type, '\nDistributed Temp =', num2str(distributedTemp),'\nOutput Location: ',regexprep(outputLocation,"\","\\\"),'\n'])

% CALLITE HEADER FOR 12 RIM INFLOWS
firstLine = '31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400';
secondLine = ['/CALLITE/I_FOLSM/FLOW-INFLOW//1MON/2020D09E/,',...
    '/CALLITE/I_MCLRE/FLOW-INFLOW//1MON/2020D09E/,'...
    '/CALLITE/I_MELON/FLOW-INFLOW//1MON/2020D09E/,'...
    '/CALLITE/I_MLRTN/FLOW-INFLOW//1MON/2020D09E/,'...
    '/CALLITE/I_MOKELUMNE/FLOW-INFLOW//1MON/2020D09E/,'...
    '/CALLITE/I_NHGAN/FLOW-INFLOW//1MON/2020D09E/,'...
    '/CALLITE/I_OROVL/FLOW-INFLOW//1MON/2020D09E/,'...
    '/CALLITE/I_PEDRO/FLOW-INFLOW//1MON/2020D09E/,'...
    '/CALLITE/I_SHSTA/FLOW-INFLOW//1MON/2020D09E/,'...
    '/CALLITE/I_TRNTY/FLOW-INFLOW//1MON/2020D09E/,'...
    '/CALLITE/I_YUBA/FLOW-INFLOW//1MON/2020D09E/,'...
    '/CALLITE/I_WKYTN/FLOW-INFLOW//1MON/2020D09E/'];
thirdLine = 'CFS,TAF,TAF,CFS,CFS,TAF,CFS,TAF,CFS,CFS,CFS,CFS';

% SET SAC-SMA-DS WORKING PATH FOR MODULE & PARAMETER REFERENCE
p = genpath(SACSMADSdir);
addpath(p)

% LABEL FOR WARMING TYPE
if distributedTemp == true
    Warming_label = 'Distributed Temp';
else
    Warming_label = 'Montonic Temp';
end

% DATE VECTORS
genyr = endYr - startYr;
strdate = [startYr 10 1];
enddate = [endYr 9 30];
datemat = datevec(datenum(strdate):datenum(enddate));
datemat_mon = [];
for i=startYr:endYr
    datemat_mon = [datemat_mon; repmat(i,12,1),(1:12)'];
end
datemat_mon = datemat_mon(10:end-3,:); % remove first 9 months and last 3 months (to create WY 895-2011)

if leap == true
    ymdLeap = datetime(datemat(:,[1 2 3]));
    yearLeap = datemat(:,1);
    monthLeap = datemat(:,2);
    domLeap = datemat(:,3);
    doyLeap = day(datetime(datemat(:,[1 2 3])),'dayofyear');
    dateCal = datenum(datemat(:,[1 2 3]));
    datematLeap=table(ymdLeap,dateCal,yearLeap,monthLeap,domLeap,doyLeap);
end

datemat(datemat(:,2) == 2 & datemat(:,3) == 29,:) = []; % remove leap days
ymd = datetime(datemat(:,[1 2 3]));
year = datemat(:,1);
month = datemat(:,2);
dom = datemat(:,3);
doy = day(datetime(datemat(:,[1 2 3])),'dayofyear');
dateCal = datenum(datemat(:,[1 2 3]));

% change year to water year
ond = datemat_mon(:,2) == 10 | datemat_mon(:,2) == 11 | datemat_mon(:,2) == 12;
datemat_mon(ond,1) = datemat_mon(ond,1) + 1;
            
% step count for temp and precip change
% tp_change = 1;

for tc = T_rangeStart:T_rangeEnd
    
    for pc = P_rangeStart:P_rangeEnd
   
        mon_datemat_day = eomday(datemat_mon(:,1),datemat_mon(:,2));
 
        mon_simflow_tot = [];

        for site = siteStart:siteEnd

            hru_par = [];
            
            gridinfo = load(['parameters/GridInfo/',set_to_run,'/GridInfo_',watershedsToRun{site},'.txt']);
            grid_lat = gridinfo(:,1);
            grid_lon = gridinfo(:,2);
            grid_area = gridinfo(:,3);
            grid_elev = gridinfo(:,4);
            grid_flen = gridinfo(:,5);
            numtotgrid = length(grid_lat);
            
            fid = fopen(char(strcat('parameters/calibration/',set_to_run,'/mpiresult_SACSMA_DS_',watershedsToRun(site))));
            
            while ~feof(fid)
                str = fgets(fid);
                if length(str)>=61 && strcmp(str(1:61),'                                         GA OPTIMUM PARAMETER')
                    fgets(fid);
                    fgets(fid);
                    for jj = 1:numtotgrid
                        str = fgets(fid);
                        hru_par = [hru_par;textscan(str,['%s %s '...
                        '%f %f %f %f %f %f %f %f %f %f '...
                        '%f %f %f %f %f %f %f %f %f %f '...
                        '%f %f %f %f %f %f %f %f %f %f %f '])];
                    end
                    break;
                end
            end
            fclose(fid);

            isOutlet = zeros(length(grid_lat),1);
            [~,minind] = min(grid_flen);
            isOutlet(minind) = 1;
            
            PR = zeros(genyr*365,length(gridinfo));
            TAS = zeros(genyr*365,length(gridinfo));
            PET = zeros(genyr*365,length(gridinfo));
            TOTAL_ET = zeros(genyr*365,length(gridinfo));
            UZTWC = zeros(genyr*365,length(gridinfo));
            UZFWC = zeros(genyr*365,length(gridinfo));
            LZTWC = zeros(genyr*365,length(gridinfo));
            LZFSC = zeros(genyr*365,length(gridinfo));
            LZFPC = zeros(genyr*365,length(gridinfo));
            ADIMC = zeros(genyr*365,length(gridinfo));
            SNOW_SWE = zeros(genyr*365,length(gridinfo));
            SNOW_OUTFLOW = zeros(genyr*365,length(gridinfo));
            SNOW_MELT = zeros(genyr*365,length(gridinfo));
            BASEFLOW = zeros(genyr*365,length(gridinfo));
            FLOW = zeros(genyr*365,length(gridinfo));

            syr = strdate(1);
            eyr = enddate(1);

            disp(['Running Temp (C) +',num2str(T_change(tc)),' and Prcp *',num2str(P_change(pc)),' for site: ',watershedsToRun{site}, ' (',num2str(length(gridinfo)),' grids)'])
            
            for i = 1:length(gridinfo)
                
                % disp (['Processing ', num2str(grid_lat(i),'%2.4f'),'_',num2str(grid_lon(i),'%2.4f')])
                % load climate data
                climate = load([climateData,climatePrefix,num2str(grid_lat(i),'%2.4f'),'_',num2str(grid_lon(i),'%2.4f')]);

                % perturb climate data
                pr = climate(end-size(datemat,1)+1:end,4) * P_change(pc);
                
                if distributedTemp == true
                    tas = climate(end-size(datemat,1)+1:end,5);
                    for nm = 1:12
                        tas(datemat(:,2)==nm) = tas(datemat(:,2)==nm)+ ( T_change(tc) * T_distributed_sensitivity(nm)/mean(T_distributed_sensitivity));
                    end
                else
                    tas = climate(end-size(datemat,1)+1:end,5) + T_change(tc);
                end
                
                selind = 0;
                for j=1:numtotgrid
                    if strcmp(hru_par{j,1}{1},num2str(grid_lat(i),'%10.6f')) && strcmp(hru_par{j,2}{1},num2str(grid_lon(i),'%10.6f'))
                        selind = j;
                        %                 break;
                    end
                end

                Coeff = hru_par{selind,19};
                SnowPar = [hru_par{selind,20},hru_par{selind,21},hru_par{selind,22},hru_par{selind,23},hru_par{selind,24},hru_par{selind,25},hru_par{selind,26},hru_par{selind,27},hru_par{selind,28},hru_par{selind,29}];
                SMA_Par = [hru_par{selind,3},hru_par{selind,4},hru_par{selind,5},hru_par{selind,6},hru_par{selind,7},hru_par{selind,8},hru_par{selind,9},hru_par{selind,10}...
                    ,hru_par{selind,11},hru_par{selind,12},hru_par{selind,13},hru_par{selind,14},hru_par{selind,15},hru_par{selind,16},hru_par{selind,17},hru_par{selind,18}];
                route_par = [hru_par{selind,30},hru_par{selind,31},hru_par{selind,32},hru_par{selind,33}];

                pet = nan(genyr*365,1);
                n = 1;
                for j = 1:genyr
                    pet(n:n+364) = pet_hamon([startYr 10 1], [startYr+1 9 30], tas(n:n+364), grid_lat(i), Coeff);
                    n = n +365;
                end

                snow_outflow = nan(genyr*365,1);
                snow_melt = nan(genyr*365,1);
                snow_swe = nan(genyr*365,1);
                snow_inistate = [0 0 0 0];
                n = 1;
                for j = 1:genyr
                    [snow_outflow(n:n+364), snow_melt(n:n+364), snow_swe(n:n+364), snow_inistate_new] = snow_snow17([startYr 10 1], [startYr+1 9 30], pr(n:n+364), tas(n:n+364), grid_elev(i), SnowPar, snow_inistate);
                    snow_inistate = snow_inistate_new;
                    n = n +365;
                end

                [inflow_direct, inflow_base, total_et, uztwc, uzfwc, lztwc, lzfpc, lzfsc, adimc] = sma_sacsma(pet,snow_outflow,SMA_Par,inistate);
                [runoff, baseflow] = routing_lohmann(inflow_direct, inflow_base, grid_flen(i), route_par, isOutlet(i));
                hru_flow = runoff;
                
                TAS(:,i) = tas(1:end);
                PET(:,i) = pet * grid_area(i)/sum(grid_area);
                PR(:,i) = pr(1:end) * grid_area(i)/sum(grid_area);
                TOTAL_ET(:,i) = total_et * grid_area(i)/sum(grid_area);
                UZTWC(:,i) = uztwc * grid_area(i)/sum(grid_area);
                UZFWC(:,i) = uzfwc * grid_area(i)/sum(grid_area);
                LZTWC(:,i) = lztwc * grid_area(i)/sum(grid_area);
                LZFSC(:,i) = lzfpc * grid_area(i)/sum(grid_area);
                LZFPC(:,i) = lzfsc * grid_area(i)/sum(grid_area);
                ADIMC(:,i) = adimc * grid_area(i)/sum(grid_area);
                SNOW_OUTFLOW(:,i) = snow_outflow * grid_area(i)/sum(grid_area);
                SNOW_MELT(:,i) = snow_melt * grid_area(i)/sum(grid_area);
                SNOW_SWE(:,i) = snow_swe * grid_area(i)/sum(grid_area);
                BASEFLOW(:,i) = baseflow * grid_area(i)/sum(grid_area);
                FLOW(:,i) = hru_flow * grid_area(i)/sum(grid_area);

            end

            % April 1 Snow Water Equivalent
            %april1st = (datemat(:,2)==4 & datemat(:,3)==1);
            %swe_1st  = SNOW_SWE(april1st==1,:);
            %swe_1st_ws = mean(swe_1st,2);
            
            % Monthly PET
            %pet_ws = mean(PET,2);
            %mon_pet_ws = grpstats(pet_ws,datemat(:,2),'mean');
            %mon_pet_ws = round(mon_pet_ws,3);
            
            % monthly streamflow
            mon_simflow = grpstats(sum(FLOW,2),{datemat(:,1),datemat(:,2)},'sum');
            if strcmp(set_to_run,'12RimInflows')
                mon_simflow = mon_simflow./mon_datemat_day/304.8/24/3600*watershedsToRunArea{site}*10.7639; % convert to CFS
                if site == 2 || site == 3 || site == 6 || site == 8 % convert these sites to TAF
                    mon_simflow = mon_simflow * 0.06;
                end
            end
            mon_simflow_tot = [mon_simflow_tot, mon_simflow]; 
            
            %watershed specific daily outputs
            if strcmp(run_type,'watershed_specific')
                %DAILY HYDROLOGY 
                flow = round(sum(FLOW,2),4); %streamflow values in mm
                % mm to af:: / 1000 * watershedsToRunArea{site} * 0.000810714; 
                flowVol = sum(FLOW,2) * watershedsToRunArea{site} / 304.8 / 24 / 3600 * 10.7639; % mm to cfs
                baseflow = round(sum(BASEFLOW,2),4); %baseflow values in mm
                pr = round(sum(PR,2),4); %precip values in mm
                tas = round(mean(TAS,2),4); %avg temp values in C
                swe = round(sum(SNOW_SWE,2),4); %SWE values in mm
                pet = round(sum(PET,2),4); %Potential ET values in mm
                et = round(sum(TOTAL_ET,2),4); %Actual ET values in mm
                snowRainFlow = round(sum(SNOW_OUTFLOW,2),4); %Gets daily snow outflow (rain+snow) values
                snowMelt = round(sum(SNOW_MELT,2),4); %Gets daily snow melt values
                uztwc = round(sum(UZTWC,2),4); %upper zone tension water values in mm
                uzfwc = round(sum(UZFWC,2),4); %Upper zone free water in mm
                lztwc = round(sum(LZTWC,2),4); %Lower zone tension water in mm
                lzfsc = round(sum(LZFSC,2),4); %Lower zone supplementary free water in mm
                lzfpc = round(sum(LZFPC,2),4); %Upper zone primary free water in mm
                adimc = round(sum(ADIMC,2),4); %Additional impervious area water in mm
                
                flowVolTable = table(dateCal,flowVol);
                
                wb = table(ymd,dateCal,year,month,dom,doy,tas,pet,et,pr,snowRainFlow,snowMelt,swe,baseflow,flow,uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc);
                
                rows = wb.doy~=366;
                wbd = grpstats(wb(rows,:),'doy','mean','DataVars',{'et','pr','tas','snowRainFlow','snowMelt','swe','flow'});
                
                if outputDailyAvgHydroPlot == true
                    % PLOT AND SAVE FIGURE
                    subplot(5,1,[1,2]);
                    hold on
                    yyaxis left
                    plot(wbd.mean_et,'DisplayName','ET')
                    plot(wbd.mean_pr,'DisplayName','Pr')
                    plot(wbd.mean_snowMelt,'DisplayName','SnowMelt')
                    plot(wbd.mean_snowRainFlow,'DisplayName','SnowOutflow')
                    ylabel('mm')
                    ylim([0 10])
                    yyaxis right
                    plot(wbd.mean_tas,'DisplayName','Temp(c)')
                    ylabel('C')
                    ylim([0 25])
                    hold off
                    title({strcat('SAC-SMA-DS Hydrologic Outputs:','  ', watershedsToRun{site}, ' (',...
                        num2str(round(watershedsToRunArea{site}*1e-6,0),'%.0f'),' km2)');...
                        strcat(Warming_label,' |  +',num2str(T_change(tc)),'C |  x',num2str(P_change(pc)),' precip')});
                    legend('show')
                    legend('boxoff')
                    hold on
                    subplot(5,1,[3,4,5]);
                    yyaxis left
                    plot(wbd.mean_flow,'DisplayName','Flow')
                    ylabel('mm')
                    ylim([0 10])
                    yyaxis right
                    plot(wbd.mean_swe,'DisplayName','SWE')
                    ylabel('mm')
                    ylim([0 300])
                    xlabel('Day of Year')
                    hold off
                    legend('show')
                    legend('boxoff')
                    print(strcat(outputLocation,'Daily/',set_to_run,'/Hydrologic-Plots',...
                        '/dailyAvg_',watershedsToRun{site},'_',P_change_label{pc},'DP',T_change_label{tc},'DT'),'-dpdf','-fillpage')
                    close
                end
                
                if leap == true
                    
                    wb = outerjoin(wb,datematLeap,'Keys','dateCal');
                    wb = [wb.tas wb.pet wb.et wb.pr, wb.snowRainFlow...
                        wb.snowMelt wb.swe wb.baseflow wb.flow];
                    wb = fillmissing(wb,'movmedian',[1,1]);
                    wb = table(ymdLeap, yearLeap, monthLeap, domLeap, doyLeap, wb(:,1),...
                        wb(:,2),wb(:,3),wb(:,4),wb(:,5),wb(:,6),wb(:,7),wb(:,8),wb(:,9),...
                        'VariableNames',{'ymd','year','month','dom','doy',...
                        'tas','pet','et','pr','snowRainFlow','snowMelt','swe','baseflow','flow'});
                    
                    flowVolTable = outerjoin(flowVolTable,datematLeap,'Keys','dateCal');
                    flowVolTable = [flowVolTable.flowVol];
                    flowVolTable = fillmissing(flowVolTable,'movmedian',[1,1]);
                    flowVolTable = table(yearLeap, monthLeap, domLeap,flowVolTable(:,1),...
                        'VariableNames',{'Year','Month','Day','Flow_cfs'});
                    
                    % flowVolTableAnnual = grpstats(flowVolTable(:,:),'Year','sum','DataVars',{'Flow_AF'});
                end
                if outputDailyHydroSequence == true
                    %Note (2191:end) removes first 6 yrs (WY895- WY900) of model runup
                    writetable(wb(2192:end,:),strcat(outputLocation,'Daily/',set_to_run,...
                        '/SACSMA_','DT',num2str(T_change(tc)),'_DP',num2str(P_change(pc)),'_',watershedsToRun{site},'_full.txt'),'Delimiter',' ');
                end
                if outputDailyFlowsOnly == true
                    %Note (2191:end) removes first 6 yrs (WY895- WY900) of model runup
                    writetable(flowVolTable(2192:end,:),strcat(outputLocation,'Daily/',set_to_run,...
                        '/SACSMA_','DT',num2str(T_change(tc)),'_DP',num2str(P_change(pc)),'_',watershedsToRun{site},'.txt'),'Delimiter',' ');
                end
            end
        end
        
        % CalLite outputs
        if strcmp(run_type,'callite')
            
            

            mkdir (outputLocation,['Monthly/',strcat(P_change_label{pc},'DP',T_change_label{tc},'DT')])
            
            if strcmp(set_to_run,'9Unimpaired')
                
                fid = fopen([outputLocation,'Monthly/',P_change_label{pc},'DP',T_change_label{tc},'DT/','SACSMA_Inflow_CalLite_9.txt'],'w');
                fprintf(fid,'%s\n','Year, Month, BearRiver, CacheCreek, CalaverasRiver, ChowchillaRiver, CosumnesRiver, FresnoRiver, MokelumneRiver, PutahCreek, StonyCreek');
                fprintf(fid,'%4.0f, %2.0f, %14.6f, %14.6f, %14.6f, %14.6f, %14.6f, %14.6f, %14.6f, %14.6f, %14.6f\n',[datemat_mon(61:end,:) mon_simflow_tot(61:end,:)]');% remove first 5 yrs (WY895- WY899)
                fclose(fid);

            elseif strcmp(set_to_run,'11ObsInflows')
                
                fid = fopen([outputLocation,'Monthly/',P_change_label{pc},'DP',T_change_label{tc},'DT/','SACSMA_Inflow_CalLite_11.txt'],'w');
                fprintf(fid,'%s\n','Year, Month, AMF, BLB, BND, FTO, MRC, SHA, SJF, SNS, TLG, TNL, YRS');
                fprintf(fid,'%4.0f, %2.0f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f\n',[datemat_mon(61:end,:) mon_simflow_tot(61:end,:)]'); % remove first 5 yrs (WY895- WY899)
                fclose(fid);
                
            elseif strcmp(set_to_run,'12RimInflows')
                
                fid = fopen([outputLocation,'Monthly/',P_change_label{pc},'DP',T_change_label{tc},'DT/','SACSMA_Inflow_CalLite_12.txt'],'w');
                fprintf(fid,'%s\n',firstLine);
                fprintf(fid,'%s\n',secondLine);
                fprintf(fid,'%s\n',thirdLine);
                fprintf(fid,'%1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f\n',mon_simflow_tot(61:end,:)');% remove first 5 yrs (WY895- WY899)
                fclose(fid);
                
            end
        end

    %swe_1st_ws_tp(:,tp_change) = swe_1st_ws;
    %tp_change = tp_change + 1;
    
    end
    
    %mon_pet_ws_t(:,tc) = mon_pet_ws;

end
toc;


