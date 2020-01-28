

clear; clc
tic;
% parpool
foldername = ls('C:\Users\Wi\Dropbox\CalLite\file transfer\WGEN_SIMULATIONS_MAY2017\T*');
trial_all = [];
for i = 1:size(foldername,1)
    trial_all = [trial_all; str2double(foldername(i,7:end))];
end
trial_all = sort(trial_all);

site = {'FOL_I','LK_MC','N_MEL','MILLE','PRD_C','N_HOG','OROVI','DPR_I','SHAST','TRINI','SMART','WKYTN'};
firstlow = '31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400,31OCT1921 2400';
secondlow = ['/CALLITE/I_FOLSM/FLOW-INFLOW//1MON/2020D09E/,',...
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
thirdlow = 'CFS,TAF,TAF,CFS,CFS,TAF,CFS,TAF,CFS,CFS,CFS,CFS';
addpath('E:\CALIFORNIA\Kevin_dwr\BatchScript\module')

% genyr = 87;
genyr = 55;

strdate = [1921 10 1];
% enddate = [2008 9 30];
enddate = [1976 9 30];
datemat = datevec(datenum(strdate):datenum(enddate));
datemat(datemat(:,2) == 2 & datemat(:,3) == 29,:) = [];

inistate = [0 0 5 5 5 0];

[~,~,areainfo] = xlsread('E:\CALIFORNIA\SubBasin\CalLiteInput_CAsubbasin.xlsx','Area');

T_change = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
P_change = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3];
T_change_label = {'0_0','0_5','1_0','1_5','2_0','2_5','3_0','3_5','4_0'};
P_change_label = {'0_7','0_8','0_9','1_0','1_1','1_2','1_3'};
for tc = 1%1:9
    
    for pc = 4%1:7
        
        mon_datemat = [];
        for i = strdate(1):enddate(1)
            mon_datemat = [mon_datemat;[repmat(i,12,1) (1:12)']];
        end
        
        mon_datemat_day = zeros(size(mon_datemat,1),1);
        for i = 1:size(mon_datemat,1)
            mon_datemat_day(i) = eomday(1921,mon_datemat(i,2));
        end
        mon_datemat_day = mon_datemat_day(10:end-3);
        
        h = waitbar(0,'Please Wait...');
        for itrial = 1:length(trial_all)

            mon_simflow_tot = [];
            for isite = 1:length(site)
                
                
                file = ls(['E:\CALIFORNIA\Kevin_dwr\Model Parameter File\calibration\12RimInflows\mpiresult_SACSMA_DS_',site{isite},'*']);
                
                areaind = 0;
                for i = 1:size(areainfo,1)
                    if strcmp(site{isite},areainfo(i,1))
                        areaind = i;
                    end
                end
                selarea = areainfo{areaind,2};
                
                
                gridinfo = load(['E:\CALIFORNIA\Kevin_dwr\Model Parameter File\GridInfo\12RimInflows\GridInfo_',site{isite},'.txt']);
                grid_lat = gridinfo(:,1);
                grid_lon = gridinfo(:,2);
                grid_area = gridinfo(:,3);
				grid_elev = gridinfo(:,4);
                grid_flen = gridinfo(:,5);
                
                
                numtotgrid = length(grid_lat);
                
                fid = fopen(['E:\CALIFORNIA\Kevin_dwr\Model Parameter File\calibration\12RimInflows\',deblank(file)]);
                
                hru_par = [];
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
                
                
                FLOW = zeros(genyr*365,length(gridinfo));
                
                syr = strdate(1);
                eyr = enddate(1);
                
                
                for i = 1:length(gridinfo)

                    climate = load(['C:\Users\Wi\Dropbox\CalLite\file transfer\WGEN_SIMULATIONS_MAY2017\TRIAL_',num2str(trial_all(itrial)),'\TRIAL_',num2str(trial_all(itrial)),'_D.P.MEAN_1_D.T.MEAN_1\data_',num2str(grid_lat(i),'%2.4f'),'_',num2str(grid_lon(i),'%2.4f')]);    
                    pr = climate(1:size(datemat,1),4)*P_change(pc);
                    tas = climate(1:size(datemat,1),5)+0.33+T_change(tc);

                    
                    selind = 0;
                    for j=1:numtotgrid
                        if strcmp(hru_par{j,1}{1},num2str(grid_lat(i),'%2.6f')) && strcmp(hru_par{j,2}{1},num2str(grid_lon(i),'%2.6f'))
                            selind = j;
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
                        pet(n:n+364) = pet_hamon([1921 10 1], [1922 9 30], tas(n:n+364), grid_lat(i), Coeff);
                        n = n +365;
                    end
					
					snow_outflow = nan(genyr*365,1);
					snow_melt = nan(genyr*365,1);
					snow_swe = nan(genyr*365,1);
                    snow_inistate = [0 0 0 0];
					n = 1;
                    for j = 1:genyr
                        [snow_outflow(n:n+364), snow_melt(n:n+364), snow_swe(n:n+364), snow_inistate_new] = snow_snow17([1921 10 1], [1922 9 30], pr(n:n+364), tas(n:n+364), grid_elev(i), SnowPar, snow_inistate);
                        snow_inistate = snow_inistate_new;
                        n = n +365;
                    end
                   
                    [inflow_direct, inflow_base] = sma_sacsma(pet,snow_outflow,SMA_Par, inistate);
                    [runoff, baseflow] = routing_lohmann(inflow_direct, inflow_base, grid_flen(i), route_par, isOutlet(i));
                    hru_flow = runoff;
                    FLOW(:,i) = hru_flow * grid_area(i)/sum(grid_area);
                    
                end
                
                mon_simflow = grpstats(sum(FLOW,2),{datemat(:,1),datemat(:,2)},'sum');              
                mon_simflow = mon_simflow./mon_datemat_day/304.8/24/3600*selarea;
                if isite == 2 || isite == 3 || isite == 6 || isite == 8
                    mon_simflow = mon_simflow * 0.06;
                end
                
                mon_simflow_tot = [mon_simflow_tot, mon_simflow(61:end)];
                                
            end
            
            fid = fopen(['C:\Users\Wi\Dropbox\CalLite\file transfer\12RimInflows\TRIAL_',num2str(trial_all(itrial)),'_SACSMA_Inflow_CalLite.txt'],'w');
            fprintf(fid,'%s\n',firstlow);
            fprintf(fid,'%s\n',secondlow);
            fprintf(fid,'%s\n',thirdlow);
            fprintf(fid,'%1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f\n',mon_simflow_tot');
            fclose(fid);
            
            waitbar(itrial/length(trial_all),h,[num2str(itrial),'/',num2str(length(trial_all))]);

        end
        close(h)
        
    end
    
end


toc;



%%


site = {'AMF','BLB','BND','FTO','MRC','SHA','SJF','SNS','TLG','TNL','YRS'};

% genyr = 87;
genyr = 55;

strdate = [1921 10 1];
% enddate = [2008 9 30];
enddate = [1976 9 30];
datemat = datevec(datenum(strdate):datenum(enddate));


datemat_mon = [];
for i=1927:1976%2008
    datemat_mon = [datemat_mon; repmat(i,12,1),[10,11,12,(1:9)]'];
end

% datemat2= datemat;

datemat(datemat(:,2) == 2 & datemat(:,3) == 29,:) = [];

for tc = 1%1:9
    
    for pc = 4%1:7
        
%        if ~(tc == 1 && pc == 4)
        
        mon_datemat = [];
        for i = strdate(1):enddate(1)
            mon_datemat = [mon_datemat;[repmat(i,12,1) (1:12)']];
        end
        
        mon_datemat_day = zeros(size(mon_datemat,1),1);
        for i = 1:size(mon_datemat,1)
            mon_datemat_day(i) = eomday(1921,mon_datemat(i,2));
        end
        mon_datemat_day = mon_datemat_day(10:end-3);
        
        h = waitbar(0,'Please Wait...');
        for itrial = 1:length(trial_all)
%             tic;
            mon_simflow_tot = [];
            for isite = 1:length(site)
                
                
                file = ls(['E:\CALIFORNIA\Kevin_dwr\Model Parameter File\calibration\11ObsInflows\mpiresult_SACSMA_DS_',site{isite},'*']);
                
                
                gridinfo = load(['E:\CALIFORNIA\Kevin_dwr\Model Parameter File\GridInfo\11ObsInflows\GridInfo_',site{isite},'.txt']);
                grid_lat = gridinfo(:,1);
                grid_lon = gridinfo(:,2);
                grid_area = gridinfo(:,3);
				grid_elev = gridinfo(:,4);
                grid_flen = gridinfo(:,5);
                
                
                numtotgrid = length(grid_lat);
                
                fid = fopen(['E:\CALIFORNIA\Kevin_dwr\Model Parameter File\calibration\11ObsInflows\',deblank(file)]);
                
                hru_par = [];
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
                
                
                FLOW = zeros(genyr*365,length(gridinfo));
                
                syr = strdate(1);
                eyr = enddate(1);
                
                
                for i = 1:length(gridinfo)
                    %for i = 1:length(gridinfo)
                    climate = load(['C:\Users\Wi\Dropbox\CalLite\file transfer\WGEN_SIMULATIONS_MAY2017\TRIAL_',num2str(trial_all(itrial)),'\TRIAL_',num2str(trial_all(itrial)),'_D.P.MEAN_1_D.T.MEAN_1\data_',num2str(grid_lat(i),'%2.4f'),'_',num2str(grid_lon(i),'%2.4f')]);
                    pr = climate(1:size(datemat,1),4)*P_change(pc);
                    tas = climate(1:size(datemat,1),5)+0.33+T_change(tc);
                    
                    selind = 0;
                    for j=1:numtotgrid
                        if strcmp(hru_par{j,1}{1},num2str(grid_lat(i),'%2.6f')) && strcmp(hru_par{j,2}{1},num2str(grid_lon(i),'%2.6f'))
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
                        pet(n:n+364) = pet_hamon([1921 10 1], [1922 9 30], tas(n:n+364), grid_lat(i), Coeff);
                        n = n +365;
                    end
					
					snow_outflow = nan(genyr*365,1);
					snow_melt = nan(genyr*365,1);
					snow_swe = nan(genyr*365,1);
                    snow_inistate = [0 0 0 0];
					n = 1;
                    for j = 1:genyr
                        [snow_outflow(n:n+364), snow_melt(n:n+364), snow_swe(n:n+364), snow_inistate_new] = snow_snow17([1921 10 1], [1922 9 30], pr(n:n+364), tas(n:n+364), grid_elev(i), SnowPar, snow_inistate);
                        snow_inistate = snow_inistate_new;
                        n = n +365;
                    end
                   
                    [inflow_direct, inflow_base] = sma_sacsma(pet,snow_outflow,SMA_Par, inistate);
                    [runoff, baseflow] = routing_lohmann(inflow_direct, inflow_base, grid_flen(i), route_par, isOutlet(i));
                    hru_flow = runoff;
                    FLOW(:,i) = hru_flow * grid_area(i)/sum(grid_area);
                    
                end
                
                mon_simflow = grpstats(sum(FLOW,2),{datemat(:,1),datemat(:,2)},'sum');
                
                mon_simflow_tot = [mon_simflow_tot, mon_simflow(61:end)];
                
                
            end
            
%             fid = fopen(['H:\WG_SAC_SMA_12_runs\',P_change_label{pc},'DP',T_change_label{tc},'DT\11ObsInflows\TRIAL_',num2str(trial_all(itrial)),'_SACSMA_Inflow_11stations.txt'],'w');
            fid = fopen(['C:\Users\Wi\Dropbox\CalLite\file transfer\11ObsInflows\TRIAL_',num2str(trial_all(itrial)),'_SACSMA_Inflow_11stations.txt'],'w');
            fprintf(fid,'%s\n','Year, Month, AMF, BLB, BND, FTO, MRC, SHA, SJF, SNS, TLG, TNL, YRS');
            fprintf(fid,'%4.0f, %2.0f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f\n',[datemat_mon mon_simflow_tot]');
            fclose(fid);
            
            waitbar(itrial/length(trial_all),h,[num2str(itrial),'/',num2str(length(trial_all))]);
%             toc;
        end
        
        close(h)
        
%        end
    end
    
end

toc;



%%

site = {'BearRiver','CacheCreek','CalaverasRiver','ChowchillaRiver','ConsumnesRiver','FresnoRiver','MokelumneRiver','PutahCreek','StonyCreek'};

area = [730.77176,2455.821437,927.537363,685.393248,1384.751552,717.246675,1408.7368,1477.609964,1926.664012];


% genyr = 87;
genyr = 55;

strdate = [1921 10 1];
% enddate = [2008 9 30];
enddate = [1976 9 30];
datemat = datevec(datenum(strdate):datenum(enddate));

datemat_mon = [];
for i=1927:1976%2008
    datemat_mon = [datemat_mon; repmat(i,12,1),[10,11,12,(1:9)]'];
end

% datemat2= datemat;

datemat(datemat(:,2) == 2 & datemat(:,3) == 29,:) = [];


for tc = 1%1:9
    
    for pc = 4%1:7
        
%        if ~(tc == 1 && pc == 4)
         
        mon_datemat = [];
        for i = strdate(1):enddate(1)
            mon_datemat = [mon_datemat;[repmat(i,12,1) (1:12)']];
        end
        
        mon_datemat_day = zeros(size(mon_datemat,1),1);
        for i = 1:size(mon_datemat,1)
            mon_datemat_day(i) = eomday(1921,mon_datemat(i,2));
        end
        mon_datemat_day = mon_datemat_day(10:end-3);
        
        h = waitbar(0,'Please Wait...');
        for itrial = 1:length(trial_all)
%             tic;
            mon_simflow_tot = [];
            for isite = 1:length(site)
                
                file = ls(['E:\CALIFORNIA\Kevin_dwr\Model Parameter File\calibration\9Unimpaired\mpiresult_SACSMA_DS_',site{isite},'*']);
                
                gridinfo = load(['E:\CALIFORNIA\Kevin_dwr\Model Parameter File\GridInfo\9Unimpaired\GridInfo_',site{isite},'.txt']);
                grid_lat = gridinfo(:,1);
                grid_lon = gridinfo(:,2);
                grid_area = gridinfo(:,3);
				grid_elev = gridinfo(:,4);
                grid_flen = gridinfo(:,5);
                
                
                numtotgrid = length(grid_lat);
                
                fid = fopen(['E:\CALIFORNIA\Kevin_dwr\Model Parameter File\calibration\9Unimpaired\',deblank(file)]);
                
                hru_par = [];
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
                
                
                FLOW = zeros(genyr*365,length(gridinfo));
                
                syr = strdate(1);
                eyr = enddate(1);
                
                
                for i = 1:length(gridinfo)
                    climate = load(['C:\Users\Wi\Dropbox\CalLite\file transfer\WGEN_SIMULATIONS_MAY2017\TRIAL_',num2str(trial_all(itrial)),'\TRIAL_',num2str(trial_all(itrial)),'_D.P.MEAN_1_D.T.MEAN_1\data_',num2str(grid_lat(i),'%2.4f'),'_',num2str(grid_lon(i),'%2.4f')]);
                    pr = climate(1:size(datemat,1),4)*P_change(pc);
                    tas = climate(1:size(datemat,1),5)+0.33+T_change(tc);
                    
                    selind = 0;
                    for j=1:numtotgrid
                        if strcmp(hru_par{j,1}{1},num2str(grid_lat(i),'%2.6f')) && strcmp(hru_par{j,2}{1},num2str(grid_lon(i),'%2.6f'))
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
                        pet(n:n+364) = pet_hamon([1921 10 1], [1922 9 30], tas(n:n+364), grid_lat(i), Coeff);
                        n = n +365;
                    end
					
					snow_outflow = nan(genyr*365,1);
					snow_melt = nan(genyr*365,1);
					snow_swe = nan(genyr*365,1);
                    snow_inistate = [0 0 0 0];
					n = 1;
                    for j = 1:genyr
                        [snow_outflow(n:n+364), snow_melt(n:n+364), snow_swe(n:n+364), snow_inistate_new] = snow_snow17([1921 10 1], [1922 9 30], pr(n:n+364), tas(n:n+364), grid_elev(i), SnowPar, snow_inistate);
                        snow_inistate = snow_inistate_new;
                        n = n +365;
                    end
                   
                    [inflow_direct, inflow_base] = sma_sacsma(pet,snow_outflow,SMA_Par, inistate);
                    [runoff, baseflow] = routing_lohmann(inflow_direct, inflow_base, grid_flen(i), route_par, isOutlet(i));
                    hru_flow = runoff;
                    FLOW(:,i) = hru_flow * grid_area(i)/sum(grid_area);
                    
                end
                
                mon_simflow = grpstats(sum(FLOW,2),{datemat(:,1),datemat(:,2)},'sum');
                
                mon_simflow = mon_simflow/1233.48184/1000*(area(isite)*1000000)/1000;
                
                mon_simflow_tot = [mon_simflow_tot, mon_simflow(61:end)];
                
            end
            
            fid = fopen(['C:\Users\Wi\Dropbox\CalLite\file transfer\9Unimpaired\TRIAL_',num2str(trial_all(itrial)),'_SACSMA_Inflow_9unimpairedstation.txt'],'w');
            fprintf(fid,'%s\n','Year, Month, BearRiver, CacheCreek, CalaverasRiver, ChowchillaRiver, CosumnesRiver, FresnoRiver, MokelumneRiver, PutahCreek, StonyCreek');
            fprintf(fid,'%4.0f, %2.0f, %14.6f, %14.6f, %14.6f, %14.6f, %14.6f, %14.6f, %14.6f, %14.6f, %14.6f\n',[datemat_mon mon_simflow_tot]');
            fclose(fid);
            
            waitbar(itrial/length(trial_all),h,[num2str(itrial),'/',num2str(length(trial_all))]);
%             toc;
        end
        close(h)
        
%        end
        
    end
end

toc;






