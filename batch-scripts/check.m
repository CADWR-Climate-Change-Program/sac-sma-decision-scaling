

strdate = [1949 1 1];
% enddate = [2008 9 30];
enddate = [2003 9 30];
datemat = datevec(datenum(strdate):datenum(enddate));

% datemat2= datemat;

% datemat(datemat(:,2) == 2 & datemat(:,3) == 29,:) = [];

inistate = [0 0 5 5 5 0];


[~,~,areainfo] = xlsread('E:\CALIFORNIA\SubBasin\CalLiteInput_CAsubbasin.xlsx','Area');



mon_datemat = [];
for i = strdate(1):enddate(1)
    mon_datemat = [mon_datemat;[repmat(i,12,1) (1:12)']];
end

mon_datemat_day = zeros(size(mon_datemat,1),1);
for i = 1:size(mon_datemat,1)
    mon_datemat_day(i) = eomday(1921,mon_datemat(i,2));
end
mon_datemat_day = mon_datemat_day(10:end-3);


isite = 1;



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


FLOW = zeros(size(datemat,1),length(gridinfo));

syr = strdate(1);
eyr = enddate(1);


for i = 1:length(gridinfo)
    %for i = 1:length(gridinfo)
    %                     climate = load(['C:\Users\Wi\Dropbox\CalLite\file transfer\WGEN_SIMULATIONS_MAY2017\TRIAL_',num2str(trial_all(itrial)),'\TRIAL_',num2str(trial_all(itrial)),'_D.P.MEAN_1_D.T.MEAN_1\data_',num2str(grid_lat(i),'%2.4f'),'_',num2str(grid_lon(i),'%2.4f')]);
    climate = load(['E:\CALIFORNIA\ClimateForcing\data_extension_prcp_tas\data_',num2str(grid_lat(i),'%2.4f'),'_',num2str(grid_lon(i),'%2.4f')]);
    %                     pr = climate(1:size(datemat,1),4)*P_change(pc);
    %                     tas = climate(1:size(datemat,1),5)+0.33+T_change(tc);
    sind = find(climate(:,1)==strdate(1)&climate(:,2)==strdate(2)&climate(:,3)==strdate(3));
    eind = find(climate(:,1)==enddate(1)&climate(:,2)==enddate(2)&climate(:,3)==enddate(3));
    pr = climate(sind:eind,4);
    tas = climate(sind:eind,5);
    
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
    
    
   
    pet = pet_hamon(strdate, enddate, tas, grid_lat(i), Coeff);
    
   
    
    [snow_outflow, snow_melt, snow_swe] = snow_snow17(strdate, enddate, pr, tas, grid_elev(i), SnowPar, [0 0 0 0]);
    
   
    [inflow_direct, inflow_base] = sma_sacsma(pet,snow_outflow,SMA_Par, inistate);
    
    [runoff, baseflow] = routing_lohmann(inflow_direct, inflow_base, grid_flen(i), route_par, isOutlet(i));
    
    hru_flow = runoff;
    
    
    FLOW(:,i) = hru_flow * grid_area(i)/sum(grid_area);
    
end

mon_simflow = grpstats(sum(FLOW,2),{datemat(:,1),datemat(:,2)},'sum');