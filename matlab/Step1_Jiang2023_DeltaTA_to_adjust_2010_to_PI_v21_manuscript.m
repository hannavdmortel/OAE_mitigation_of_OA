
% ----------
% Jiang2023.m 
%
% Maps for OAE mitigation of OA manuscript by van de Mortel et al. (2025)
%
% Greg Pelletier (gjpelletier@gmail.com) 
% 13-May-2025
% ----------

% OceanSODA-ETHZ downloaded from https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0220059

clear all
clc

% ------------


addpath(genpath([pwd,'/greg/']))  
line_colors;	% run the line_colors.m script in matlab/greg to define RGB values for matab line colors in plotting
addpath(genpath([pwd,'/m_map/']))  
addpath(genpath([pwd,'/tight_subplot/']))  

outdir = [pwd,'\png'];

fn_temp_hist = [pwd,'\nc\median\Temperature_median_historical.nc'];
fn_sal_hist = [pwd,'\nc\median\Salinity_median_historical.nc'];
fn_talk_hist = [pwd,'\nc\median\TA_median_historical.nc'];
fn_dic_hist = [pwd,'\nc\median\DIC_median_historical.nc'];
fn_phtot_hist = [pwd,'\nc\median\pHT_median_historical.nc'];
fn_omara_hist = [pwd,'\nc\median\Aragonite_median_historical.nc'];
fn_omcal_hist = [pwd,'\nc\median\Calcite_median_historical.nc'];
fn_co3_hist = [pwd,'\nc\median\CO3_median_historical.nc'];
fn_pco2_hist = [pwd,'\nc\median\pCO2_median_historical.nc'];

fn_temp_ssp585 = [pwd,'\nc\median\Temperature_median_ssp585.nc'];
fn_sal_ssp585 = [pwd,'\nc\median\Salinity_median_ssp585.nc'];
fn_talk_ssp585 = [pwd,'\nc\median\TA_median_ssp585.nc'];
fn_dic_ssp585 = [pwd,'\nc\median\DIC_median_ssp585.nc'];
fn_co3_ssp585 = [pwd,'\nc\median\CO3_median_ssp585.nc'];
fn_omara_ssp585 = [pwd,'\nc\median\Aragonite_median_ssp585.nc'];
fn_phtot_ssp585 = [pwd,'\nc\median\pHT_median_ssp585.nc'];
fn_pco2_ssp585 = [pwd,'\nc\median\pCO2_median_ssp585.nc'];


% - - -
% GLODAPv2.2023

load([pwd,'\GLODAPv2_2023_in_LME_v20250220.mat']);

% - - -
% interpolate dist2coast to GLODAP station locations
load([pwd,'\dist2coast_GLODAPv2_2023.mat']);
G2.dist2coast = Vq;
G2.in_dist100 = G2.dist2coast<=100;
G2.in_dist300 = G2.dist2coast<=300;

% ----------

% ----------
% ----------

% read data from Jiang et al 2023 nc files

% ----------
% ----------

% ----------

talk_hist = ncread(fn_talk_hist,'TA');
dic_hist = ncread(fn_dic_hist,'DIC');
temp_hist = ncread(fn_temp_hist,'temperature');
sal_hist = ncread(fn_sal_hist,'salinity');
phtot_hist = ncread(fn_phtot_hist,'pHT');
omara_hist = ncread(fn_omara_hist,'aragonite');
omcal_hist = ncread(fn_omcal_hist,'calcite');
co3_hist = ncread(fn_co3_hist,'carbonate');
pco2_hist = ncread(fn_pco2_hist,'pCO2');

talk_ssp585 = ncread(fn_talk_ssp585,'TA');
dic_ssp585 = ncread(fn_dic_ssp585,'DIC');
co3_ssp585 = ncread(fn_co3_ssp585,'carbonate');
omara_ssp585 = ncread(fn_omara_ssp585,'aragonite');
phtot_ssp585 = ncread(fn_phtot_ssp585,'pHT');
pco2_ssp585 = ncread(fn_pco2_ssp585,'pCO2');
temp_ssp585 = ncread(fn_temp_ssp585,'temperature');
sal_ssp585 = ncread(fn_sal_ssp585,'salinity');

time_hist = ncread(fn_omara_hist,'time');
time_ssp585 = ncread(fn_omara_ssp585,'time');

lon = ncread(fn_omara_hist,'longitude');	% 360x180 array of longitudes of the grid cells
lat = ncread(fn_omara_hist,'latitude');		% 360x180 array of latitudes of the grid cells

% - - -
% 1750
talk_1750 = squeeze(talk_hist(:,:,1));
dic_1750 = squeeze(dic_hist(:,:,1));
temp_1750 = squeeze(temp_hist(:,:,1));
sal_1750 = squeeze(sal_hist(:,:,1));

phtot_1750 = squeeze(phtot_hist(:,:,1));
pco2_1750 = squeeze(pco2_hist(:,:,1));
omara_1750 = squeeze(omara_hist(:,:,1));
co3_1750 = squeeze(co3_hist(:,:,1));

% - - -
% 2010
talk_2010 = squeeze(talk_hist(:,:,18));
dic_2010 = squeeze(dic_hist(:,:,18));
temp_2010 = squeeze(temp_hist(:,:,18));
sal_2010 = squeeze(sal_hist(:,:,18));

phtot_2010 = squeeze(phtot_hist(:,:,18));
pco2_2010 = squeeze(pco2_hist(:,:,18));
omara_2010 = squeeze(omara_hist(:,:,18));
co3_2010 = squeeze(co3_hist(:,:,18));

% - - -
% 2050
talk_2050 = squeeze(talk_ssp585(:,:,4));
dic_2050 = squeeze(dic_ssp585(:,:,4));
temp_2050 = squeeze(temp_ssp585(:,:,4));
sal_2050 = squeeze(sal_ssp585(:,:,4));

co3_2050 = squeeze(co3_ssp585(:,:,4));
omara_2050 = squeeze(omara_ssp585(:,:,4));
phtot_2050 = squeeze(phtot_ssp585(:,:,4));
pco2_2050 = squeeze(pco2_ssp585(:,:,4));

% - - -
% 2100
talk_2100 = squeeze(talk_ssp585(:,:,9));
dic_2100 = squeeze(dic_ssp585(:,:,9));
temp_2100 = squeeze(temp_ssp585(:,:,9));
sal_2100 = squeeze(sal_ssp585(:,:,9));

co3_2100 = squeeze(co3_ssp585(:,:,9));
omara_2100 = squeeze(omara_ssp585(:,:,9));
phtot_2100 = squeeze(phtot_ssp585(:,:,9));
pco2_2100 = squeeze(pco2_ssp585(:,:,9));

talk2dic_1750 = talk_1750 ./ dic_1750;
talk2dic_2010 = talk_2010 ./ dic_2010;
talk2dic_2050 = talk_2050 ./ dic_2050;
talk2dic_2100 = talk_2100 ./ dic_2100;


% get dimensions
[NX NY] = size(talk_1750);


% % - - -
% % read bathymetry data
% fn_etopo = 'C:\z\data\Bathymetry_Etopo5\Etopo5_1deg.nc';

% bath = ncread(fn_etopo,'bath');
% lon_bath = ncread(fn_etopo,'lon');
% lat_bath = ncread(fn_etopo,'lat');

% % convert to lon axis 20.5:379.5 (same format as Jiang et al arrays)
% bath2 = nan(numel(lon_bath),numel(lat_bath));
% bath2(1:340,:) = bath(21:360,:);
% bath2(341:360,:) = bath(1:20,:);

% - - -
% read distance to nearest coast
fn_dist2coast = [pwd,'\dist2coast_1deg.nc'];

dist2coast = ncread(fn_dist2coast,'dist2coast');
lon_dist2coast = ncread(fn_dist2coast,'lon');
lat_dist2coast = ncread(fn_dist2coast,'lat');

% convert to lon axis 20.5:379.5 (same format as Jiang et al arrays)
dist2 = nan(numel(lon_dist2coast),numel(lat_dist2coast));
dist2(1:340,:) = dist2coast(21:360,:);
dist2(341:360,:) = dist2coast(1:20,:);

dist.lon = lon;
dist.lat = lat;
dist.dist2coast = dist2;



% - - -
% - - -
% - - -

% regions: dist300, global, and polar

% - - -
% - - -
% - - -

% - - -
% dist300 - 300km distance from the coast for the coastal mask

% distance from coast map
% mapview(dist2)

lim_dist = 300;
mask_KT = squeeze(mean(talk_hist,3));
mask_KT(~isnan(mask_KT)) = 1;
mask_KT_dist300 = mask_KT;
mask_KT_dist300(dist2>lim_dist)=nan;
mapview(mask_KT_dist300)
in_dist300 = ~isnan(mask_KT_dist300);

% - - -
% polar - north of 60N or south of 60S

mask_KT = squeeze(mean(talk_hist,3));
mask_KT(~isnan(mask_KT)) = 1;
mask_KT_polar = mask_KT;
mask_KT_polar(lat>-60 & lat<60)=nan;
mapview(mask_KT_polar)
in_polar = ~isnan(mask_KT_polar);


% - - -
% global - all ocean surface grid cells

mask_KT = squeeze(mean(talk_hist,3));
mask_KT(~isnan(mask_KT)) = 1;
mask_KT_global = mask_KT;
mapview(mask_KT_global) 
in_global = ~isnan(mask_KT_global);





% - - -
% WOA2018 silicate and phosphate (used by Jiang et al)

load('woa2018_po4.mat');
load('woa2018_sio3.mat');


% - - -
% - - -
% - - -
% CO2SYS calcs
% - - -
% - - -
% - - -
disp('CO2SYS calcs of 1750, 2010, and 2100 variables')

% - - -
% 1750

% STEP 1: control conditions
% input vectors for co2sys function
vec_dic = reshape(dic_1750,[],1);
vec_talk = reshape(talk_1750,[],1);
vec_sal = reshape(sal_1750,[],1);
vec_temp = reshape(temp_1750,[],1);
vec_pres = zeros(size(vec_dic));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic,vec_talk,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% control conditions of driver variables are defined at this step
vec_dic_control = vec_dic;
vec_talk_control = vec_talk;
vec_pco2_control = vec_pco2;
% reshape back to NX x NY (360 x 180) grids
co2sys.phtot_1750 = reshape(vec_phtot,NX,[]);
co2sys.pco2_1750 = reshape(vec_pco2,NX,[]);
co2sys.fco2_1750 = reshape(vec_fco2,NX,[]);
co2sys.hco3_1750 = reshape(vec_hco3,NX,[]);
co2sys.co3_1750 = reshape(vec_co3,NX,[]);
co2sys.revelle_1750 = reshape(vec_revelle,NX,[]);
co2sys.omcal_1750 = reshape(vec_omcal,NX,[]);
co2sys.omara_1750 = reshape(vec_omara,NX,[]);

% ----------
% STEP 2a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 1; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_1750_dtalk_1 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_1750_dtalk_1 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_1750_dtalk_1 = reshape(vec_etamax,NX,[]);
co2sys.talk_1750_dtalk_1_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_1750_dtalk_1_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_1750_dtalk_1_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_1750_dtalk_1_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_1750_dtalk_1_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_1750_dtalk_1_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_1750_dtalk_1_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_1750_dtalk_1_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_1750_dtalk_1_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_1750_dtalk_1_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 2b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_1750_dtalk_1_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_1750_dtalk_1_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_1750_dtalk_1_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_1750_dtalk_1_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_1750_dtalk_1_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_1750_dtalk_1_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_1750_dtalk_1_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_1750_dtalk_1_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_1750_dtalk_1_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_1750_dtalk_1_cdreff_80 = reshape(vec_omara,NX,[]);

% ----------
% STEP 3a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 10; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_1750_dtalk_10 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_1750_dtalk_10 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_1750_dtalk_10 = reshape(vec_etamax,NX,[]);
co2sys.talk_1750_dtalk_10_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_1750_dtalk_10_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_1750_dtalk_10_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_1750_dtalk_10_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_1750_dtalk_10_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_1750_dtalk_10_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_1750_dtalk_10_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_1750_dtalk_10_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_1750_dtalk_10_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_1750_dtalk_10_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 3b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_1750_dtalk_10_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_1750_dtalk_10_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_1750_dtalk_10_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_1750_dtalk_10_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_1750_dtalk_10_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_1750_dtalk_10_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_1750_dtalk_10_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_1750_dtalk_10_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_1750_dtalk_10_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_1750_dtalk_10_cdreff_80 = reshape(vec_omara,NX,[]);

% ----------
% STEP 4a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 100; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_1750_dtalk_100 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_1750_dtalk_100 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_1750_dtalk_100 = reshape(vec_etamax,NX,[]);
co2sys.talk_1750_dtalk_100_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_1750_dtalk_100_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_1750_dtalk_100_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_1750_dtalk_100_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_1750_dtalk_100_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_1750_dtalk_100_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_1750_dtalk_100_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_1750_dtalk_100_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_1750_dtalk_100_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_1750_dtalk_100_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 4b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_1750_dtalk_100_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_1750_dtalk_100_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_1750_dtalk_100_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_1750_dtalk_100_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_1750_dtalk_100_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_1750_dtalk_100_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_1750_dtalk_100_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_1750_dtalk_100_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_1750_dtalk_100_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_1750_dtalk_100_cdreff_80 = reshape(vec_omara,NX,[]);


% - - -
% 2010

% STEP 1: control conditions
% input vectors for co2sys function
vec_dic = reshape(dic_2010,[],1);
vec_talk = reshape(talk_2010,[],1);
vec_sal = reshape(sal_2010,[],1);
vec_temp = reshape(temp_2010,[],1);
vec_pres = zeros(size(vec_dic));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic,vec_talk,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% control conditions of driver variables are defined at this step
vec_dic_control = vec_dic;
vec_talk_control = vec_talk;
vec_pco2_control = vec_pco2;
% reshape back to NX x NY (360 x 180) grids
co2sys.phtot_2010 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2010 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2010 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2010 = reshape(vec_hco3,NX,[]);
co2sys.co3_2010 = reshape(vec_co3,NX,[]);
co2sys.revelle_2010 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2010 = reshape(vec_omcal,NX,[]);
co2sys.omara_2010 = reshape(vec_omara,NX,[]);

% ----------
% STEP 2a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 1; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_2010_dtalk_1 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_2010_dtalk_1 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_2010_dtalk_1 = reshape(vec_etamax,NX,[]);
co2sys.talk_2010_dtalk_1_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2010_dtalk_1_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_2010_dtalk_1_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2010_dtalk_1_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_2010_dtalk_1_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2010_dtalk_1_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_2010_dtalk_1_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_2010_dtalk_1_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2010_dtalk_1_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_2010_dtalk_1_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 2b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_2010_dtalk_1_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2010_dtalk_1_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_2010_dtalk_1_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2010_dtalk_1_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2010_dtalk_1_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2010_dtalk_1_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_2010_dtalk_1_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_2010_dtalk_1_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2010_dtalk_1_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_2010_dtalk_1_cdreff_80 = reshape(vec_omara,NX,[]);

% ----------
% STEP 3a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 10; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_2010_dtalk_10 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_2010_dtalk_10 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_2010_dtalk_10 = reshape(vec_etamax,NX,[]);
co2sys.talk_2010_dtalk_10_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2010_dtalk_10_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_2010_dtalk_10_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2010_dtalk_10_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_2010_dtalk_10_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2010_dtalk_10_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_2010_dtalk_10_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_2010_dtalk_10_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2010_dtalk_10_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_2010_dtalk_10_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 3b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_2010_dtalk_10_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2010_dtalk_10_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_2010_dtalk_10_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2010_dtalk_10_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2010_dtalk_10_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2010_dtalk_10_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_2010_dtalk_10_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_2010_dtalk_10_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2010_dtalk_10_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_2010_dtalk_10_cdreff_80 = reshape(vec_omara,NX,[]);

% ----------
% STEP 4a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 100; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_2010_dtalk_100 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_2010_dtalk_100 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_2010_dtalk_100 = reshape(vec_etamax,NX,[]);
co2sys.talk_2010_dtalk_100_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2010_dtalk_100_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_2010_dtalk_100_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2010_dtalk_100_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_2010_dtalk_100_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2010_dtalk_100_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_2010_dtalk_100_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_2010_dtalk_100_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2010_dtalk_100_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_2010_dtalk_100_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 4b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_2010_dtalk_100_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2010_dtalk_100_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_2010_dtalk_100_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2010_dtalk_100_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2010_dtalk_100_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2010_dtalk_100_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_2010_dtalk_100_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_2010_dtalk_100_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2010_dtalk_100_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_2010_dtalk_100_cdreff_80 = reshape(vec_omara,NX,[]);

% - - -
% 2050

% STEP 1: control conditions
% input vectors for co2sys function
vec_dic = reshape(dic_2050,[],1);
vec_talk = reshape(talk_2050,[],1);
vec_sal = reshape(sal_2050,[],1);
vec_temp = reshape(temp_2050,[],1);
vec_pres = zeros(size(vec_dic));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic,vec_talk,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% control conditions of driver variables are defined at this step
vec_dic_control = vec_dic;
vec_talk_control = vec_talk;
vec_pco2_control = vec_pco2;
% reshape back to NX x NY (360 x 180) grids
co2sys.phtot_2050 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2050 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2050 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2050 = reshape(vec_hco3,NX,[]);
co2sys.co3_2050 = reshape(vec_co3,NX,[]);
co2sys.revelle_2050 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2050 = reshape(vec_omcal,NX,[]);
co2sys.omara_2050 = reshape(vec_omara,NX,[]);

% ----------
% STEP 2a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 1; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_2050_dtalk_1 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_2050_dtalk_1 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_2050_dtalk_1 = reshape(vec_etamax,NX,[]);
co2sys.talk_2050_dtalk_1_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2050_dtalk_1_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_2050_dtalk_1_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2050_dtalk_1_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_2050_dtalk_1_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2050_dtalk_1_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_2050_dtalk_1_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_2050_dtalk_1_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2050_dtalk_1_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_2050_dtalk_1_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 2b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_2050_dtalk_1_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2050_dtalk_1_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_2050_dtalk_1_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2050_dtalk_1_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2050_dtalk_1_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2050_dtalk_1_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_2050_dtalk_1_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_2050_dtalk_1_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2050_dtalk_1_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_2050_dtalk_1_cdreff_80 = reshape(vec_omara,NX,[]);

% ----------
% STEP 3a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 10; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_2050_dtalk_10 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_2050_dtalk_10 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_2050_dtalk_10 = reshape(vec_etamax,NX,[]);
co2sys.talk_2050_dtalk_10_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2050_dtalk_10_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_2050_dtalk_10_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2050_dtalk_10_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_2050_dtalk_10_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2050_dtalk_10_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_2050_dtalk_10_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_2050_dtalk_10_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2050_dtalk_10_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_2050_dtalk_10_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 3b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_2050_dtalk_10_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2050_dtalk_10_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_2050_dtalk_10_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2050_dtalk_10_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2050_dtalk_10_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2050_dtalk_10_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_2050_dtalk_10_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_2050_dtalk_10_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2050_dtalk_10_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_2050_dtalk_10_cdreff_80 = reshape(vec_omara,NX,[]);

% ----------
% STEP 4a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 100; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_2050_dtalk_100 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_2050_dtalk_100 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_2050_dtalk_100 = reshape(vec_etamax,NX,[]);
co2sys.talk_2050_dtalk_100_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2050_dtalk_100_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_2050_dtalk_100_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2050_dtalk_100_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_2050_dtalk_100_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2050_dtalk_100_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_2050_dtalk_100_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_2050_dtalk_100_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2050_dtalk_100_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_2050_dtalk_100_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 4b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_2050_dtalk_100_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2050_dtalk_100_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_2050_dtalk_100_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2050_dtalk_100_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2050_dtalk_100_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2050_dtalk_100_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_2050_dtalk_100_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_2050_dtalk_100_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2050_dtalk_100_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_2050_dtalk_100_cdreff_80 = reshape(vec_omara,NX,[]);

% - - -
% 2100

% STEP 1: control conditions
% input vectors for co2sys function
vec_dic = reshape(dic_2100,[],1);
vec_talk = reshape(talk_2100,[],1);
vec_sal = reshape(sal_2100,[],1);
vec_temp = reshape(temp_2100,[],1);
vec_pres = zeros(size(vec_dic));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic,vec_talk,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% control conditions of driver variables are defined at this step
vec_dic_control = vec_dic;
vec_talk_control = vec_talk;
vec_pco2_control = vec_pco2;
% reshape back to NX x NY (360 x 180) grids
co2sys.phtot_2100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2100 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2100 = reshape(vec_hco3,NX,[]);
co2sys.co3_2100 = reshape(vec_co3,NX,[]);
co2sys.revelle_2100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2100 = reshape(vec_omcal,NX,[]);
co2sys.omara_2100 = reshape(vec_omara,NX,[]);

% ----------
% STEP 2a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 1; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_2100_dtalk_1 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_2100_dtalk_1 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_2100_dtalk_1 = reshape(vec_etamax,NX,[]);
co2sys.talk_2100_dtalk_1_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2100_dtalk_1_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_2100_dtalk_1_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2100_dtalk_1_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_2100_dtalk_1_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2100_dtalk_1_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_2100_dtalk_1_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_2100_dtalk_1_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2100_dtalk_1_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_2100_dtalk_1_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 2b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_2100_dtalk_1_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2100_dtalk_1_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_2100_dtalk_1_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2100_dtalk_1_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2100_dtalk_1_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2100_dtalk_1_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_2100_dtalk_1_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_2100_dtalk_1_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2100_dtalk_1_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_2100_dtalk_1_cdreff_80 = reshape(vec_omara,NX,[]);

% ----------
% STEP 3a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 10; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_2100_dtalk_10 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_2100_dtalk_10 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_2100_dtalk_10 = reshape(vec_etamax,NX,[]);
co2sys.talk_2100_dtalk_10_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2100_dtalk_10_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_2100_dtalk_10_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2100_dtalk_10_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_2100_dtalk_10_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2100_dtalk_10_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_2100_dtalk_10_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_2100_dtalk_10_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2100_dtalk_10_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_2100_dtalk_10_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 3b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_2100_dtalk_10_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2100_dtalk_10_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_2100_dtalk_10_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2100_dtalk_10_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2100_dtalk_10_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2100_dtalk_10_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_2100_dtalk_10_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_2100_dtalk_10_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2100_dtalk_10_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_2100_dtalk_10_cdreff_80 = reshape(vec_omara,NX,[]);

% ----------
% STEP 4a: with delta_talk=1, use vec_pco2_control and talk_treatment to find cdrpot, etamax, and other variables with CDReff = 100% 
delta_talk = 100; 
vec_talk_treatment = vec_talk_control + delta_talk;		
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% pH total scale for CO2SYS "input" condition
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ delta_talk;
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.diceq_2100_dtalk_100 = reshape(vec_diceq,NX,[]);
co2sys.cdrpot_2100_dtalk_100 = reshape(vec_cdrpot,NX,[]);
co2sys.etamax_2100_dtalk_100 = reshape(vec_etamax,NX,[]);
co2sys.talk_2100_dtalk_100_cdreff_100 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2100_dtalk_100_cdreff_100 = reshape(vec_dic_control,NX,[]);	% assume treamtment with NaOH does not add DIC
co2sys.phtot_2100_dtalk_100_cdreff_100 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2100_dtalk_100_cdreff_100 = reshape(vec_pco2,NX,[]);	% pco2_treatment = pco2_control at CDReff=0
co2sys.fco2_2100_dtalk_100_cdreff_100 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2100_dtalk_100_cdreff_100 = reshape(vec_hco3,NX,[]);
co2sys.co3_2100_dtalk_100_cdreff_100 = reshape(vec_co3,NX,[]);
co2sys.revelle_2100_dtalk_100_cdreff_100 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2100_dtalk_100_cdreff_100 = reshape(vec_omcal,NX,[]);
co2sys.omara_2100_dtalk_100_cdreff_100 = reshape(vec_omara,NX,[]);

% STEP 4b: with delta_talk=1 and assuming CDReff=80%, 
cdreff = 0.8;
vec_dic_treatment = vec_dic_control + cdreff * vec_cdrpot;
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_dic_treatment,vec_talk_treatment,2,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
% carbonate system variables at treatment conditions with assumed cdreff
vec_phtot = RESULT(:,3); 	% pH total scale for CO2SYS "input" condition
vec_pco2 = RESULT(:,4); 	% pCO2 uatm for CO2SYS "input" condition
vec_fco2 = RESULT(:,5); 	% fCO2 uatm for CO2SYS "input" condition
vec_hco3 = RESULT(:,6); 	% HCO3^- umol/kg
vec_co3 = RESULT(:,7); 	% CO3^2- umol/kg
vec_revelle = RESULT(:,14); % Revelle factor for CO2SYS "input" condition
vec_omcal = RESULT(:,15); 	% Calcite saturation state Omega for CO2SYS "input" condition
vec_omara = RESULT(:,16); 	% Aragonite saturation state Omega for CO2SYS "input" condition
% reshape back to NX x NY (360 x 180) grids
co2sys.talk_2100_dtalk_100_cdreff_80 = reshape(vec_talk_treatment,NX,[]);
co2sys.dic_2100_dtalk_100_cdreff_80 = reshape(vec_dic_treatment,NX,[]);	% assume added DIC is due to realized CDR at assumed CDReff
co2sys.phtot_2100_dtalk_100_cdreff_80 = reshape(vec_phtot,NX,[]);
co2sys.pco2_2100_dtalk_100_cdreff_80 = reshape(vec_pco2,NX,[]);
co2sys.fco2_2100_dtalk_100_cdreff_80 = reshape(vec_fco2,NX,[]);
co2sys.hco3_2100_dtalk_100_cdreff_80 = reshape(vec_hco3,NX,[]);
co2sys.co3_2100_dtalk_100_cdreff_80 = reshape(vec_co3,NX,[]);
co2sys.revelle_2100_dtalk_100_cdreff_80 = reshape(vec_revelle,NX,[]);
co2sys.omcal_2100_dtalk_100_cdreff_80 = reshape(vec_omcal,NX,[]);
co2sys.omara_2100_dtalk_100_cdreff_80 = reshape(vec_omara,NX,[]);

% - - -
% dTAtrt = treatment delta_talk needed to change alkstar at time t to pre-industrial alkstar 
% - - -

% alkstar = TA-DIC
alkstar_1750 = talk_1750 - dic_1750; 
alkstar_2010 = talk_2010 - dic_2010; 
alkstar_2050 = talk_2050 - dic_2050; 
alkstar_2100 = talk_2100 - dic_2100; 

% - - -
% un-equilibrated delta_talk required for alkstar in 2010, 2050, 2100 to equal alkstar_1750

dTAtrt.unequil_2010_NaOH = alkstar_1750 - alkstar_2010;
dTAtrt.unequil_2050_NaOH = alkstar_1750 - alkstar_2050;
dTAtrt.unequil_2100_NaOH = alkstar_1750 - alkstar_2100;
dTAtrt.unequil_2010_Na2CO3 = 2 * (alkstar_1750 - alkstar_2010);
dTAtrt.unequil_2050_Na2CO3 = 2 * (alkstar_1750 - alkstar_2050);
dTAtrt.unequil_2100_Na2CO3 = 2 * (alkstar_1750 - alkstar_2100);

% mapview(dTAtrt.unequil_2010_NaOH)
% mapview(dTAtrt.unequil_2010_Na2CO3)

% - - -
% NaOH treatment with CDReff 80% equilibrated delta_talk required for alkstar in 2010, 2050, 2100 to equal alkstar_1750

% 2010
cdreff = 0.8;
% step 1
dTAtrt.equil_etamax_2010_NaOH_step1 = co2sys.etamax_2010_dtalk_1;
dTAtrt.equil_2010_NaOH_step1 = (alkstar_1750 - alkstar_2010) ./ (1 - cdreff * dTAtrt.equil_etamax_2010_NaOH_step1);
% mapview(dTAtrt.equil_etamax_2010_NaOH_step1)
% mapview(dTAtrt.equil_2010_NaOH_step1)
% step 2
vec_delta_talk = reshape(dTAtrt.equil_2010_NaOH_step1,[],1); 
vec_talk_control = reshape(talk_2010,[],1);
vec_talk_treatment = vec_talk_control + vec_delta_talk;		
% vec_pco2_control = reshape(pco2_2010,[],1);
vec_pco2_control = reshape(co2sys.pco2_2010,[],1);
vec_dic_control = reshape(dic_2010,[],1);
vec_sal = reshape(sal_2010,[],1);
vec_temp = reshape(temp_2010,[],1);
vec_pres = zeros(size(vec_sal));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% DIC at equilibrium for CO2SYS "input" condition umol/kg
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ vec_delta_talk;
% reshape back to NX x NY (360 x 180) grids
dTAtrt.cdreff80_talk_2010_NaOH_step2 = reshape(vec_talk_treatment,NX,[]);
% dTAtrt.cdreff80_dic_2010_NaOH_step2 = reshape(vec_diceq,NX,[]);
dTAtrt.cdreff80_dic_2010_NaOH_step2 = reshape(vec_dic_control + cdreff*vec_cdrpot,NX,[]);
dTAtrt.cdreff80_alkstar_2010_NaOH_step2 = dTAtrt.cdreff80_talk_2010_NaOH_step2 - dTAtrt.cdreff80_dic_2010_NaOH_step2; 
dTAtrt.cdreff80_etamax_2010_NaOH_step2 = reshape(vec_etamax,NX,[]);
dTAtrt.cdreff80_2010_NaOH_step2 = (alkstar_1750 - alkstar_2010) ./ (1 - cdreff * dTAtrt.cdreff80_etamax_2010_NaOH_step2);

% step 3
vec_delta_talk = reshape(dTAtrt.cdreff80_2010_NaOH_step2,[],1); 
vec_talk_control = reshape(talk_2010,[],1);
vec_talk_treatment = vec_talk_control + vec_delta_talk;		
% vec_pco2_control = reshape(pco2_2010,[],1);
vec_pco2_control = reshape(co2sys.pco2_2010,[],1);
vec_dic_control = reshape(dic_2010,[],1);
vec_sal = reshape(sal_2010,[],1);
vec_temp = reshape(temp_2010,[],1);
vec_pres = zeros(size(vec_sal));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% DIC at equilibrium for CO2SYS "input" condition umol/kg
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ vec_delta_talk;
% reshape back to NX x NY (360 x 180) grids
dTAtrt.cdreff80_talk_2010_NaOH_step3 = reshape(vec_talk_treatment,NX,[]);
% dTAtrt.cdreff80_dic_2010_NaOH_step3 = reshape(vec_diceq,NX,[]);
dTAtrt.cdreff80_dic_2010_NaOH_step3 = reshape(vec_dic_control + cdreff*vec_cdrpot,NX,[]);
dTAtrt.cdreff80_alkstar_2010_NaOH_step3 = dTAtrt.cdreff80_talk_2010_NaOH_step3 - dTAtrt.cdreff80_dic_2010_NaOH_step3; 
dTAtrt.cdreff80_etamax_2010_NaOH_step3 = reshape(vec_etamax,NX,[]);
dTAtrt.cdreff80_2010_NaOH_step3 = (alkstar_1750 - alkstar_2010) ./ (1 - cdreff * dTAtrt.cdreff80_etamax_2010_NaOH_step3);

% step 4
vec_delta_talk = reshape(dTAtrt.cdreff80_2010_NaOH_step3,[],1); 
vec_talk_control = reshape(talk_2010,[],1);
vec_talk_treatment = vec_talk_control + vec_delta_talk;		
% vec_pco2_control = reshape(pco2_2010,[],1);
vec_pco2_control = reshape(co2sys.pco2_2010,[],1);
vec_dic_control = reshape(dic_2010,[],1);
vec_sal = reshape(sal_2010,[],1);
vec_temp = reshape(temp_2010,[],1);
vec_pres = zeros(size(vec_sal));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% DIC at equilibrium for CO2SYS "input" condition umol/kg
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ vec_delta_talk;
% reshape back to NX x NY (360 x 180) grids
dTAtrt.cdreff80_talk_2010_NaOH_step4 = reshape(vec_talk_treatment,NX,[]);
% dTAtrt.cdreff80_dic_2010_NaOH_step4 = reshape(vec_diceq,NX,[]);
dTAtrt.cdreff80_dic_2010_NaOH_step4 = reshape(vec_dic_control + cdreff*vec_cdrpot,NX,[]);
dTAtrt.cdreff80_alkstar_2010_NaOH_step4 = dTAtrt.cdreff80_talk_2010_NaOH_step4 - dTAtrt.cdreff80_dic_2010_NaOH_step4; 
dTAtrt.cdreff80_etamax_2010_NaOH_step4 = reshape(vec_etamax,NX,[]);
dTAtrt.cdreff80_2010_NaOH_step4 = (alkstar_1750 - alkstar_2010) ./ (1 - cdreff * dTAtrt.cdreff80_etamax_2010_NaOH_step4);

% step 5
vec_delta_talk = reshape(dTAtrt.cdreff80_2010_NaOH_step4,[],1); 
vec_talk_control = reshape(talk_2010,[],1);
vec_talk_treatment = vec_talk_control + vec_delta_talk;		
% vec_pco2_control = reshape(pco2_2010,[],1);
vec_pco2_control = reshape(co2sys.pco2_2010,[],1);
vec_dic_control = reshape(dic_2010,[],1);
vec_sal = reshape(sal_2010,[],1);
vec_temp = reshape(temp_2010,[],1);
vec_pres = zeros(size(vec_sal));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% DIC at equilibrium for CO2SYS "input" condition umol/kg
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ vec_delta_talk;
% reshape back to NX x NY (360 x 180) grids
dTAtrt.cdreff80_talk_2010_NaOH_step5 = reshape(vec_talk_treatment,NX,[]);
% dTAtrt.cdreff80_dic_2010_NaOH_step5 = reshape(vec_diceq,NX,[]);
dTAtrt.cdreff80_dic_2010_NaOH_step5 = reshape(vec_dic_control + cdreff*vec_cdrpot,NX,[]);
dTAtrt.cdreff80_alkstar_2010_NaOH_step5 = dTAtrt.cdreff80_talk_2010_NaOH_step5 - dTAtrt.cdreff80_dic_2010_NaOH_step5; 
dTAtrt.cdreff80_etamax_2010_NaOH_step5 = reshape(vec_etamax,NX,[]);
dTAtrt.cdreff80_2010_NaOH_step5 = (alkstar_1750 - alkstar_2010) ./ (1 - cdreff * dTAtrt.cdreff80_etamax_2010_NaOH_step5);

% step2 check (should be 0)
alkstar_test = dTAtrt.cdreff80_alkstar_2010_NaOH_step2;
alkstar_dif_test = (alkstar_1750 - alkstar_test);
mapview(alkstar_dif_test)

% step3 check (should be 0)
alkstar_test = dTAtrt.cdreff80_alkstar_2010_NaOH_step3;
alkstar_dif_test = (alkstar_1750 - alkstar_test);
mapview(alkstar_dif_test)

% step4 check (should be 0)
alkstar_test = dTAtrt.cdreff80_alkstar_2010_NaOH_step4;
alkstar_dif_test = (alkstar_1750 - alkstar_test);
mapview(alkstar_dif_test)

% step5 check (should be 0)
alkstar_test = dTAtrt.cdreff80_alkstar_2010_NaOH_step5;
alkstar_dif_test = (alkstar_1750 - alkstar_test);
mapview(alkstar_dif_test)
difsq = alkstar_dif_test .^ 2;
rmse = sqrt(nanmean2(difsq(:)))
rpd = 100 * rmse / nanmean2(alkstar_1750(:));

% - - -
% Na2CO3 treatment with CDReff 80% equilibrated delta_talk required for alkstar in 2010, 2050, 2100 to equal alkstar_1750

% 2010
cdreff = 0.8;
% step 1
dTAtrt.equil_etamax_2010_Na2CO3_step1 = co2sys.etamax_2010_dtalk_1;
dTAtrt.equil_2010_Na2CO3_step1 = (alkstar_1750 - alkstar_2010) ./ (0.5 - cdreff * (dTAtrt.equil_etamax_2010_Na2CO3_step1-0.5));
% mapview(dTAtrt.equil_etamax_2010_Na2CO3_step1)
% mapview(dTAtrt.equil_2010_Na2CO3_step1)
% step 2
vec_delta_talk = reshape(dTAtrt.equil_2010_Na2CO3_step1,[],1); 
vec_talk_control = reshape(talk_2010,[],1);
vec_talk_treatment = vec_talk_control + vec_delta_talk;		
% vec_pco2_control = reshape(pco2_2010,[],1);
vec_pco2_control = reshape(co2sys.pco2_2010,[],1);
vec_dic_control = reshape(dic_2010,[],1);
vec_sal = reshape(sal_2010,[],1);
vec_temp = reshape(temp_2010,[],1);
vec_pres = zeros(size(vec_sal));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% DIC at equilibrium for CO2SYS "input" condition umol/kg
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ vec_delta_talk;
vec_dic = vec_dic_control + 0.5 * vec_delta_talk + cdreff .* (vec_etamax - 0.5) .* vec_delta_talk;
% reshape back to NX x NY (360 x 180) grids
dTAtrt.cdreff80_talk_2010_Na2CO3_step2 = reshape(vec_talk_treatment,NX,[]);
% dTAtrt.cdreff80_dic_2010_Na2CO3_step2 = reshape(vec_diceq,NX,[]);
% dTAtrt.cdreff80_dic_2010_Na2CO3_step2 = reshape(vec_dic_control + cdreff*vec_cdrpot,NX,[]);
dTAtrt.cdreff80_dic_2010_Na2CO3_step2 = reshape(vec_dic,NX,[]);
dTAtrt.cdreff80_alkstar_2010_Na2CO3_step2 = dTAtrt.cdreff80_talk_2010_Na2CO3_step2 - dTAtrt.cdreff80_dic_2010_Na2CO3_step2; 
dTAtrt.cdreff80_etamax_2010_Na2CO3_step2 = reshape(vec_etamax,NX,[]);
dTAtrt.cdreff80_2010_Na2CO3_step2 = (alkstar_1750 - alkstar_2010) ./ (0.5 - cdreff * (dTAtrt.cdreff80_etamax_2010_Na2CO3_step2-0.5));

% step 3
vec_delta_talk = reshape(dTAtrt.cdreff80_2010_Na2CO3_step2,[],1); 
vec_talk_control = reshape(talk_2010,[],1);
vec_talk_treatment = vec_talk_control + vec_delta_talk;		
% vec_pco2_control = reshape(pco2_2010,[],1);
vec_pco2_control = reshape(co2sys.pco2_2010,[],1);
vec_dic_control = reshape(dic_2010,[],1);
vec_sal = reshape(sal_2010,[],1);
vec_temp = reshape(temp_2010,[],1);
vec_pres = zeros(size(vec_sal));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% DIC at equilibrium for CO2SYS "input" condition umol/kg
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ vec_delta_talk;
vec_dic = vec_dic_control + 0.5 * vec_delta_talk + cdreff .* (vec_etamax - 0.5) .* vec_delta_talk;
% reshape back to NX x NY (360 x 180) grids
dTAtrt.cdreff80_talk_2010_Na2CO3_step3 = reshape(vec_talk_treatment,NX,[]);
% dTAtrt.cdreff80_dic_2010_Na2CO3_step3 = reshape(vec_diceq,NX,[]);
% dTAtrt.cdreff80_dic_2010_Na2CO3_step3 = reshape(vec_dic_control + cdreff*vec_cdrpot,NX,[]);
dTAtrt.cdreff80_dic_2010_Na2CO3_step3 = reshape(vec_dic,NX,[]);
dTAtrt.cdreff80_alkstar_2010_Na2CO3_step3 = dTAtrt.cdreff80_talk_2010_Na2CO3_step3 - dTAtrt.cdreff80_dic_2010_Na2CO3_step3; 
dTAtrt.cdreff80_etamax_2010_Na2CO3_step3 = reshape(vec_etamax,NX,[]);
dTAtrt.cdreff80_2010_Na2CO3_step3 = (alkstar_1750 - alkstar_2010) ./ (0.5 - cdreff * (dTAtrt.cdreff80_etamax_2010_Na2CO3_step3-0.5));

% step 4
vec_delta_talk = reshape(dTAtrt.cdreff80_2010_Na2CO3_step3,[],1); 
vec_talk_control = reshape(talk_2010,[],1);
vec_talk_treatment = vec_talk_control + vec_delta_talk;		
% vec_pco2_control = reshape(pco2_2010,[],1);
vec_pco2_control = reshape(co2sys.pco2_2010,[],1);
vec_dic_control = reshape(dic_2010,[],1);
vec_sal = reshape(sal_2010,[],1);
vec_temp = reshape(temp_2010,[],1);
vec_pres = zeros(size(vec_sal));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% DIC at equilibrium for CO2SYS "input" condition umol/kg
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ vec_delta_talk;
vec_dic = vec_dic_control + 0.5 * vec_delta_talk + cdreff .* (vec_etamax - 0.5) .* vec_delta_talk;
% reshape back to NX x NY (360 x 180) grids
dTAtrt.cdreff80_talk_2010_Na2CO3_step4 = reshape(vec_talk_treatment,NX,[]);
% dTAtrt.cdreff80_dic_2010_Na2CO3_step4 = reshape(vec_diceq,NX,[]);
% dTAtrt.cdreff80_dic_2010_Na2CO3_step4 = reshape(vec_dic_control + cdreff*vec_cdrpot,NX,[]);
dTAtrt.cdreff80_dic_2010_Na2CO3_step4 = reshape(vec_dic,NX,[]);
dTAtrt.cdreff80_alkstar_2010_Na2CO3_step4 = dTAtrt.cdreff80_talk_2010_Na2CO3_step4 - dTAtrt.cdreff80_dic_2010_Na2CO3_step4; 
dTAtrt.cdreff80_etamax_2010_Na2CO3_step4 = reshape(vec_etamax,NX,[]);
dTAtrt.cdreff80_2010_Na2CO3_step4 = (alkstar_1750 - alkstar_2010) ./ (0.5 - cdreff * (dTAtrt.cdreff80_etamax_2010_Na2CO3_step4-0.5));

% step 5
vec_delta_talk = reshape(dTAtrt.cdreff80_2010_Na2CO3_step4,[],1); 
vec_talk_control = reshape(talk_2010,[],1);
vec_talk_treatment = vec_talk_control + vec_delta_talk;		
% vec_pco2_control = reshape(pco2_2010,[],1);
vec_pco2_control = reshape(co2sys.pco2_2010,[],1);
vec_dic_control = reshape(dic_2010,[],1);
vec_sal = reshape(sal_2010,[],1);
vec_temp = reshape(temp_2010,[],1);
vec_pres = zeros(size(vec_sal));
vec_si = reshape(sio3,[],1);
vec_po4 = reshape(po4,[],1);
[RESULT,HEADERS,NICEHEADERS]=CO2SYS_v21_PF87( ...
	vec_pco2_control,vec_talk_treatment,4,1, ...
	vec_sal, ...
	vec_temp, ones(size(vec_temp)) .* 20, ...
	vec_pres, ones(size(vec_pres)) .* 0, ...
	vec_si,vec_po4, ...
	1,10,3,2);  % 1=pH total scale, 10=K1K2 of Lueker et al 2000, 3=KSO4 of Dickson 1990a & TB of Lee 2010, 2= KF of Perez and Fraga 1987  
vec_diceq = RESULT(:,2); 	% DIC at equilibrium for CO2SYS "input" condition umol/kg
vec_cdrpot = vec_diceq - vec_dic_control;
vec_etamax = vec_cdrpot ./ vec_delta_talk;
vec_dic = vec_dic_control + 0.5 * vec_delta_talk + cdreff .* (vec_etamax - 0.5) .* vec_delta_talk;
% reshape back to NX x NY (360 x 180) grids
dTAtrt.cdreff80_talk_2010_Na2CO3_step5 = reshape(vec_talk_treatment,NX,[]);
% dTAtrt.cdreff80_dic_2010_Na2CO3_step5 = reshape(vec_diceq,NX,[]);
% dTAtrt.cdreff80_dic_2010_Na2CO3_step5 = reshape(vec_dic_control + cdreff*vec_cdrpot,NX,[]);
dTAtrt.cdreff80_dic_2010_Na2CO3_step5 = reshape(vec_dic,NX,[]);
dTAtrt.cdreff80_alkstar_2010_Na2CO3_step5 = dTAtrt.cdreff80_talk_2010_Na2CO3_step5 - dTAtrt.cdreff80_dic_2010_Na2CO3_step5; 
dTAtrt.cdreff80_etamax_2010_Na2CO3_step5 = reshape(vec_etamax,NX,[]);
dTAtrt.cdreff80_2010_Na2CO3_step5 = (alkstar_1750 - alkstar_2010) ./ (0.5 - cdreff * (dTAtrt.cdreff80_etamax_2010_Na2CO3_step5-0.5));

% step2 check (should be 0)
alkstar_test = dTAtrt.cdreff80_alkstar_2010_Na2CO3_step2;
alkstar_dif_test = (alkstar_1750 - alkstar_test);
mapview(alkstar_dif_test)

% step3 check (should be 0)
alkstar_test = dTAtrt.cdreff80_alkstar_2010_Na2CO3_step3;
alkstar_dif_test = (alkstar_1750 - alkstar_test);
mapview(alkstar_dif_test)

% step4 check (should be 0)
alkstar_test = dTAtrt.cdreff80_alkstar_2010_Na2CO3_step4;
alkstar_dif_test = (alkstar_1750 - alkstar_test);
mapview(alkstar_dif_test)

% step5 check (should be 0)
alkstar_test = dTAtrt.cdreff80_alkstar_2010_Na2CO3_step5;
alkstar_dif_test = (alkstar_1750 - alkstar_test);
mapview(alkstar_dif_test)




% - - -
% - - -
% - - -

% Make global maps of OAE dTA umol/kg required to restore TA-DIC in 2010 to 1750 conditions

Step2_Jiang2023_DeltaTA_maps_v21_manuscript

% - - -
% - - -
% - - -








% - - -
% - - -
% - - -
% Scatter plots of 1750-2010 co3 vs talk2dic
% - - -
% - - -
% - - -










% - - -
% - - -
% - - -
% - - -
% - - -
% - - -
 
% global - Regression of GLODAP observation data vs TA-DIC

% - - -
% - - -
% - - -
% - - -
% - - -
% - - -




% - - -
% - - -
% - - -

% global - Regression of CO3^2- vs TA-DIC

% - - -
% - - -
% - - -

deplim = 10;
markersize = 5;
markersize2 = 5;
fontsize_title = 16;
fontsize_lab = 12;
fontsize_ax = 12;
fontsize_eqn = 10;
fontsize_stats = 10;

figure(1)
hFig=gcf;
clf(hFig);

tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

hax1 = nexttile;

% Jiang data
Y1 = co3_1750;
Y1(~in_global) = nan;
Y2 = co3_2010;
Y2(~in_global) = nan;
Y3 = co3_2100;
Y3(~in_global) = nan;

X1 = talk_1750 -  dic_1750;
X1(~in_global) = nan;
X2 = talk_2010 - dic_2010;
X2(~in_global) = nan;
X3 = talk_2100 - dic_2100;
X3(~in_global) = nan;

X1 = reshape(X1,[],1);
X2 = reshape(X2,[],1);
X3 = reshape(X3,[],1);
Y1 = reshape(Y1,[],1);
Y2 = reshape(Y2,[],1);
Y3 = reshape(Y3,[],1);

% GLODAP obs
idx = find(G2.in_global & G2.depth<deplim);
Y4 = G2.co2sys_co3(idx);
X4 = G2.talk(idx) - G2.dic(idx);

	% regression with 95%CI
	x = X4;
	y = Y4;
	% eliminate nans from x and y
	idx = find(~isnan(x) & ~isnan(y));
	x = x(idx);
	y = y(idx);
	% create a function for the nonlinear regression model equation
	f = @(param,xval) param(1) + param(2) .*  xval + param(3) .* xval .^2 + param(4) .* xval .^ 3;
	% specify initial parameter values to use for the nonlinear regression optimization
	p_init = [nanmean2(y(:)) 0 0 0];
	% use nlinfit to find the optimum parameter values (popt) and covariance matrix (pcov)
	[popt,~,~,pcov] = nlinfit(x,y,f,p_init);					  
	% evenly spaced new observations that extend from below min(x) to above max(x)
	nlinspace = 100;
	x_new = linspace(nanmin2(x(:)), nanmax2(x(:)), nlinspace);
	% probability level for the prediction limits (e.g. alpha=0.05 for 95% prediction limits and 95% confidence limits)
	alpha=0.05;
	% use the delta_method.m function to calculate confidence intervals and prediction intervals
	d = delta_method(pcov,popt,x_new,f,x,y,alpha);

hold on

	patch([d.x_new' fliplr(d.x_new')],[d.upr_pred' fliplr(d.lwr_pred')],[.8 .8 .8],'edgecolor','none','facealpha',0.5)
	plot(x,y,'k.','markersize',5)
	plot(d.x_new,d.y_new,'k-','linewidth',1)
	b1 = num2str(popt(1),'%.2s');
	b2 = num2str(popt(2),'%.2s');
	b3 = num2str(popt(3),'%.2s');
	b4 = num2str(popt(4),'%.2s');
	eqn = ['y = ',b1,' + ',b2,' x + ',b3,' x^2 + ',b4,' x^3'];
	pstr = num2str(d.pvalue);
	rsq = num2str(d.rsquared,"%.3f");
	ser = num2str(d.syx,"%.1f");
	nobs = num2str(d.nobs);
	text(25,370,'GLODAP regression statistics','fontsize',fontsize_eqn);    
	text(25,350,eqn,'fontsize',fontsize_eqn);    
	text(25,325,['r^2=',rsq,', p=',pstr,', rmse=',ser,', n=',nobs],'fontsize',fontsize_stats);     

% plot(X1,Y1,'.','color',cmap_lines(1,:),'markersize',markersize)
% plot(X2,Y2,'o','color',cmap_lines(2,:),'markersize',markersize/3)
% plot(X3,Y3,'.','color',cmap_lines(5,:),'markersize',.75*markersize)
% plot(X4,Y4,'k.','markersize',markersize2)
ylabel(['CO_3^{2-} \mumol kg^{-1}'])
xlabel(['TA-DIC \mumol kg^{-1}'])
xlim([0 500])
ylim([0 400])
% legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','Jiang 1750','Jiang 2010','Jiang 2100 (SSP585)','location','southeast')
legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','location','southeast')
title(['CO_3^{2-} vs TA-DIC'])
hold off

% - - -
k=.8;
set(gcf, 'PaperPosition', [0 0 k*6 k*6])   
print(gcf, [pwd '/png/regression_global_co3_vs_alkstar_0-10m_GLODAP_v21.png'], '-dpng', '-r600' );  



% - - -
% - - -
% - - -

% global - Regression of omara vs TA-DIC

% - - -
% - - -
% - - -

deplim = 10;
markersize = 5;
markersize2 = 5;
fontsize_title = 16;
fontsize_lab = 12;
fontsize_ax = 12;
fontsize_eqn = 10;
fontsize_stats = 10;

figure(1)
hFig=gcf;
clf(hFig);

tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

hax1 = nexttile;

% Jiang data
Y1 = omara_1750;
Y1(~in_global) = nan;
Y2 = omara_2010;
Y2(~in_global) = nan;
Y3 = omara_2100;
Y3(~in_global) = nan;

X1 = talk_1750 -  dic_1750;
X1(~in_global) = nan;
X2 = talk_2010 - dic_2010;
X2(~in_global) = nan;
X3 = talk_2100 - dic_2100;
X3(~in_global) = nan;

X1 = reshape(X1,[],1);
X2 = reshape(X2,[],1);
X3 = reshape(X3,[],1);
Y1 = reshape(Y1,[],1);
Y2 = reshape(Y2,[],1);
Y3 = reshape(Y3,[],1);

% GLODAP obs
idx = find(G2.in_global & G2.depth<deplim);
Y4 = G2.co2sys_omara(idx);
X4 = G2.talk(idx) - G2.dic(idx);

	% regression with 95%CI
	x = X4;
	y = Y4;
	% eliminate nans from x and y
	idx = find(~isnan(x) & ~isnan(y));
	x = x(idx);
	y = y(idx);
	% create a function for the nonlinear regression model equation
	f = @(param,xval) param(1) + param(2) .*  xval + param(3) .* xval .^2 + param(4) .* xval .^ 3;
	% specify initial parameter values to use for the nonlinear regression optimization
	p_init = [nanmean2(y(:)) 0 0 0];
	% use nlinfit to find the optimum parameter values (popt) and covariance matrix (pcov)
	[popt,~,~,pcov] = nlinfit(x,y,f,p_init);					  
	% evenly spaced new observations that extend from below min(x) to above max(x)
	nlinspace = 100;
	x_new = linspace(nanmin2(x(:)), nanmax2(x(:)), nlinspace);
	% probability level for the prediction limits (e.g. alpha=0.05 for 95% prediction limits and 95% confidence limits)
	alpha=0.05;
	% use the delta_method.m function to calculate confidence intervals and prediction intervals
	d = delta_method(pcov,popt,x_new,f,x,y,alpha);

hold on

	patch([d.x_new' fliplr(d.x_new')],[d.upr_pred' fliplr(d.lwr_pred')],[.8 .8 .8],'edgecolor','none','facealpha',0.5)
	plot(x,y,'k.','markersize',5)
	plot(d.x_new,d.y_new,'k-','linewidth',1)
	b1 = num2str(popt(1),'%.2s');
	b2 = num2str(popt(2),'%.2s');
	b3 = num2str(popt(3),'%.2s');
	b4 = num2str(popt(4),'%.2s');
	eqn = ['y = ',b1,' + ',b2,' x + ',b3,' x^2 + ',b4,' x^3'];
	pstr = num2str(d.pvalue);
	rsq = num2str(d.rsquared,"%.3f");
	ser = num2str(d.syx,"%.2f");
	nobs = num2str(d.nobs);
	text(25,6,'GLODAP regression statistics','fontsize',fontsize_eqn);    
	text(25,5.6,eqn,'fontsize',fontsize_eqn);    
	text(25,5.2,['r^2=',rsq,', p=',pstr,', rmse=',ser,', n=',nobs],'fontsize',fontsize_stats);     

% plot(X1,Y1,'.','color',cmap_lines(1,:),'markersize',markersize)
% plot(X2,Y2,'o','color',cmap_lines(2,:),'markersize',markersize/3)
% plot(X3,Y3,'.','color',cmap_lines(5,:),'markersize',.75*markersize)
% plot(X4,Y4,'k.','markersize',markersize2)
ylabel(['\Omega_{ara}'])
xlabel(['TA-DIC \mumol kg^{-1}'])
xlim([0 500])
ylim([0 7])
% legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','Jiang 1750','Jiang 2010','Jiang 2100 (SSP585)','location','southeast')
legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','location','southeast')
title(['\Omega_{ara} vs TA-DIC'])
hold off

% - - -
k=.8;
set(gcf, 'PaperPosition', [0 0 k*6 k*6])   
print(gcf, [pwd '/png/regression_global_omara_vs_alkstar_0-10m_GLODAP_v21.png'], '-dpng', '-r600' );  



% - - -
% - - -
% - - -

% global - Regression of phtot vs TA-DIC

% - - -
% - - -
% - - -

deplim = 10;
markersize = 5;
markersize2 = 5;
fontsize_title = 16;
fontsize_lab = 12;
fontsize_ax = 12;
fontsize_eqn = 10;
fontsize_stats = 10;

figure(1)
hFig=gcf;
clf(hFig);

tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

hax1 = nexttile;

% Jiang data
Y1 = phtot_1750;
Y1(~in_global) = nan;
Y2 = phtot_2010;
Y2(~in_global) = nan;
Y3 = phtot_2100;
Y3(~in_global) = nan;

X1 = talk_1750 -  dic_1750;
X1(~in_global) = nan;
X2 = talk_2010 - dic_2010;
X2(~in_global) = nan;
X3 = talk_2100 - dic_2100;
X3(~in_global) = nan;

X1 = reshape(X1,[],1);
X2 = reshape(X2,[],1);
X3 = reshape(X3,[],1);
Y1 = reshape(Y1,[],1);
Y2 = reshape(Y2,[],1);
Y3 = reshape(Y3,[],1);

% GLODAP obs
idx = find(G2.in_global & G2.depth<deplim);
Y4 = G2.co2sys_phtot(idx);
X4 = G2.talk(idx) - G2.dic(idx);

	% regression with 95%CI
	x = X4;
	y = Y4;
	% eliminate nans from x and y
	idx = find(~isnan(x) & ~isnan(y));
	x = x(idx);
	y = y(idx);
	% create a function for the nonlinear regression model equation
	f = @(param,xval) param(1) + param(2) .*  xval + param(3) .* xval .^2 + param(4) .* xval .^ 3;
	% specify initial parameter values to use for the nonlinear regression optimization
	p_init = [nanmean2(y(:)) 0 0 0];
	% use nlinfit to find the optimum parameter values (popt) and covariance matrix (pcov)
	[popt,~,~,pcov] = nlinfit(x,y,f,p_init);					  
	% evenly spaced new observations that extend from below min(x) to above max(x)
	nlinspace = 100;
	x_new = linspace(nanmin2(x(:)), nanmax2(x(:)), nlinspace);
	% probability level for the prediction limits (e.g. alpha=0.05 for 95% prediction limits and 95% confidence limits)
	alpha=0.05;
	% use the delta_method.m function to calculate confidence intervals and prediction intervals
	d = delta_method(pcov,popt,x_new,f,x,y,alpha);

hold on

	patch([d.x_new' fliplr(d.x_new')],[d.upr_pred' fliplr(d.lwr_pred')],[.8 .8 .8],'edgecolor','none','facealpha',0.5)
	plot(x,y,'k.','markersize',5)
	plot(d.x_new,d.y_new,'k-','linewidth',1)
	b1 = num2str(popt(1),'%.2s');
	b2 = num2str(popt(2),'%.2s');
	b3 = num2str(popt(3),'%.2s');
	b4 = num2str(popt(4),'%.2s');
	eqn = ['y = ',b1,' + ',b2,' x + ',b3,' x^2 + ',b4,' x^3'];
	pstr = num2str(d.pvalue);
	rsq = num2str(d.rsquared,"%.3f");
	ser = num2str(d.syx,"%.2f");
	nobs = num2str(d.nobs);
	text(25,8.85,'GLODAP regression statistics','fontsize',fontsize_eqn);    
	text(25,8.75,eqn,'fontsize',fontsize_eqn);    
	text(25,8.65,['r^2=',rsq,', p=',pstr,', rmse=',ser,', n=',nobs],'fontsize',fontsize_stats);     

% plot(X1,Y1,'.','color',cmap_lines(1,:),'markersize',markersize)
% plot(X2,Y2,'o','color',cmap_lines(2,:),'markersize',markersize/3)
% plot(X3,Y3,'.','color',cmap_lines(5,:),'markersize',.75*markersize)
% plot(X4,Y4,'k.','markersize',markersize2)
ylabel(['pH (total)'])
xlabel(['TA-DIC \mumol kg^{-1}'])
xlim([0 500])
ylim([7.4 9])
% legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','Jiang 1750','Jiang 2010','Jiang 2100 (SSP585)','location','southeast')
legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','location','southeast')
title(['pH vs TA-DIC'])
hold off

% - - -
k=.8;
set(gcf, 'PaperPosition', [0 0 k*6 k*6])   
print(gcf, [pwd '/png/regression_global_phtot_vs_alkstar_0-10m_GLODAP_v21.png'], '-dpng', '-r600' );  


% - - -
% - - -
% - - -
% - - -
% - - -
% - - -
 
% polar - Regression of GLODAP observation data vs TA-DIC

% - - -
% - - -
% - - -
% - - -
% - - -
% - - -




% - - -
% - - -
% - - -

% polar - Regression of CO3^2- vs TA-DIC

% - - -
% - - -
% - - -

deplim = 10;
markersize = 5;
markersize2 = 5;
fontsize_title = 16;
fontsize_lab = 12;
fontsize_ax = 12;
fontsize_eqn = 10;
fontsize_stats = 10;

figure(1)
hFig=gcf;
clf(hFig);

tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

hax1 = nexttile;

% Jiang data
Y1 = co3_1750;
Y1(~in_polar) = nan;
Y2 = co3_2010;
Y2(~in_polar) = nan;
Y3 = co3_2100;
Y3(~in_polar) = nan;

X1 = talk_1750 -  dic_1750;
X1(~in_polar) = nan;
X2 = talk_2010 - dic_2010;
X2(~in_polar) = nan;
X3 = talk_2100 - dic_2100;
X3(~in_polar) = nan;

X1 = reshape(X1,[],1);
X2 = reshape(X2,[],1);
X3 = reshape(X3,[],1);
Y1 = reshape(Y1,[],1);
Y2 = reshape(Y2,[],1);
Y3 = reshape(Y3,[],1);

% GLODAP obs
idx = find(G2.in_polar & G2.depth<deplim);
Y4 = G2.co2sys_co3(idx);
X4 = G2.talk(idx) - G2.dic(idx);

	% regression with 95%CI
	x = X4;
	y = Y4;
	% eliminate nans from x and y
	idx = find(~isnan(x) & ~isnan(y));
	x = x(idx);
	y = y(idx);
	% create a function for the nonlinear regression model equation
	f = @(param,xval) param(1) + param(2) .*  xval + param(3) .* xval .^2 + param(4) .* xval .^ 3;
	% specify initial parameter values to use for the nonlinear regression optimization
	p_init = [nanmean2(y(:)) 0 0 0];
	% use nlinfit to find the optimum parameter values (popt) and covariance matrix (pcov)
	[popt,~,~,pcov] = nlinfit(x,y,f,p_init);					  
	% evenly spaced new observations that extend from below min(x) to above max(x)
	nlinspace = 100;
	x_new = linspace(nanmin2(x(:)), nanmax2(x(:)), nlinspace);
	% probability level for the prediction limits (e.g. alpha=0.05 for 95% prediction limits and 95% confidence limits)
	alpha=0.05;
	% use the delta_method.m function to calculate confidence intervals and prediction intervals
	d = delta_method(pcov,popt,x_new,f,x,y,alpha);

hold on

	patch([d.x_new' fliplr(d.x_new')],[d.upr_pred' fliplr(d.lwr_pred')],[.8 .8 .8],'edgecolor','none','facealpha',0.5)
	plot(x,y,'k.','markersize',5)
	plot(d.x_new,d.y_new,'k-','linewidth',1)
	b1 = num2str(popt(1),'%.2s');
	b2 = num2str(popt(2),'%.2s');
	b3 = num2str(popt(3),'%.2s');
	b4 = num2str(popt(4),'%.2s');
	eqn = ['y = ',b1,' + ',b2,' x + ',b3,' x^2 + ',b4,' x^3'];
	pstr = num2str(d.pvalue);
	rsq = num2str(d.rsquared,"%.3f");
	ser = num2str(d.syx,"%.1f");
	nobs = num2str(d.nobs);
	text(25,370,'GLODAP regression statistics','fontsize',fontsize_eqn);    
	text(25,350,eqn,'fontsize',fontsize_eqn);    
	text(25,325,['r^2=',rsq,', p=',pstr,', rmse=',ser,', n=',nobs],'fontsize',fontsize_stats);     

% plot(X1,Y1,'.','color',cmap_lines(1,:),'markersize',markersize)
% plot(X2,Y2,'o','color',cmap_lines(2,:),'markersize',markersize/3)
% plot(X3,Y3,'.','color',cmap_lines(5,:),'markersize',.75*markersize)
% plot(X4,Y4,'k.','markersize',markersize2)
ylabel(['CO_3^{2-} \mumol kg^{-1}'])
xlabel(['TA-DIC \mumol kg^{-1}'])
xlim([0 500])
ylim([0 400])
% legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','Jiang 1750','Jiang 2010','Jiang 2100 (SSP585)','location','southeast')
legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','location','southeast')
title(['CO_3^{2-} vs TA-DIC'])
hold off

% - - -
k=.8;
set(gcf, 'PaperPosition', [0 0 k*6 k*6])   
print(gcf, [pwd '/png/regression_polar_co3_vs_alkstar_0-10m_GLODAP_v21.png'], '-dpng', '-r600' );  



% - - -
% - - -
% - - -

% polar - Regression of omara vs TA-DIC

% - - -
% - - -
% - - -

deplim = 10;
markersize = 5;
markersize2 = 5;
fontsize_title = 16;
fontsize_lab = 12;
fontsize_ax = 12;
fontsize_eqn = 10;
fontsize_stats = 10;

figure(1)
hFig=gcf;
clf(hFig);

tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

hax1 = nexttile;

% Jiang data
Y1 = omara_1750;
Y1(~in_polar) = nan;
Y2 = omara_2010;
Y2(~in_polar) = nan;
Y3 = omara_2100;
Y3(~in_polar) = nan;

X1 = talk_1750 -  dic_1750;
X1(~in_polar) = nan;
X2 = talk_2010 - dic_2010;
X2(~in_polar) = nan;
X3 = talk_2100 - dic_2100;
X3(~in_polar) = nan;

X1 = reshape(X1,[],1);
X2 = reshape(X2,[],1);
X3 = reshape(X3,[],1);
Y1 = reshape(Y1,[],1);
Y2 = reshape(Y2,[],1);
Y3 = reshape(Y3,[],1);

% GLODAP obs
idx = find(G2.in_polar & G2.depth<deplim);
Y4 = G2.co2sys_omara(idx);
X4 = G2.talk(idx) - G2.dic(idx);

	% regression with 95%CI
	x = X4;
	y = Y4;
	% eliminate nans from x and y
	idx = find(~isnan(x) & ~isnan(y));
	x = x(idx);
	y = y(idx);
	% create a function for the nonlinear regression model equation
	f = @(param,xval) param(1) + param(2) .*  xval + param(3) .* xval .^2 + param(4) .* xval .^ 3;
	% specify initial parameter values to use for the nonlinear regression optimization
	p_init = [nanmean2(y(:)) 0 0 0];
	% use nlinfit to find the optimum parameter values (popt) and covariance matrix (pcov)
	[popt,~,~,pcov] = nlinfit(x,y,f,p_init);					  
	% evenly spaced new observations that extend from below min(x) to above max(x)
	nlinspace = 100;
	x_new = linspace(nanmin2(x(:)), nanmax2(x(:)), nlinspace);
	% probability level for the prediction limits (e.g. alpha=0.05 for 95% prediction limits and 95% confidence limits)
	alpha=0.05;
	% use the delta_method.m function to calculate confidence intervals and prediction intervals
	d = delta_method(pcov,popt,x_new,f,x,y,alpha);

hold on

	patch([d.x_new' fliplr(d.x_new')],[d.upr_pred' fliplr(d.lwr_pred')],[.8 .8 .8],'edgecolor','none','facealpha',0.5)
	plot(x,y,'k.','markersize',5)
	plot(d.x_new,d.y_new,'k-','linewidth',1)
	b1 = num2str(popt(1),'%.2s');
	b2 = num2str(popt(2),'%.2s');
	b3 = num2str(popt(3),'%.2s');
	b4 = num2str(popt(4),'%.2s');
	eqn = ['y = ',b1,' + ',b2,' x + ',b3,' x^2 + ',b4,' x^3'];
	pstr = num2str(d.pvalue);
	rsq = num2str(d.rsquared,"%.3f");
	ser = num2str(d.syx,"%.2f");
	nobs = num2str(d.nobs);
	text(25,6,'GLODAP regression statistics','fontsize',fontsize_eqn);    
	text(25,5.6,eqn,'fontsize',fontsize_eqn);    
	text(25,5.2,['r^2=',rsq,', p=',pstr,', rmse=',ser,', n=',nobs],'fontsize',fontsize_stats);     

% plot(X1,Y1,'.','color',cmap_lines(1,:),'markersize',markersize)
% plot(X2,Y2,'o','color',cmap_lines(2,:),'markersize',markersize/3)
% plot(X3,Y3,'.','color',cmap_lines(5,:),'markersize',.75*markersize)
% plot(X4,Y4,'k.','markersize',markersize2)
ylabel(['\Omega_{ara}'])
xlabel(['TA-DIC \mumol kg^{-1}'])
xlim([0 500])
ylim([0 7])
% legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','Jiang 1750','Jiang 2010','Jiang 2100 (SSP585)','location','southeast')
legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','location','southeast')
title(['\Omega_{ara} vs TA-DIC'])
hold off

% - - -
k=.8;
set(gcf, 'PaperPosition', [0 0 k*6 k*6])   
print(gcf, [pwd '/png/regression_polar_omara_vs_alkstar_0-10m_GLODAP_v21.png'], '-dpng', '-r600' );  



% - - -
% - - -
% - - -

% polar - Regression of phtot vs TA-DIC

% - - -
% - - -
% - - -

deplim = 10;
markersize = 5;
markersize2 = 5;
fontsize_title = 16;
fontsize_lab = 12;
fontsize_ax = 12;
fontsize_eqn = 10;
fontsize_stats = 10;

figure(1)
hFig=gcf;
clf(hFig);

tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

hax1 = nexttile;

% Jiang data
Y1 = phtot_1750;
Y1(~in_polar) = nan;
Y2 = phtot_2010;
Y2(~in_polar) = nan;
Y3 = phtot_2100;
Y3(~in_polar) = nan;

X1 = talk_1750 -  dic_1750;
X1(~in_polar) = nan;
X2 = talk_2010 - dic_2010;
X2(~in_polar) = nan;
X3 = talk_2100 - dic_2100;
X3(~in_polar) = nan;

X1 = reshape(X1,[],1);
X2 = reshape(X2,[],1);
X3 = reshape(X3,[],1);
Y1 = reshape(Y1,[],1);
Y2 = reshape(Y2,[],1);
Y3 = reshape(Y3,[],1);

% GLODAP obs
idx = find(G2.in_polar & G2.depth<deplim);
Y4 = G2.co2sys_phtot(idx);
X4 = G2.talk(idx) - G2.dic(idx);

	% regression with 95%CI
	x = X4;
	y = Y4;
	% eliminate nans from x and y
	idx = find(~isnan(x) & ~isnan(y));
	x = x(idx);
	y = y(idx);
	% create a function for the nonlinear regression model equation
	f = @(param,xval) param(1) + param(2) .*  xval + param(3) .* xval .^2 + param(4) .* xval .^ 3;
	% specify initial parameter values to use for the nonlinear regression optimization
	p_init = [nanmean2(y(:)) 0 0 0];
	% use nlinfit to find the optimum parameter values (popt) and covariance matrix (pcov)
	[popt,~,~,pcov] = nlinfit(x,y,f,p_init);					  
	% evenly spaced new observations that extend from below min(x) to above max(x)
	nlinspace = 100;
	x_new = linspace(nanmin2(x(:)), nanmax2(x(:)), nlinspace);
	% probability level for the prediction limits (e.g. alpha=0.05 for 95% prediction limits and 95% confidence limits)
	alpha=0.05;
	% use the delta_method.m function to calculate confidence intervals and prediction intervals
	d = delta_method(pcov,popt,x_new,f,x,y,alpha);

hold on

	patch([d.x_new' fliplr(d.x_new')],[d.upr_pred' fliplr(d.lwr_pred')],[.8 .8 .8],'edgecolor','none','facealpha',0.5)
	plot(x,y,'k.','markersize',5)
	plot(d.x_new,d.y_new,'k-','linewidth',1)
	b1 = num2str(popt(1),'%.2s');
	b2 = num2str(popt(2),'%.2s');
	b3 = num2str(popt(3),'%.2s');
	b4 = num2str(popt(4),'%.2s');
	eqn = ['y = ',b1,' + ',b2,' x + ',b3,' x^2 + ',b4,' x^3'];
	pstr = num2str(d.pvalue);
	rsq = num2str(d.rsquared,"%.3f");
	ser = num2str(d.syx,"%.2f");
	nobs = num2str(d.nobs);
	text(25,8.85,'GLODAP regression statistics','fontsize',fontsize_eqn);    
	text(25,8.75,eqn,'fontsize',fontsize_eqn);    
	text(25,8.65,['r^2=',rsq,', p=',pstr,', rmse=',ser,', n=',nobs],'fontsize',fontsize_stats);     

% plot(X1,Y1,'.','color',cmap_lines(1,:),'markersize',markersize)
% plot(X2,Y2,'o','color',cmap_lines(2,:),'markersize',markersize/3)
% plot(X3,Y3,'.','color',cmap_lines(5,:),'markersize',.75*markersize)
% plot(X4,Y4,'k.','markersize',markersize2)
ylabel(['pH (total)'])
xlabel(['TA-DIC \mumol kg^{-1}'])
xlim([0 500])
ylim([7.4 9])
% legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','Jiang 1750','Jiang 2010','Jiang 2100 (SSP585)','location','southeast')
legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','location','southeast')
title(['pH vs TA-DIC'])
hold off

% - - -
k=.8;
set(gcf, 'PaperPosition', [0 0 k*6 k*6])   
print(gcf, [pwd '/png/regression_polar_phtot_vs_alkstar_0-10m_GLODAP_v21.png'], '-dpng', '-r600' );  














% - - -
% - - -
% - - -
% - - -
% - - -
% - - -
 
% dist300 - Regression of GLODAP observation data vs TA-DIC

% - - -
% - - -
% - - -
% - - -
% - - -
% - - -




% - - -
% - - -
% - - -

% dist300 - Regression of CO3^2- vs TA-DIC

% - - -
% - - -
% - - -

deplim = 10;
markersize = 5;
markersize2 = 5;
fontsize_title = 16;
fontsize_lab = 12;
fontsize_ax = 12;
fontsize_eqn = 10;
fontsize_stats = 10;

figure(1)
hFig=gcf;
clf(hFig);

tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

hax1 = nexttile;

% Jiang data
Y1 = co3_1750;
Y1(~in_dist300) = nan;
Y2 = co3_2010;
Y2(~in_dist300) = nan;
Y3 = co3_2100;
Y3(~in_dist300) = nan;

X1 = talk_1750 -  dic_1750;
X1(~in_dist300) = nan;
X2 = talk_2010 - dic_2010;
X2(~in_dist300) = nan;
X3 = talk_2100 - dic_2100;
X3(~in_dist300) = nan;

X1 = reshape(X1,[],1);
X2 = reshape(X2,[],1);
X3 = reshape(X3,[],1);
Y1 = reshape(Y1,[],1);
Y2 = reshape(Y2,[],1);
Y3 = reshape(Y3,[],1);

% GLODAP obs
idx = find(G2.in_dist300 & G2.depth<deplim);
Y4 = G2.co2sys_co3(idx);
X4 = G2.talk(idx) - G2.dic(idx);

	% regression with 95%CI
	x = X4;
	y = Y4;
	% eliminate nans from x and y
	idx = find(~isnan(x) & ~isnan(y));
	x = x(idx);
	y = y(idx);
	% create a function for the nonlinear regression model equation
	f = @(param,xval) param(1) + param(2) .*  xval + param(3) .* xval .^2 + param(4) .* xval .^ 3;
	% specify initial parameter values to use for the nonlinear regression optimization
	p_init = [nanmean2(y(:)) 0 0 0];
	% use nlinfit to find the optimum parameter values (popt) and covariance matrix (pcov)
	[popt,~,~,pcov] = nlinfit(x,y,f,p_init);					  
	% evenly spaced new observations that extend from below min(x) to above max(x)
	nlinspace = 100;
	x_new = linspace(nanmin2(x(:)), nanmax2(x(:)), nlinspace);
	% probability level for the prediction limits (e.g. alpha=0.05 for 95% prediction limits and 95% confidence limits)
	alpha=0.05;
	% use the delta_method.m function to calculate confidence intervals and prediction intervals
	d = delta_method(pcov,popt,x_new,f,x,y,alpha);

hold on

	patch([d.x_new' fliplr(d.x_new')],[d.upr_pred' fliplr(d.lwr_pred')],[.8 .8 .8],'edgecolor','none','facealpha',0.5)
	plot(x,y,'k.','markersize',5)
	plot(d.x_new,d.y_new,'k-','linewidth',1)
	b1 = num2str(popt(1),'%.2s');
	b2 = num2str(popt(2),'%.2s');
	b3 = num2str(popt(3),'%.2s');
	b4 = num2str(popt(4),'%.2s');
	eqn = ['y = ',b1,' + ',b2,' x + ',b3,' x^2 + ',b4,' x^3'];
	pstr = num2str(d.pvalue);
	rsq = num2str(d.rsquared,"%.3f");
	ser = num2str(d.syx,"%.1f");
	nobs = num2str(d.nobs);
	text(25,370,'GLODAP regression statistics','fontsize',fontsize_eqn);    
	text(25,350,eqn,'fontsize',fontsize_eqn);    
	text(25,325,['r^2=',rsq,', p=',pstr,', rmse=',ser,', n=',nobs],'fontsize',fontsize_stats);     

% plot(X1,Y1,'.','color',cmap_lines(1,:),'markersize',markersize)
% plot(X2,Y2,'o','color',cmap_lines(2,:),'markersize',markersize/3)
% plot(X3,Y3,'.','color',cmap_lines(5,:),'markersize',.75*markersize)
% plot(X4,Y4,'k.','markersize',markersize2)
ylabel(['CO_3^{2-} \mumol kg^{-1}'])
xlabel(['TA-DIC \mumol kg^{-1}'])
xlim([0 500])
ylim([0 400])
% legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','Jiang 1750','Jiang 2010','Jiang 2100 (SSP585)','location','southeast')
legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','location','southeast')
title(['CO_3^{2-} vs TA-DIC'])
hold off

% - - -
k=.8;
set(gcf, 'PaperPosition', [0 0 k*6 k*6])   
print(gcf, [pwd '/png/regression_dist300_co3_vs_alkstar_0-10m_GLODAP_v21.png'], '-dpng', '-r600' );  



% - - -
% - - -
% - - -

% dist300 - Regression of omara vs TA-DIC

% - - -
% - - -
% - - -

deplim = 10;
markersize = 5;
markersize2 = 5;
fontsize_title = 16;
fontsize_lab = 12;
fontsize_ax = 12;
fontsize_eqn = 10;
fontsize_stats = 10;

figure(1)
hFig=gcf;
clf(hFig);

tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

hax1 = nexttile;

% Jiang data
Y1 = omara_1750;
Y1(~in_dist300) = nan;
Y2 = omara_2010;
Y2(~in_dist300) = nan;
Y3 = omara_2100;
Y3(~in_dist300) = nan;

X1 = talk_1750 -  dic_1750;
X1(~in_dist300) = nan;
X2 = talk_2010 - dic_2010;
X2(~in_dist300) = nan;
X3 = talk_2100 - dic_2100;
X3(~in_dist300) = nan;

X1 = reshape(X1,[],1);
X2 = reshape(X2,[],1);
X3 = reshape(X3,[],1);
Y1 = reshape(Y1,[],1);
Y2 = reshape(Y2,[],1);
Y3 = reshape(Y3,[],1);

% GLODAP obs
idx = find(G2.in_dist300 & G2.depth<deplim);
Y4 = G2.co2sys_omara(idx);
X4 = G2.talk(idx) - G2.dic(idx);

	% regression with 95%CI
	x = X4;
	y = Y4;
	% eliminate nans from x and y
	idx = find(~isnan(x) & ~isnan(y));
	x = x(idx);
	y = y(idx);
	% create a function for the nonlinear regression model equation
	f = @(param,xval) param(1) + param(2) .*  xval + param(3) .* xval .^2 + param(4) .* xval .^ 3;
	% specify initial parameter values to use for the nonlinear regression optimization
	p_init = [nanmean2(y(:)) 0 0 0];
	% use nlinfit to find the optimum parameter values (popt) and covariance matrix (pcov)
	[popt,~,~,pcov] = nlinfit(x,y,f,p_init);					  
	% evenly spaced new observations that extend from below min(x) to above max(x)
	nlinspace = 100;
	x_new = linspace(nanmin2(x(:)), nanmax2(x(:)), nlinspace);
	% probability level for the prediction limits (e.g. alpha=0.05 for 95% prediction limits and 95% confidence limits)
	alpha=0.05;
	% use the delta_method.m function to calculate confidence intervals and prediction intervals
	d = delta_method(pcov,popt,x_new,f,x,y,alpha);

hold on

	patch([d.x_new' fliplr(d.x_new')],[d.upr_pred' fliplr(d.lwr_pred')],[.8 .8 .8],'edgecolor','none','facealpha',0.5)
	plot(x,y,'k.','markersize',5)
	plot(d.x_new,d.y_new,'k-','linewidth',1)
	b1 = num2str(popt(1),'%.2s');
	b2 = num2str(popt(2),'%.2s');
	b3 = num2str(popt(3),'%.2s');
	b4 = num2str(popt(4),'%.2s');
	eqn = ['y = ',b1,' + ',b2,' x + ',b3,' x^2 + ',b4,' x^3'];
	pstr = num2str(d.pvalue);
	rsq = num2str(d.rsquared,"%.3f");
	ser = num2str(d.syx,"%.2f");
	nobs = num2str(d.nobs);
	text(25,6,'GLODAP regression statistics','fontsize',fontsize_eqn);    
	text(25,5.6,eqn,'fontsize',fontsize_eqn);    
	text(25,5.2,['r^2=',rsq,', p=',pstr,', rmse=',ser,', n=',nobs],'fontsize',fontsize_stats);     

% plot(X1,Y1,'.','color',cmap_lines(1,:),'markersize',markersize)
% plot(X2,Y2,'o','color',cmap_lines(2,:),'markersize',markersize/3)
% plot(X3,Y3,'.','color',cmap_lines(5,:),'markersize',.75*markersize)
% plot(X4,Y4,'k.','markersize',markersize2)
ylabel(['\Omega_{ara}'])
xlabel(['TA-DIC \mumol kg^{-1}'])
xlim([0 500])
ylim([0 7])
% legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','Jiang 1750','Jiang 2010','Jiang 2100 (SSP585)','location','southeast')
legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','location','southeast')
title(['\Omega_{ara} vs TA-DIC'])
hold off

% - - -
k=.8;
set(gcf, 'PaperPosition', [0 0 k*6 k*6])   
print(gcf, [pwd '/png/regression_dist300_omara_vs_alkstar_0-10m_GLODAP_v21.png'], '-dpng', '-r600' );  



% - - -
% - - -
% - - -

% dist300 - Regression of phtot vs TA-DIC

% - - -
% - - -
% - - -

deplim = 10;
markersize = 5;
markersize2 = 5;
fontsize_title = 16;
fontsize_lab = 12;
fontsize_ax = 12;
fontsize_eqn = 10;
fontsize_stats = 10;

figure(1)
hFig=gcf;
clf(hFig);

tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

hax1 = nexttile;

% Jiang data
Y1 = phtot_1750;
Y1(~in_dist300) = nan;
Y2 = phtot_2010;
Y2(~in_dist300) = nan;
Y3 = phtot_2100;
Y3(~in_dist300) = nan;

X1 = talk_1750 -  dic_1750;
X1(~in_dist300) = nan;
X2 = talk_2010 - dic_2010;
X2(~in_dist300) = nan;
X3 = talk_2100 - dic_2100;
X3(~in_dist300) = nan;

X1 = reshape(X1,[],1);
X2 = reshape(X2,[],1);
X3 = reshape(X3,[],1);
Y1 = reshape(Y1,[],1);
Y2 = reshape(Y2,[],1);
Y3 = reshape(Y3,[],1);

% GLODAP obs
idx = find(G2.in_dist300 & G2.depth<deplim);
Y4 = G2.co2sys_phtot(idx);
X4 = G2.talk(idx) - G2.dic(idx);

	% regression with 95%CI
	x = X4;
	y = Y4;
	% eliminate nans from x and y
	idx = find(~isnan(x) & ~isnan(y));
	x = x(idx);
	y = y(idx);
	% create a function for the nonlinear regression model equation
	f = @(param,xval) param(1) + param(2) .*  xval + param(3) .* xval .^2 + param(4) .* xval .^ 3;
	% specify initial parameter values to use for the nonlinear regression optimization
	p_init = [nanmean2(y(:)) 0 0 0];
	% use nlinfit to find the optimum parameter values (popt) and covariance matrix (pcov)
	[popt,~,~,pcov] = nlinfit(x,y,f,p_init);					  
	% evenly spaced new observations that extend from below min(x) to above max(x)
	nlinspace = 100;
	x_new = linspace(nanmin2(x(:)), nanmax2(x(:)), nlinspace);
	% probability level for the prediction limits (e.g. alpha=0.05 for 95% prediction limits and 95% confidence limits)
	alpha=0.05;
	% use the delta_method.m function to calculate confidence intervals and prediction intervals
	d = delta_method(pcov,popt,x_new,f,x,y,alpha);

hold on

	patch([d.x_new' fliplr(d.x_new')],[d.upr_pred' fliplr(d.lwr_pred')],[.8 .8 .8],'edgecolor','none','facealpha',0.5)
	plot(x,y,'k.','markersize',5)
	plot(d.x_new,d.y_new,'k-','linewidth',1)
	b1 = num2str(popt(1),'%.2s');
	b2 = num2str(popt(2),'%.2s');
	b3 = num2str(popt(3),'%.2s');
	b4 = num2str(popt(4),'%.2s');
	eqn = ['y = ',b1,' + ',b2,' x + ',b3,' x^2 + ',b4,' x^3'];
	pstr = num2str(d.pvalue);
	rsq = num2str(d.rsquared,"%.3f");
	ser = num2str(d.syx,"%.2f");
	nobs = num2str(d.nobs);
	text(25,8.85,'GLODAP regression statistics','fontsize',fontsize_eqn);    
	text(25,8.75,eqn,'fontsize',fontsize_eqn);    
	text(25,8.65,['r^2=',rsq,', p=',pstr,', rmse=',ser,', n=',nobs],'fontsize',fontsize_stats);     

% plot(X1,Y1,'.','color',cmap_lines(1,:),'markersize',markersize)
% plot(X2,Y2,'o','color',cmap_lines(2,:),'markersize',markersize/3)
% plot(X3,Y3,'.','color',cmap_lines(5,:),'markersize',.75*markersize)
% plot(X4,Y4,'k.','markersize',markersize2)
ylabel(['pH (total)'])
xlabel(['TA-DIC \mumol kg^{-1}'])
xlim([0 500])
ylim([7.4 9])
% legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','Jiang 1750','Jiang 2010','Jiang 2100 (SSP585)','location','southeast')
legend('GLODAP 95%PI','GLODAP obs','GLODAP regression','location','southeast')
title(['pH vs TA-DIC'])
hold off

% - - -
k=.8;
set(gcf, 'PaperPosition', [0 0 k*6 k*6])   
print(gcf, [pwd '/png/regression_dist300_phtot_vs_alkstar_0-10m_GLODAP_v21.png'], '-dpng', '-r600' );  












