function Annual_record = umodel_core (Path,Global_temperature_setting,Soil_respiration_equation,C_soil_option)



%% Data input
%NPP
% P model 10 years mean, provided by Shirley
NPP_ref=ncread([Path,'NPP\NPP_from_p_model_10years_360720.nc'],'npp');% (KgC/m2)
        NPP_ref=NPP_ref/365/24;
% CO2 record
Atm_Co2=csvread([Path,'CO2\CO2increase.csv']);
% u equation
Location_of_u=[Path,'NPP\half_saturation_u_360720.nc'];
u_value=ncread(Location_of_u,'u');
% Cveg in 1993
load vegetation\Totalbiomass_1993_360720.mat;
Cveg=Total_biomass_360720/10000*(10^3); %from MgC/ha to KgC/m2
%Csoil at the momment
switch C_soil_option
    case 'Assume_1500PgC'
     LT.Csoil=1500; %PgC   
    case 'Soilgrids'
        warning off
        %[LT.Csoil,~]=geotiffread('Soil\Csoil360720.tif');
        load Soil\Csoil_360720_gapfill.mat
        LT.Csoil=Csoil;
        LT.Csoil=LT.Csoil*1000/10000; %from Ton/ha to KgC/m2
        warning on
end

% Temperature record
 Temperature_annually=load([Path,'CRU_temperature\cru_temperature_1901_2017_05degree_annually.mat']);%Near surface Temp, 'degrees Celsius'
Temperature_annually=Temperature_annually.Temperature_annually(:,:,3:end);
%Temperature_annually(Temperature_annually<0)=0;
Temperature=Temperature_annually(:,:,1);
%Temperature(Temperature<3)=3;
% clear pixel that NaN in any dataset
Uniform_missing_pixel=LT.Csoil+Cveg+u_value+NPP_ref+Temperature;
Uniform_missing_pixel(~isnan(Uniform_missing_pixel))=0;
for year=1:size(Temperature_annually,3)
    if strcmp(Global_temperature_setting,'constant_temperature')
            Temperature_annually(:,:,year)=Temperature_annually(:,:,1)+Uniform_missing_pixel;
    else
    Temperature_annually(:,:,year)=Temperature_annually(:,:,year)+Uniform_missing_pixel;

    end
end
%NPP_ref=NPP_ref+Uniform_missing_pixel;
Cveg=Cveg+Uniform_missing_pixel;
u_value=u_value+Uniform_missing_pixel;
Temperature=Temperature_annually(:,:,1);
% Other Coefficient
co2_ref=372;
u=@(x) (x./(x+u_value)).*((co2_ref+u_value)./co2_ref).*NPP_ref;%u function
LT.E0=308.56;
LT.Tref=15;
LT.T0=46.02; %unit oC
switch Soil_respiration_equation
    case 'LT'
LT.func=@(x) nthroot(exp(LT.E0.*(1./(LT.Tref+LT.T0)-1./(x+LT.T0))),3);
    case 'Q10'
LT.func=@(x) 1.4.^((x-LT.Tref)/10);
end
% Map of area per grid
load('Area_WGS_1984_05degree.mat') %this is a map of area per quadrangle collapse 




%% vegetation box: calculate vegetation turnover rate
CO2_1993=Atm_Co2(Atm_Co2(:,1)==1993,2);
NPP_1993=u(CO2_1993);
tau_veg=Cveg./(24*365*NPP_1993); %unit: year
tau_veg(tau_veg<0.5)=0.5;% get rid of some extremely small value
%tau_veg(tau_veg<0.5)=0.5;% anything smalled than 0.5 would cause problem (when timestep is a year)
clearvars Location_of_u Location_of_Cveg CO2_1901 CO2_1993 Cveg Csoil

%% Soil Box: calculate K2
CO2_1903=Atm_Co2(Atm_Co2(:,1)==1903,2);% calculate total NPP in 1901
NPP_1903=u(CO2_1903)*24*365; %unit KgC/m^2

switch C_soil_option
    case 'Soilgrids'
    %PgC
    K_2 =NPP_1903./LT.Csoil./LT.func(Temperature);%note that K_2 and LT.func dont have unit
C_soil=NPP_1903./K_2./LT.func(Temperature);%KgC/m2
    otherwise
    K_2 =sum(sum(Area_WGS_1984_05degree.*NPP_1903./LT.func(Temperature),'omitnan'),'omitnan' )./(LT.Csoil*10^12);%note that K_2 and LT.func dont have unit
    C_soil=NPP_1903./K_2./LT.func(Temperature);%KgC/m2
end
%% model sign up
C_veg=tau_veg.*(NPP_1903);
Start_year=1903;
End_year=2017;

%% model run

for year=Start_year:End_year
CO2=Atm_Co2(Atm_Co2(:,1)==year,2);
NPP=u(CO2)*24*365; %unit KgC/m^2
Flux_veg_to_soil=C_veg./tau_veg;%unit KgC/m^2
C_veg_new=C_veg+NPP-Flux_veg_to_soil;   
Soil_Respiration=K_2.*C_soil.*LT.func(Temperature_annually(:,:,year-Start_year+1)); %unit KgC/m^2    
C_soil_new=C_soil+Flux_veg_to_soil-Soil_Respiration;
Land_sink=NPP-Soil_Respiration; 
C_veg=C_veg_new;
C_soil=C_soil_new; 
Annual_record.midflux_map(:,:,year-Start_year+1)=Flux_veg_to_soil;
Annual_record.NPP_map(:,:,year-Start_year+1)=NPP;
Annual_record.Respiration_map(:,:,year-Start_year+1)=Soil_Respiration;
Annual_record.Csoil_map(:,:,year-Start_year+1)=C_soil;
Annual_record.C_veg_map(:,:,year-Start_year+1)=C_veg;
Annual_record.Land_sinkmap(:,:,year-Start_year+1)=Land_sink;
Annual_record.Flux_veg_to_soil(year-Start_year+1)=sum(sum(Area_WGS_1984_05degree.*Flux_veg_to_soil/(10^12),'omitnan'),'omitnan');
Annual_record.Soil_Respiration(year-Start_year+1)=sum(sum(Area_WGS_1984_05degree.*Soil_Respiration/(10^12),'omitnan'),'omitnan');
Annual_record.Land_sink(year-Start_year+1)=sum(sum(Area_WGS_1984_05degree.*Land_sink/(10^12),'omitnan'),'omitnan');
Annual_record.NPP(year-Start_year+1)=sum(sum(Area_WGS_1984_05degree.*NPP/(10^12),'omitnan'),'omitnan');
Annual_record.C_soil(year-Start_year+1)=sum(sum(Area_WGS_1984_05degree.*C_soil/(10^12),'omitnan'),'omitnan');
Annual_record.C_veg(year-Start_year+1)=sum(sum(Area_WGS_1984_05degree.*C_veg/(10^12),'omitnan'),'omitnan');

end
Annual_record.K_2=K_2;
%% 6.Draw something...


