clc
clear
close all

%% run the model
Path = 'J:\Imperial_college\Master_Project\New_model\u_model_masterfile.mlx'; %the directory of the script
Path = [fileparts(Path),filesep'];%get the directory of the current script
cd(Path)
C_soil_option='Assume_1500PgC';
Soil_respiration_equation='Q10';
Annual_record{1} = umodel_core (Path,'constant_temperature',Soil_respiration_equation,C_soil_option); %this is to trigger the model with different setting.
Annual_record{2} = umodel_core (Path,'yearly_temperature',Soil_respiration_equation,C_soil_option);
load('Area_WGS_1984_05degree.mat') %this is a map of area per quadrangle collapse 

%% compare spatial pattern with inversion models
Threshold_latitude=@(x) 360-(x*2+180);
Land_sink=Area_WGS_1984_05degree.*Annual_record{1, 1}.Land_sinkmap;
Land_sink1=Land_sink(:,:,93:108);
Land_sink1_global=sum(sum(Land_sink1(:,:,:)/(10^9),'omitnan'),'omitnan');
Land_sink1_global=reshape(Land_sink1_global,[],1);
Land_sink2(:,1)=Land_sink1_global;

Land_sink1_55N=sum(sum(Land_sink1(1:Threshold_latitude(55),:,:)/(10^9),'omitnan'),'omitnan');
Land_sink1_55N=reshape(Land_sink1_55N,[],1);
Land_sink2(:,2)=Land_sink1_55N;

Land_sink1_35N55N=sum(sum(Land_sink1(Threshold_latitude(55):Threshold_latitude(35),:,:)/(10^9),'omitnan'),'omitnan');
Land_sink1_35N55N=reshape(Land_sink1_35N55N,[],1);
Land_sink2(:,3)=Land_sink1_35N55N;

Land_sink1_15N35N=sum(sum(Land_sink1(Threshold_latitude(35):Threshold_latitude(15),:,:)/(10^9),'omitnan'),'omitnan');
Land_sink1_15N35N=reshape(Land_sink1_15N35N,[],1);
Land_sink2(:,4)=Land_sink1_15N35N;

Land_sink1_15N15S=sum(sum(Land_sink1(Threshold_latitude(15):Threshold_latitude(-15),:,:)/(10^9),'omitnan'),'omitnan');
Land_sink1_15N15S=reshape(Land_sink1_15N15S,[],1);
Land_sink2(:,5)=Land_sink1_15N15S;

Land_sink1_15S35S=sum(sum(Land_sink1(Threshold_latitude(-15):Threshold_latitude(-35),:,:)/(10^9),'omitnan'),'omitnan');
Land_sink1_15S35S=reshape(Land_sink1_15S35S,[],1);
Land_sink2(:,6)=Land_sink1_15S35S;

Land_sink1_35S55S=sum(sum(Land_sink1(Threshold_latitude(-35):Threshold_latitude(-55),:,:)/(10^9),'omitnan'),'omitnan');
Land_sink1_35S55S=reshape(Land_sink1_35S55S,[],1);
Land_sink2(:,7)=Land_sink1_35S55S;
%% map of land use change
LandUseChange=fliplr(rot90(LandUseChange,3));
Threshold_latitude=@(x) 360-(x*2+180);
Land_sink1=Area_WGS_1984_05degree.*LandUseChange(:,:);
Land_sink1_global=sum(sum(Land_sink1(:,:)/(10^12),'omitnan'),'omitnan');
Land_sink1_global=reshape(Land_sink1_global,[],1);
Land_sink2(:,1)=Land_sink1_global;

Land_sink1_55N=sum(sum(Land_sink1(1:Threshold_latitude(55),:)/(10^12),'omitnan'),'omitnan');
Land_sink1_55N=reshape(Land_sink1_55N,[],1);
Land_sink2(:,2)=Land_sink1_55N;

Land_sink1_35N55N=sum(sum(Land_sink1(Threshold_latitude(55):Threshold_latitude(35),:)/(10^12),'omitnan'),'omitnan');
Land_sink1_35N55N=reshape(Land_sink1_35N55N,[],1);
Land_sink2(:,3)=Land_sink1_35N55N;

Land_sink1_15N35N=sum(sum(Land_sink1(Threshold_latitude(35):Threshold_latitude(15),:)/(10^12),'omitnan'),'omitnan');
Land_sink1_15N35N=reshape(Land_sink1_15N35N,[],1);
Land_sink2(:,4)=Land_sink1_15N35N;

Land_sink1_15N15S=sum(sum(Land_sink1(Threshold_latitude(15):Threshold_latitude(-15),:)/(10^12),'omitnan'),'omitnan');
Land_sink1_15N15S=reshape(Land_sink1_15N15S,[],1);
Land_sink2(:,5)=Land_sink1_15N15S;

Land_sink1_15S35S=sum(sum(Land_sink1(Threshold_latitude(-15):Threshold_latitude(-35),:)/(10^12),'omitnan'),'omitnan');
Land_sink1_15S35S=reshape(Land_sink1_15S35S,[],1);
Land_sink2(:,6)=Land_sink1_15S35S;

Land_sink1_35S55S=sum(sum(Land_sink1(Threshold_latitude(-35):Threshold_latitude(-55),:)/(10^12),'omitnan'),'omitnan');
Land_sink1_35S55S=reshape(Land_sink1_35S55S,[],1);
Land_sink2(:,7)=Land_sink1_35S55S;