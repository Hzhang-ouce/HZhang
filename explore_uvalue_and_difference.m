clc
clear
load I:\Imperial_college\Master_Project\umodel_experiment\Experiment_output_C.mat
u_value=ncread('NPP\half_saturation_u_360720_new.nc','u');
Difference=Annual_record{1, 2}.Difference_map(:,:,114);
Uniform_missing_pixel=u_value+Difference;
Uniform_missing_pixel(~isnan(Uniform_missing_pixel))=0;
Difference=Difference+Uniform_missing_pixel;
u_value=u_value+Uniform_missing_pixel;
Num_iter=sum(sum(~isnan(u_value)));
Record_pool=zeros(Num_iter,4);
Order=1;
for rows=1:size(Difference,1)
    for columns=1:size(Difference,2)
        if isnan(Difference(rows,columns))
        else
            Record_pool(Order,1)=rows;
            Record_pool(Order,2)=columns;
            Record_pool(Order,3)=Difference(rows,columns);
            Record_pool(Order,4)=u_value(rows,columns);
            Order=Order+1;
            if mod(Order,10)==0
                disp(Order)
            end
        end
    end
end