function [jieguo,error]=climate_processing_fun(chazhi,caiyang,test,h);
if isempty(test)
    jieguo=climate_processing_iteration.m(chazhi,caiyang,h); %%%Fixed number of iterations HASM calculation
    error='NaN';
else
    [jieguo,error,time]=climate_validation(chazhi,caiyang,test,h);%%%HASM calculations with validation data
end