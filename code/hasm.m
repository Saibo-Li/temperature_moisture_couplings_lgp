function [jieguo,error]=hasm(chazhi,caiyang,test,h);
if isempty(test)
    jieguo=hasmguding(chazhi,caiyang,h); %%%固定迭代次数HASM计算
    error='NaN';
else
    [jieguo,error,time]=hasmtest(chazhi,caiyang,test,h);%%%有验证数据的HASM计算
end