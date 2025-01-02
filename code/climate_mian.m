clear all;
global wuchaid CGerror CGmaxtimes guimo caiyangid lamdazhi lamdaid pinghuaid hasmtime hasmerror yuchuliid banjing songchi raodong minzhi maxzhi luid;
CGerror=1e-9;  %%% PCG iterative error termination parameter
CGmaxtimes=5000;  %% Maximum number of PCG iterations
guimo=2000000;  %%%% According to the calculation scale, we can judge whether the system of algebraic equations is solved directly or by PCG iteration.
lamdaid=0;%%% The weight of each sampling point, 0 indicates that all sampling points are the same, 1 represents the variance of each point around the sampling point, and 2 represents the difference of the mean value between the sampling point and the surrounding point.
hasmerror=1e-12;  %%%Iterative error termination parameter
yuchuliid=1; %%% The method of pretreatment is 1 or 3.
wuchaid=1;
raodong=0;%% Determine the disturbance of the value of the grid with sampling points on the sampling value [0 1)

canshus=[0.2 0.0 1 1 3 1 20 1];     
lamdazhi=canshus(1,1);%%% The value of lamda, which is greater than 0, is not too large (<10)**********
songchi=canshus(1,2);%% When determining the upper and lower bound of a grid point, the relaxation coefficient is determined according to the extreme value of its adjacent point, which ranges from 0 to 1.[0.1 0.2 0.3 0.4 0.5 1]*******
jizhiid=canshus(1,3);%%%The extreme value ID determines how the global extreme value should be relaxed. The value is[1 2 3]
pinghuaid=canshus(1,4);%%%Take 1 or 0 to decide whether to use smoothing measures.***********
caiyangid=canshus(1,5);%%% Sampling processing mode, take one of [1 3 5] to represent different truncation orders
hasmtime=canshus(1,6);%%% 
banjing=canshus(1,7);%%%The number of neighbor points in the search is generally less than 20 at that time.***********
luid=canshus(1,8);%% Whether to use the upper and lower bound control. A value of 0 means no use, and a value of 1 means using***********
if pinghuaid==1
    lamdaid=2;
end
mulu='data\';
caiyangshuju=strcat(mulu,'hasm_res.xls');  %%%% Sampled data
testshuju=[];   %%%% Test data     
chuzhishuju=strcat(mulu,'initial.txt'); 
shuchumulu=strcat(mulu,'output.txt');  %

[x0,y0,H,mx,my,h,caiyangs,test,chazhi]=shujuxls(caiyangshuju,testshuju,chuzhishuju);%%
%[x0,y0,H,mx,my,h,caiyangs,test,chazhi]=shuju(caiyangshuju,testshuju,chuzhishuju);%%%%

caiyang=[caiyangs(:,1) caiyangs(:,2) caiyangs(:,3)];
%test=[tests(:,1) tests(:,2) tests(:,4)];
minzzz=min(caiyang(:,3));
maxzzz=max(caiyang(:,3));
if jizhiid==1
    minzhi=minzzz;
    maxzhi=maxzzz;
elseif jizhiid==2
    if minzzz<0
        minzhi=1.1*minzzz;
    else
        minzhi=0.9*minzzz;
    end
    if maxzzz<0
        maxzhi=0.9*maxzzz;
    else
        maxzhi=1.1*maxzzz;
    end
else
    if minzzz<0
        minzhi=1.2*minzzz;
    else
        minzhi=0.8*minzzz;
    end
    if maxzzz<0
        maxzhi=0.8*maxzzz;
    else
        maxzhi=1.2*maxzzz;
    end
end
clear maxzzz minzzz;
tic;
jieguo=climate_processing_fun(chazhi,caiyang,test,h);
%[jieguo,error,time]=hasmtest(chazhi,caiyang,test,h);
toc;
%[X0,Y0,H,mx,my,exwgs]=shuju2(strcat(mulu,'exwgs.txt'));
%jieguo(find(exwgs==-9999))=-9999;
jieguoshuchu(jieguo,shuchumulu,mx,my,x0,y0,H);




