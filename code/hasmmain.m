clear all;
global wuchaid CGerror CGmaxtimes guimo caiyangid lamdazhi lamdaid pinghuaid hasmtime hasmerror yuchuliid banjing songchi raodong minzhi maxzhi luid;
CGerror=1e-9;  %%% PCG迭代误差终止参数
CGmaxtimes=5000;  %% PCG迭代最大迭代次数
guimo=2000000;  %%%% 根据计算规模，判断代数方程组是直接求解还是PCG迭代
lamdaid=0;%%% 各采样点权重，0表示所有采样点一样，1表示本采样点周边各点的方差，2表示本采样点与周围点均值的差，
hasmerror=1e-12;  %%% HASM 迭代误差终止参数
yuchuliid=1; %%% 预处理的方法 取1或者3
wuchaid=1;
raodong=0;%% 确定在有采样点的网格，其值在采样值上的扰动[0 1)

canshus=[0.2 0.0 1 1 3 1 20 1];     
lamdazhi=canshus(1,1);%%% lamda 的值，大于0，一般不取太大(<10)**********
songchi=canshus(1,2);%% 松弛系数，决定某网格点上下界时，如何根据其邻点极值进行松弛决定,取值范围为0到1之间[0.1 0.2 0.3 0.4 0.5 1]*******
jizhiid=canshus(1,3);%%%极值ID确定总体极值应该如何放松取值为[1 2 3]
pinghuaid=canshus(1,4);%%%取1或者0，决定是否采用平滑措施***********
caiyangid=canshus(1,5);%%% 采样处理方式，取[1 3 5]之一，表示不同的截断阶数
hasmtime=canshus(1,6);%%% HASM迭代最大迭代数***********
banjing=canshus(1,7);%%%计算上下届时搜索的邻点数 一般大于3，20以内***********
luid=canshus(1,8);%%是否用上下界控制，取值为0表示不用，取值为1表示用***********
if pinghuaid==1
    lamdaid=2;
end
mulu='D:\01文件\07 研究生\25HASM培训\guokeda\wugongshandata\';
caiyangshuju=strcat(mulu,'hasm_res.xls');  %%%% 采样数据
testshuju=[];   %%%% 检验数据        
chuzhishuju=strcat(mulu,'idw.txt');  %%%% 初值
shuchumulu=strcat(mulu,'gwr_1000_res_hasm_new8.txt');  %%%% 输出目录

[x0,y0,H,mx,my,h,caiyangs,test,chazhi]=shujuxls(caiyangshuju,testshuju,chuzhishuju);%%%%如果采样数据为excel格式，用该程序
%[x0,y0,H,mx,my,h,caiyangs,test,chazhi]=shuju(caiyangshuju,testshuju,chuzhishuju);%%%%如果采样数据为txt格式，用该程序

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
jieguo=hasm(chazhi,caiyang,test,h);
%[jieguo,error,time]=hasmtest(chazhi,caiyang,test,h);
toc;
%[X0,Y0,H,mx,my,exwgs]=shuju2(strcat(mulu,'exwgs.txt'));
%jieguo(find(exwgs==-9999))=-9999;
jieguoshuchu(jieguo,shuchumulu,mx,my,x0,y0,H);




