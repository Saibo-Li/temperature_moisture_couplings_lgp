function [jieguo,error,time]=hasmtest(chazhi,caiyang,test,h);
global caiyangid guimo hasmtime hasmerror yuchuliid luid quanzhong;
%%%%  固定迭代次数，无误差检验的HASM迭代程序
[mx,my]=size(chazhi);
caiyang(caiyang(:,1)<(0.5*h),:)=[];
caiyang(caiyang(:,1)>((mx-2)*h+0.5*h),:)=[];
caiyang(caiyang(:,2)<(0.5*h),:)=[];
caiyang(caiyang(:,2)>((my-2)*h+0.5*h),:)=[]; %%%把超出范围的采样点去掉
test(test(:,1)<(0.5*h),:)=[];
test(test(:,1)>((mx-2)*h+0.5*h),:)=[];
test(test(:,2)<(0.5*h),:)=[];
test(test(:,2)>((my-2)*h+0.5*h),:)=[]; %%%把超出范围的采样点去掉

if (isempty(caiyang(:,1))==1)
    disp('区域内没有采样点，请检查程序和数据是否有错！');
    jieguo=chazhi;
    error=wuchajianyan(jieguo,test,h);
    time=0;
else
    chazhi0=chazhi;
    time=0;
    error(2)=wuchajianyan(chazhi,test,h);
    error(1)=error(2)*2;
    lamda=lisangdu(caiyang(:,1),caiyang(:,2),caiyang(:,3));
    id=find(caiyang(:,4)==0);
    lamda(id)=quanzhong*lamda(id);
    if ismember(caiyangid,[1 3 5])==0
        disp('caiyangid赋值错误(只能是1、3或者5)，程序自动改成3');
        caiyangid=3;
    end
    [GG,GH]=taylorcaiyang(caiyang(:,1),caiyang(:,2),caiyang(:,3),chazhi,lamda,h);
    if ((mx-2)*(my-2))<guimo
        xishus=taylorxishu(mx,my,GG);
        if luid==1
            [l,u]=shangxiajie2(caiyang,mx,my,h);
        end
        while diedaizhi(error(time+1),error(time+2),time,hasmtime,hasmerror)==1
            time=time+1;
            chazhi0=chazhi;
            b=quyumoni(chazhi,h)+GH;  
            x=xishus\b;
            clear b;
            chazhi(2:(mx-1),2:(my-1))=(reshape(x,(my-2),(mx-2)))';
            clear x;
            if luid==1
                for i=1:(mx-2)
                    for j=1:(my-2)
                        if chazhi(i+1,j+1)<l(i,j)
                            chazhi(i+1,j+1)=l(i,j);
                        elseif chazhi(i+1,j+1)>u(i,j)
                            chazhi(i+1,j+1)=u(i,j);
                        end
                    end
                end
            end
            error(time+2)=wuchajianyan(chazhi,test,h);
        end
    else 
        M=duijiaoxian(mx,my); 
        yuxishu=yuchulixishu(mx,my,GG,M,yuchuliid); 
        if luid==1
            [l,u]=shangxiajie2(caiyang,mx,my,h);
        end
        while diedaizhi(error(time+1),error(time+2),time,hasmtime,hasmerror)==1
                chazhi0=chazhi;
                time=time+1;
                b=quyumoni(chazhi,h)+GH;  
                xx=reshape((chazhi(2:(mx-1),2:(my-1)))',(mx-2)*(my-2),1);
                xx=tayloryuchuliCG(mx,my,b,xx,yuxishu,GG,M);                  
                clear b;
                chazhi(2:(mx-1),2:(my-1))=(reshape(xx,(my-2),(mx-2)))';
                clear xx;   
                if luid==1
                    for i=1:(mx-2)
                        for j=1:(my-2)
                            if chazhi(i+1,j+1)<l(i,j)
                                chazhi(i+1,j+1)=l(i,j);
                            elseif chazhi(i+1,j+1)>u(i,j)
                                chazhi(i+1,j+1)=u(i,j);
                            end
                        end
                    end
                end
                error(time+2)=wuchajianyan(chazhi,test,h);
 %           xinxi=strcat('正在进行HASM迭代的第',num2str(time),'次迭代');
  %          disp(xinxi);
         end
    end 
    jieguo=chazhi;
    error(1)=[];
    if error(length(error))<error(length(error)-1)
        jieguo=chazhi;
    else
        if time>1
            jieguo=chazhi0;
            error(length(error))=[];
            time=time-1;
        end
    end
    clear chazhi chazhi0;
end

function id=diedaizhi(error1,error2,time,hasmtime,hasmerror);
id=1;
if error1<error2
    id=0;
else
    if error2<hasmerror
        id=0;
    else
        if time>hasmtime
            id=0;
        else
            id=1;
        end
    end
end

function yuxishu=yuchulixishu(mx,my,GG,M,yuchuliid); %%%%生成三对角阵，作为预处理子，本矩阵为A'A＋B'B+G'G 的中间三对角 ******
if yuchuliid==3
    n=(mx-2)*(my-2);
    e=ones(n,1);
    pxia=e;
    pshang=e;
    clear e;
    mx=mx-2;
    my=my-2;
    pxia(my:my:(n-my))=0;
    pshang((my+1):my:(n-my+1))=0;
    yuxishu=spdiags([-4*pxia M -4*pshang],-1:1,n,n);
    GG=GG+yuxishu;
    clear yuxishu;
    ss=[0 (diag(GG,1))']';
    sx=[(diag(GG,-1))' 0]';
    yuxishu=spdiags([sx diag(GG) ss],-1:1,n,n); %%%% A'A+B'B+GG的三对角，作为预处理子
    clear ss sx;
elseif yuchuliid==1
    n=(mx-2)*(my-2);
    yuxishu=zeros(n,1);
    for i=1:n
        yuxishu(i)=M(i)+GG(i,i);
    end
end

function M=duijiaoxian(mx,my); 
m=(mx-2)*(my-2);
M=12*ones(m,1);
M(1:(my-2):(m+3-my))=11;
M((my-2):(my-2):m)=11;
M(1:(my-2))=M(1:(my-2))-1;
M((m-my+3):m)=M((m-my+3):m)-1;

function jieguo=taylorxishu(mx,my,GG); %%%%生成全系数矩阵，可以直接求逆，本矩阵为A'A＋B'B +GG ******
n=(mx-2)*(my-2);
e=ones(n,1);
pxia=e;
pshang=e;
mx=mx-2;
my=my-2;
pxia(my:my:(n-my))=0;
pshang((my+1):my:(n-my+1))=0;
xishu1=spdiags([pxia -2*e pshang],-1:1,n,n);
xishu2=spdiags([e -2*e e],[-my 0 my],n,n);
clear pxia pshang e;
jieguo=xishu1'*xishu1+xishu2'*xishu2+GG;
clear xishu1 xishu2;


    
function x=tayloryuchuliCG(mx,my,b,x0,yuxishu,GG,M);   %%%%%对角线预处理CG迭代方法 yuxishu 为三对角预处理阵，而不是原方程系数矩阵
global CGerror CGmaxtimes;
m=(mx-2)*(my-2);
mob=sqrt(sum(b.^2)/m);%%% b的模
r=b-AP(GG,x0,mx,my,M);
clear b;
%newr=r;
ruo=sum(r.^2)/m;
x=x0;
clear x0;
k=0;
if  min(size(yuxishu))==1
    z=r./yuxishu;
else
    z=yuxishu\r;
end
while ((sqrt(ruo)>(CGerror*mob))&(k<CGmaxtimes))
    %newz=newr./M; 
    k=k+1;
    if k==1
        p=z;
    else
        beita=newrz/rz;
        p=z+beita*p;
    end
    rz=r'*z; 
    w=AP(GG,p,mx,my,M);
    alfa=rz/(p'*w);  %    w=AP(M,p,mx,my,lamda);
    x=x+alfa*p;
    r=r-alfa*w;
    if  min(size(yuxishu))==1
        z=r./yuxishu;
    else
        z=yuxishu\r;
    end
    newrz=r'*z;
    ruo=sum(r.^2)/m;
end %%%%%迭代过程中需要存储 r z x p w 5个矩阵 做6*m次乘法
clear r p z w yuxishu;

function pp=AP(GG,p,mx,my,M);
pp=zeros((mx-2)*(my-2),1);
pp=xiangliangweiyi2(p,mx,my,2)+xiangliangweiyi2(p,mx,my,-2)-4*xiangliangweiyi2(p,mx,my,1)-4*xiangliangweiyi2(p,mx,my,-1);
pp=pp-4*xiangliangweiyi1(p,my-2)+xiangliangweiyi1(p,2*(my-2))-4*xiangliangweiyi1(p,-(my-2))+xiangliangweiyi1(p,-2*(my-2));
pp=pp+GG*p+M.*p;

function temp=xiangliangweiyi1(p,a);
m=max(size(p));
temp=zeros(m,1);
if a>=0
    temp(1:(m-a))=p((a+1):m); %%%%对角线上方
elseif a<0
    temp((-a+1):m)=p(1:(m+a)); %%%%对角线下方
end

function temp=xiangliangweiyi2(p,mx,my,a);
m=(mx-2)*(my-2);
temp=zeros(m,1);
if a==1
    p(1:(my-2):(m-my+3))=0;
elseif a==2
    p((my-1):(my-2):(m-my+3))=0;
    p(my:(my-2):(m-my+4))=0;
elseif a==(-1)
    p((my-2):(my-2):(m-my+2))=0;
elseif a==(-2)
    p((my-2):(my-2):m)=0;
    p((my-3):(my-2):(m-1))=0;
end
temp=xiangliangweiyi1(p,a);