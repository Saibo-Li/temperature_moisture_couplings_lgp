function main
mx = 1119; %列
my = 1317; %行
x = struct2array(load('x.mat'));
y = struct2array(load('y.mat'));
z = struct2array(load('z.mat'));
h = 8.944543828264759e-04;
raodong=0;
songchi=0;
banjing=20;
pinghuaid=1;
maxzhi = 3.320000000000000;
minzhi = -3.600000000000000;

xid=round2(x./h); %%% 距离采样点(x,y)最近的网格点在所有内点网格中的x序号
yid=round2(y./h); %%% 距离采样点(x,y)最近的网格点在所有内点网格中的y序号
caiyangxulie=(xid-1).*(my-2)+yid; %%%在所有内点中的排序
l=zeros(mx-2,my-2);
u=zeros(mx-2,my-2);
start_time = tic;
for i=1:(mx-2)
    for j=1:(my-2)
%         disp([num2str(i), '，',num2str(j)]);
        xuhao=(i-1)*(my-2)+j;
        juli=(x-i*h).^2+(y-j*h).^2;
        [sortjuli,ID]=sort(juli); clear sortjuli;   
        clear sortjuli;
        if ismember(xuhao,caiyangxulie)==1 %%%本栅格内有采样点，则取离网格中心最近的那一个采样值
            kongzhi=raodong;
            newz1=z(ID(1));
            tt=ID(2:(banjing+1));
            newz2=z(tt);
            newz3=saixuan([newz2' newz1]');
            if pinghuaid==2
                newz3=saixuan(newz3);
            end
            if ismember(newz1,newz3)==1
                newz=newz1;
            else  %%若离网格中心最近的点的值为奇异，则取周围点来计算上下界
                for k=4:banjing
                    va(k-3)=var(newz2(1:k));
                end
                [minva,idk]=min(va);clear minva;
                sousuobanjing=idk+3;
                tt=ID(2:sousuobanjing);
                newz=z(tt);
                if pinghuaid>=1
                    newz=saixuan(newz);%%%删除一个对方差影响最大的值
                    if pinghuaid==2
                        newz=saixuan(newz);%%%再删1个
                    end
                end
                clear tt;        
             end
        else%%%本栅格内没有采样点，则取离这个栅格中心最近的sousuobanjing个采样点
            kongzhi=songchi;
            tt=ID(1:banjing);
            newz=z(tt);
            clear tt;
            for k=4:banjing
                va(k-3)=var(newz(1:k));
            end
            [minva,idk]=min(va); clear minva;
            sousuobanjing=idk+3;
            tt=ID(1:sousuobanjing);
            newz=z(tt);
            if pinghuaid>=1
                newz=saixuan(newz);%%%删除一个对方差影响最大的值
                if pinghuaid==2
                    newz=saixuan(newz);
                end
            end
            clear tt;        
        end
        minz=min(newz);
        maxz=max(newz);
        if minz<0
            minzzz=minz*(1+kongzhi);
        else
            minzzz=minz*(1-kongzhi);
        end
        if maxz<0
            maxzzz=maxz*(1-kongzhi);
        else
            maxzzz=maxz*(1+kongzhi);
        end
        u(i,j)=min(maxzhi,maxzzz);
        l(i,j)=max(minzhi,minzzz);

        clear ID tt newz;
    end
end
% 记录结束时间
end_time = toc(start_time);
% 打印排序时间
disp(['Sorting time for iteration ' ': ' num2str(end_time) ' seconds']);


function b=round2(x)
%%%%%%%最接近的整数，如果跟两个整数距离相等，则取小的那个
a=abs(x-round(x));
if a==0.5
    b=round(x)-1;
else
    b=round(x);
end

function newz=saixuan(z); %%%把对一系列样点中，方差影响最大的点删除
fangcha=var(z);%%%所有点的方差
for i=1:length(z)
    newz=z;
    newz(i)=[];
    newf(i)=abs(var(newz)-fangcha); %%%%该样点删了以后余下点的方差，与全部样点的方差的对比
end
[nf,id]=max(newf); clear nf;
if length(id)==length(newf)
    newz=z;
else
    z(id)=[];
    newz=z;
end