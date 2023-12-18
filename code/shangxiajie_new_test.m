
start_time = tic;

site_row =  struct2array(load('y.mat'));
site_col =  struct2array(load('x.mat'));
site_val =  struct2array(load('z.mat'));
site_n = length(site_row); % 站点数

nrow = 149-2; % 行数
ncol = 126-2;% 列数
h = 0.040000000000000;
arr = ones(ncol,nrow);
arr_row = cumsum(arr,1); %按照维度1，累积求和
arr_col = cumsum(arr,2);

arr_site_distance = zeros(ncol,nrow, site_n); %储存站点到每个栅格的距离
arr_site_value = zeros(ncol,nrow, site_n); %储存每个栅格距离最近站点的value值
arr_row_col_site_value_t = zeros(ncol,nrow); %临时变量
banjing=20; %搜索半径
songchi=0;
pinghuaid=1;
maxzhi = 3.320000000000000;
minzhi = -3.600000000000000;
l=zeros(ncol,nrow);
u=zeros(ncol,nrow);

%计算每个站点到栅格的距离，储存在三维的数组中，第三维是站点维
for i = 1:site_n
%     fprintf('%d\n', i);
    site_row_i = site_row(i);
    site_col_i = site_col(i);
    
    arr_row_site = (arr_row*h-site_row_i).^2;
    arr_col_site = (arr_col*h-site_col_i).^2;
    arr_site_distance(:,:,i) = arr_row_site+arr_col_site;
end
% arr_row_col_site_min = min(arr_row_col_site,[],3); %按照维度3，求最小值
% arr_row_col_site==arr_row_col_site_min 找到最小值的位置
% mean(arr_row_col_site,3)~=arr_row_col_site_min
% 如果只找到最小值的位置，赋值为NaN，但是如果某个点的值都相等，则找到的最小值也是这个值，最终会导致这些位置输出结果是NaN
% arr_row_col_site(mean(arr_row_col_site,3)~=min(arr_row_col_site,[],3) & arr_row_col_site==min(arr_row_col_site,[],3)) = NaN;
% ttt = min(arr_row_col_site,[],3);

arr_site_distance_min = zeros(ncol,nrow); %储存最小值
arr_site_distance_minX = zeros(ncol,nrow,banjing);%储存u个最小值
for v =1:banjing
    if v==1
        arr_site_distance_min = min(arr_site_distance,[],3);
    else 
        arr_site_distance(nanmean(arr_site_distance,3)~=arr_site_distance_min & arr_site_distance==arr_site_distance_min) = NaN;
        arr_site_distance_min = min(arr_site_distance,[],3);
    end
    arr_site_distance_minX(:,:,v) = arr_site_distance_min;
    
    for j =1:site_n
        arr_row_col_site_value_t((arr_site_distance(:,:,j)-arr_site_distance_minX(:,:,v))==0)=site_val(j);%将每个栅格赋值为离他最近的值
    end
    arr_site_value(:,:,v) = arr_row_col_site_value_t;
end

%方差处理
kongzhi=songchi;
for k=4:banjing
    va(:,:,k-3)=var(arr_site_value(:,:,1:k),[],3);
end
[minva,idk]=min(va,[],3); clear minva;
sousuobanjing=idk+3;
newz=arr_site_value(:,:,1:sousuobanjing);

if pinghuaid>=1
    newz=saixuan(newz);%%%删除一个对方差影响最大的值
    if pinghuaid==2
        newz=saixuan(newz);
    end
end

minz=min(newz,[],3);
maxz=max(newz,[],3);

minz(minz<0) = minz(minz<0)+kongzhi;
minz(minz>=0) = minz(minz>=0)-kongzhi;
maxz(maxz<0) = maxz(maxz<0)-kongzhi;
maxz(maxz>=0) = maxz(maxz>=0)+kongzhi;
u=min(maxzhi,maxz);
l=max(minzhi,minz);

% 记录结束时间
end_time = toc(start_time);
% 打印排序时间
disp(['Sorting time for iteration ' ': ' num2str(end_time) ' seconds']);

function newz=saixuan(z); %%%把对一系列样点中，方差影响最大的点删除
fangcha=var(z,[],3);%%%所有点的方差

for i=1:size(z,3)
    newz=z;
    newz(:,:,i)=[];
    newf(:,:,i)=abs(var(z,[],3)-fangcha); %%%%该样点删了以后余下点的方差，与全部样点的方差的对比
end
[nf,id]=max(newf,[],3); clear nf;
if size(id,3)==size(newf,3)
    newz=z;
else
    z(:,:,id)=[];
    newz=z;
end
end