function [l,u]=shangxiajie(caiyang,mx,my,h);%%%%得到各栅格点的上下界
global banjing songchi raodong minzhi maxzhi pinghuaid; %%banjing这里为最大搜索数
m=(mx-2)*(my-2);
x=caiyang(:,1);
y=caiyang(:,2);
z=caiyang(:,3);
cyminzhi=min(z);%%采样点中的极值
cymaxzhi=max(z);
if minzhi>cyminzhi
    disp('minzhi设置过大，程序自动调整为采样集中最小值！');
    minzhi=cyminzhi;
end
if maxzhi<cymaxzhi
    disp('maxzhi设置过小，程序自动调整为采样集中最大值！');
    mazhi=cymaxzhi;
end
if (minzhi<cyminzhi*(1-songchi))
    disp('songshi或许设置太大？');
elseif (maxzhi>cymaxzhi*(1+songchi))
    disp('songshi或许设置太大？');
end
if raodong>songchi
    raodong=songchi/10;
    disp('raodong可能设置过大，程序自动调整raodong为songchi的十分之一!');
end

l=zeros(mx-2,my-2);
u=zeros(mx-2,my-2);
%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
site_row =  y;
site_col =  x;
site_val =  z;
site_n = length(site_row); % 站点数

nrow = my-2; % 行数
ncol = mx-2;% 列数

arr = ones(ncol,nrow);
arr_row = cumsum(arr,2); %按照维度2，累积求和
arr_col = cumsum(arr,1);

arr_site_distance = zeros(ncol,nrow,site_n); %储存站点到每个栅格的距离
arr_site_value = zeros(ncol,nrow, banjing); %储存每个栅格距离最近站点的value值

%计算每个站点到栅格的距离，储存在三维的数组中，第三维是站点维
for i = 1:site_n
    site_row_i = site_row(i);
    site_col_i = site_col(i);
    
    arr_row_site = (arr_row*h-site_row_i).^2;
    arr_col_site = (arr_col*h-site_col_i).^2;
    arr_site_distance(:,:,i) = arr_row_site+arr_col_site;
end

if ncol*nrow < 25000000
    %%%%%%%%行,或列<5000%%%%%%%%
    start_time = tic;
    arr_row_col_site_value_t = zeros(ncol,nrow); %临时变量
    arr_site_distance_minX = zeros(ncol,nrow,banjing);%储存banjing个最小值
    arr_site_distance_mean = nanmean(arr_site_distance,3);%储存arr_site_distance在第三个维度的平均值
    
    for v =1:banjing
        if v==1
            arr_site_distance_min = min(arr_site_distance,[],3);
        else 
            arr_site_distance(arr_site_distance_mean~=arr_site_distance_min & arr_site_distance==arr_site_distance_min) = NaN;
            arr_site_distance_min = min(arr_site_distance,[],3);
        end
        arr_site_distance_minX(:,:,v) = arr_site_distance_min;

        for j =1:site_n
            arr_row_col_site_value_t((arr_site_distance(:,:,j)-arr_site_distance_minX(:,:,v))==0)=site_val(j);%将每个栅格赋值为离他最近的值
        end
        arr_site_value(:,:,v) = arr_row_col_site_value_t;
    end
    % 记录结束时间
    end_time = toc(start_time);
    disp(['Panel time for iteration ' ': ' num2str(end_time) ' seconds']);
else
    %%%%%%%%行,或列>5000%%%%%%%%
    % 分割矩阵为块，并对每个块单独处理
    start_time = tic;
    matrix = arr_site_distance;
    % 定义每个块的大小
    blockSize = 1000;
    % 获取矩阵的维度
    [numRows, numCols, numPages] = size(matrix);
    sortedMatrix = zeros(numRows, numCols, banjing);
    for i = 1:blockSize:numRows
        for j = 1:blockSize:numCols
            % 获取当前块的范围
            rowRange = i:min(i+blockSize-1, numRows);
            colRange = j:min(j+blockSize-1, numCols);
            pageRange = 1:site_n;

            % 提取当前块，并对其进行排序
            block = matrix(rowRange, colRange, pageRange);
            sortedBlock =zeros(length(rowRange),length(colRange), banjing); %储存每个栅格距离最近站点的value值
            arr_site_distance_minX = zeros(length(rowRange),length(colRange),banjing);%储存banjing个最小值
            arr_site_distance_mean = nanmean(block,3);%储存arr_site_distance在第三个维度的平均值
            arr_row_col_site_value_t = zeros(length(rowRange),length(colRange)); %临时变量

            for v =1:banjing
                
                if v==1
                    arr_site_distance_min = min(block,[],3);
                else 
                    block(arr_site_distance_mean~=arr_site_distance_min & block==arr_site_distance_min) = NaN;
                    arr_site_distance_min = min(block,[],3);
                end
                arr_site_distance_minX(:,:,v) = arr_site_distance_min;
                
                for u =1:site_n
                    arr_row_col_site_value_t((block(:,:,u)-arr_site_distance_minX(:,:,v))==0)=site_val(u);%将每个栅格赋值为离他最近的值
                end
                sortedBlock(:,:,v) = arr_row_col_site_value_t;
                
            end
            % 将排序后的块放回原始位置
            sortedMatrix(rowRange, colRange, 1:banjing) = reshape(sortedBlock, numel(rowRange), numel(colRange), banjing);
        end
    end
    arr_site_value = sortedMatrix;
    % 记录结束时间
    end_time = toc(start_time);
    disp(['Block time for iteration ' ': ' num2str(end_time) ' seconds']);
end

%方差处理
kongzhi=songchi;
for k=4:banjing
    va(:,:,k-3)=var(arr_site_value(:,:,1:k),[],3);
end

[minva,idk]=min(va,[],3); clear minva;
sousuobanjing=idk+3;
%超出sousuobanjing的赋值为NaN
newz_var = ones(ncol,nrow,banjing);
newz_var_cumsum = cumsum(newz_var,3);
arr_site_value(newz_var_cumsum>sousuobanjing) = NaN;
newz = arr_site_value;

if pinghuaid>=1
    newz=saixuan(newz,newz_var_cumsum);%%%删除一个对方差影响最大的值
    if pinghuaid==2
        newz=saixuan(newz,newz_var_cumsum);
    end
end

minz=min(newz,[],3);
maxz=max(newz,[],3);

minz(minz<0) = minz(minz<0)*(1+kongzhi);
minz(minz>=0) = minz(minz>=0)*(1-kongzhi);
maxz(maxz<0) = maxz(maxz<0)*(1-kongzhi);
maxz(maxz>=0) = maxz(maxz>=0)*(1+kongzhi);
u=min(maxzhi,maxz);
l=max(minzhi,minz);
% 记录结束时间
end_time = toc(start_time);
disp(['Shangxiejie time for iteration ' ': ' num2str(end_time) ' seconds']);
%%%%%%%%%%%%%%%%%%%%%%%%%%
if pinghuaid==1
    l=(smoothts(l,'b',round(my/10))+(smoothts(l','b',round(mx/10)))')/2;
    u=(smoothts(u,'b',round(my/10))+(smoothts(u','b',round(mx/10)))')/2;
end

function newz=saixuan(z,newz_var_cumsum); %%%把对一系列样点中，方差影响最大的点删除
fangcha=nanvar(z,[],3);%%%所有点的方差

for i=1:size(z,3)
    newz=z;
    newz(:,:,i)=[];
    newf(:,:,i)=abs(nanvar(newz,[],3)-fangcha); %%%%该样点删了以后余下点的方差，与全部样点的方差的对比
end
[nf,id]=max(newf,[],3);clear nf;
if size(id,3)==size(newf,3)
    newz=z;
else
    z(newz_var_cumsum==id)=NaN;
    newz=z;
end