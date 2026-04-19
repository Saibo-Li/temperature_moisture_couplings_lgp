function [l,u]=range_raster_value(caiyang,mx,my,h);%%%%Get the upper and lower bounds of each raster point
global banjing songchi raodong minzhi maxzhi pinghuaid; 
m=(mx-2)*(my-2);
x=caiyang(:,1);
y=caiyang(:,2);
z=caiyang(:,3);
cyminzhi=min(z);%%Extreme values in sampling points
cymaxzhi=max(z);
if minzhi>cyminzhi
    disp('minzhi The setting is too large and the program automatically adjusts to the minimum value in the sampling set!');
    minzhi=cyminzhi;
end
if maxzhi<cymaxzhi
    disp('maxzhi The setting is too small and the program automatically adjusts to the maximum value in the sampling set!');
    mazhi=cymaxzhi;
end
if (minzhi<cyminzhi*(1-songchi))
    disp('songshi Maybe the settings are too big?');
elseif (maxzhi>cymaxzhi*(1+songchi))
    disp('songshi Maybe the settings are too big?');
end
if raodong>songchi
    raodong=songchi/10;
    disp('raodong It may be set too large, the program automatically adjusts raodong to one tenth of songchi!');
end

l=zeros(mx-2,my-2);
u=zeros(mx-2,my-2);
%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
site_row =  y;
site_col =  x;
site_val =  z;
site_n = length(site_row); % Number of sites

nrow = my-2; % 
ncol = mx-2;% 

arr = ones(ncol,nrow);
arr_row = cumsum(arr,2); %Cumulative summation according to dimension 2
arr_col = cumsum(arr,1);

arr_site_distance = zeros(ncol,nrow,site_n); %Store the distance from the site to each raster
arr_site_value = zeros(ncol,nrow, banjing); %Stores the VALUE value of each raster's distance from the nearest station

%The distance from each site to the raster is calculated and stored in a three-dimensional array, with the third dimension being the site dimension
for i = 1:site_n
    site_row_i = site_row(i);
    site_col_i = site_col(i);
    
    arr_row_site = (arr_row*h-site_row_i).^2;
    arr_col_site = (arr_col*h-site_col_i).^2;
    arr_site_distance(:,:,i) = arr_row_site+arr_col_site;
end

if ncol*nrow < 25000000
    start_time = tic;
    arr_row_col_site_value_t = zeros(ncol,nrow); 
    arr_site_distance_minX = zeros(ncol,nrow,banjing);
    arr_site_distance_mean = nanmean(arr_site_distance,3);
    
    for v =1:banjing
        if v==1
            arr_site_distance_min = min(arr_site_distance,[],3);
        else 
            arr_site_distance(arr_site_distance_mean~=arr_site_distance_min & arr_site_distance==arr_site_distance_min) = NaN;
            arr_site_distance_min = min(arr_site_distance,[],3);
        end
        arr_site_distance_minX(:,:,v) = arr_site_distance_min;

        for j =1:site_n
            arr_row_col_site_value_t((arr_site_distance(:,:,j)-arr_site_distance_minX(:,:,v))==0)=site_val(j);
        end
        arr_site_value(:,:,v) = arr_row_col_site_value_t;
    end
    end_time = toc(start_time);
    disp(['Panel time for iteration ' ': ' num2str(end_time) ' seconds']);
else
    start_time = tic;
    matrix = arr_site_distance;
    blockSize = 1000;
    [numRows, numCols, numPages] = size(matrix);
    sortedMatrix = zeros(numRows, numCols, banjing);
    for i = 1:blockSize:numRows
        for j = 1:blockSize:numCols
            rowRange = i:min(i+blockSize-1, numRows);
            colRange = j:min(j+blockSize-1, numCols);
            pageRange = 1:site_n;

            block = matrix(rowRange, colRange, pageRange);
            sortedBlock =zeros(length(rowRange),length(colRange), banjing); 
            arr_site_distance_minX = zeros(length(rowRange),length(colRange),banjing);
            arr_site_distance_mean = nanmean(block,3);
            arr_row_col_site_value_t = zeros(length(rowRange),length(colRange)); 

            for v =1:banjing
                
                if v==1
                    arr_site_distance_min = min(block,[],3);
                else 
                    block(arr_site_distance_mean~=arr_site_distance_min & block==arr_site_distance_min) = NaN;
                    arr_site_distance_min = min(block,[],3);
                end
                arr_site_distance_minX(:,:,v) = arr_site_distance_min;
                
                for u =1:site_n
                    arr_row_col_site_value_t((block(:,:,u)-arr_site_distance_minX(:,:,v))==0)=site_val(u);
                end
                sortedBlock(:,:,v) = arr_row_col_site_value_t;
                
            end
            sortedMatrix(rowRange, colRange, 1:banjing) = reshape(sortedBlock, numel(rowRange), numel(colRange), banjing);
        end
    end
    arr_site_value = sortedMatrix;
    end_time = toc(start_time);
    disp(['Block time for iteration ' ': ' num2str(end_time) ' seconds']);
end

kongzhi=songchi;
for k=4:banjing
    va(:,:,k-3)=var(arr_site_value(:,:,1:k),[],3);
end

[minva,idk]=min(va,[],3); clear minva;
sousuobanjing=idk+3;
newz_var = ones(ncol,nrow,banjing);
newz_var_cumsum = cumsum(newz_var,3);
arr_site_value(newz_var_cumsum>sousuobanjing) = NaN;
newz = arr_site_value;

if pinghuaid>=1
    newz=saixuan(newz,newz_var_cumsum);
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
% Record end time
end_time = toc(start_time);
disp(['Shangxiejie time for iteration ' ': ' num2str(end_time) ' seconds']);
%%%%%%%%%%%%%%%%%%%%%%%%%%
if pinghuaid==1
    l=(smoothts(l,'b',round(my/10))+(smoothts(l','b',round(mx/10)))')/2;
    u=(smoothts(u,'b',round(my/10))+(smoothts(u','b',round(mx/10)))')/2;
end

function newz=saixuan(z,newz_var_cumsum); %%%Remove the point that has the greatest effect on the variance in a series of sample points
fangcha=nanvar(z,[],3);%%%Variance at all points

for i=1:size(z,3)
    newz=z;
    newz(:,:,i)=[];
    newf(:,:,i)=abs(nanvar(newz,[],3)-fangcha); %%%%Comparison of the variance of the remaining points after deletion of the sample point with the variance of all the sample points.
end
[nf,id]=max(newf,[],3);clear nf;
if size(id,3)==size(newf,3)
    newz=z;
else
    z(newz_var_cumsum==id)=NaN;
    newz=z;
end