u_new = struct2array(load('u_new.mat'));
l_new = struct2array(load('l_new.mat'));

u_old = struct2array(load('u_old.mat'));
l_old = struct2array(load('l_old.mat'));

c = cumsum(ones(3,3,3),3);
idd = [1,2,3;3,2,2;1,1,1];

c(:,:,idd)=NaN

iddd = [1,1,1;1,1,1;1,1,1];
a = [1,3,4;]

for i=1:100
   c(:,:,i) =randi(10,3000,3000);
end

c=randi(10,2000,2000,100);
[c_value,c_index] = sort(c,3);
% 生成一个随机的4000*4000*200矩阵
matrix = rand(4000, 4000, 200);

% 定义每个块的大小
blockSize = 1000;

% 获取矩阵的维度
[numRows, numCols, numPages] = size(matrix);

% 分割矩阵为块，并对每个块进行排序
sortedMatrix = zeros(numRows, numCols, numPages);
for i = 1:blockSize:numRows
    for j = 1:blockSize:numCols
        for k = 1:blockSize:numPages
            % 获取当前块的范围
            rowRange = i:min(i+blockSize-1, numRows);
            colRange = j:min(j+blockSize-1, numCols);
            pageRange = k:min(k+blockSize-1, numPages);
            
            % 提取当前块，并对其进行排序
            block = matrix(rowRange, colRange, pageRange);
            sortedBlock = sort(block(:));
            
            % 将排序后的块放回原始位置
            sortedMatrix(rowRange, colRange, pageRange) = reshape(sortedBlock, numel(rowRange), numel(colRange), numel(pageRange));
        end
    end
end

% 合并排序后的块
sortedMatrix = reshape(sortedMatrix, numRows, numCols, numPages);

% 检查排序结果是否正确
isSorted = all(all(diff(sortedMatrix, 1, 1) >= 0) & all(diff(sortedMatrix, 1, 2) >= 0) & all(diff(sortedMatrix, 1, 3) >= 0));
disp(['排序结果是否正确：', num2str(isSorted)]);
 