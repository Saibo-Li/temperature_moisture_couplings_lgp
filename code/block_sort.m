% 生成一个随机的4000*4000*200矩阵
matrix = randi(10,10,10,3);

% 定义每个块的大小
blockSize = 5;

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
            sortedBlock = sort(block,3);
            
            % 将排序后的块放回原始位置
            sortedMatrix(rowRange, colRange, pageRange) = reshape(sortedBlock, numel(rowRange), numel(colRange), numel(pageRange));
        end
    end
end

% 合并排序后的块
sortedMatrix = reshape(sortedMatrix, numRows, numCols, numPages);

% 检查排序结果是否正确
isSorted = all(all(diff(sortedMatrix, 1, 3) >= 0));
disp(['排序结果是否正确：', num2str(isSorted)]);