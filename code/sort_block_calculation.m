%  Generate a random 4000*4000*200 matrix
matrix = randi(10,10,10,3);

% Define the size of each block
blockSize = 5;

% Get the dimensions of the matrix
[numRows, numCols, numPages] = size(matrix);

% Split the matrix into blocks and sort each block
sortedMatrix = zeros(numRows, numCols, numPages);
for i = 1:blockSize:numRows
    for j = 1:blockSize:numCols
        for k = 1:blockSize:numPages
            % Get the range of the current block
            rowRange = i:min(i+blockSize-1, numRows);
            colRange = j:min(j+blockSize-1, numCols);
            pageRange = k:min(k+blockSize-1, numPages);
            
            % Extracts the current block and sorts it
            block = matrix(rowRange, colRange, pageRange);
            sortedBlock = sort(block,3);
            
            % Putting the sorted block back in its original position
            sortedMatrix(rowRange, colRange, pageRange) = reshape(sortedBlock, numel(rowRange), numel(colRange), numel(pageRange));
        end
    end
end

% Merge sorted blocks
sortedMatrix = reshape(sortedMatrix, numRows, numCols, numPages);

% Check if the sorting result is correct
isSorted = all(all(diff(sortedMatrix, 1, 3) >= 0));
disp(['Check if the sorting result is correct£º', num2str(isSorted)]);