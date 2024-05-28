function A = readMatrixFromFile(filename)
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Cannot open file: %s', filename);
    end

    % tline = fgetl(fileID);
    % while ischar(tline)
    %     disp(tline)
    %     tline = fgetl(fileID);
    % end

    sizeA = fscanf(fileID, '%d %d', 2);

    A = fscanf(fileID, '%f', sizeA(1) * sizeA(2));

    % Reshape the data into a matrix
    A = reshape(A, [sizeA(2), sizeA(1)])';

    fclose(fileID);
end