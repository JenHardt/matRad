function dataOut2 = matRad_readBinData(binFile,cubeDim)
fID = fopen(binFile);
data = fread(fID,inf,'double');
fclose(fID);

if rem(numel(data),prod(cubeDim))==0
    numOfFields = numel(data)/prod(cubeDim);
    for i = 1:numOfFields
        dataOut{i} = reshape(data(i:numOfFields:end),cubeDim(2),cubeDim(1),cubeDim(3));
        dataOut{i} = permute(dataOut{i},[2 1 3]);
    end
else
    error('bin data contains an odd number of entries.')
end
dataOut2 = dataOut{1};

end
