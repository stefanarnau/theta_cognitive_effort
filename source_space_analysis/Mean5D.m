function DataOut=Mean5D(Data,MeanDims)
% If you have data with several dimensions and you want to calculate means over all or several dimensions, it can be annoying to handle all the brackets
% and squeezes. Therefore, this function shall assist this procedure, in order to more quickly calculate means over several dimensions of data.

% function computes mean over all dimensions contained in MeanDims, starting from 5th dimension
% MeanDims ist vector with integers , e.g. MeanDims=[2 5 4 3 1] means that a mean shall be calculated over the 2nd, 5th, 4th, 3rd and 1st
% dimension of the data.
% Data can contain NANs.


% always starts with biggest MeanDim:
Dims=sort(MeanDims,'descend');

for dim=Dims
    if size(Data,dim)>1
        Data=nanmean(Data,dim);
    end
end
DataOut=squeeze(Data);
end