% Author: Paul Antony 20140408

function [MaskFiltered] = f_RemoveBigObjects(Mask, AreaThreshold)
%This function removes big objects from a mask
%   Mask: binary image
%   AreaThreshold: objects with more than 'AreaThreshold pixels will be removed'
%   MaskFiltered: returned filtered mask

cc = bwconncomp(Mask);
ccprops = regionprops(cc, 'area', 'PixelIdxList');
 %Handle the case where the length of ccprops is 1
    if length(ccprops) == 1
        ccTable = struct2cell(ccprops);
        ccTable = ccTable';
        ccTable = cell2table(ccTable);
        ccTable.Properties.VariableNames = {'Area', 'PixelIdxList'};
    else
        ccTable = struct2table(ccprops);
    end
ccFilteredRows = ccTable.Area < AreaThreshold;
ccFilteredVars = {'Area','PixelIdxList'};
ccFiltered = ccTable(ccFilteredRows,ccFilteredVars);
Objects = zeros(size(Mask));

TableDims = size(ccFiltered);
MaskIndex = [];
for i = 1:TableDims(1)
    MaskIndex = vertcat(MaskIndex, ccFiltered.PixelIdxList{i});
end
Objects(MaskIndex) = 1;
MaskFiltered = Objects;

end

