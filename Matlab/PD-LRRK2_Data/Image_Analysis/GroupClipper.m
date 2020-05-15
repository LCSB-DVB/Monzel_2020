function [xClipped] = GroupClipper(x)
%Remove empty rows and columns from cell array

    CellsWithData = cellfun(@(c) ~isempty(c), x);
    ColumnsToKeep = find(sum(CellsWithData, 1) ~= 0);
    RowsToKeep = find(sum(CellsWithData, 2) ~= 0);
    xClipped = x(RowsToKeep, ColumnsToKeep);

end