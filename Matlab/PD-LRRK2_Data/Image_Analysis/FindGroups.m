function [GroupsTable, GroupsIm5DCellArray] = FindGroups(InfoTable)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    XPositions = unique(InfoTable.PositionX);
    YPositions = unique(InfoTable.PositionY);

    XDeltas = round(Deltas(XPositions), 10); % unique(XDeltas)
    XMinDelta = min(XDeltas);

    YDeltas = round(Deltas(YPositions), 10);
    YMinDelta = min(YDeltas);

    % Convert Positions to integers
    InfoTable.PosXInteger = int16(InfoTable.PositionX / XMinDelta);
    InfoTable.PosYInteger = int16(InfoTable.PositionY / YMinDelta);
    XMin = min(InfoTable.PosXInteger);
    XMax = max(InfoTable.PosXInteger);
    YMin = min(InfoTable.PosYInteger);
    YMax = max(InfoTable.PosYInteger);

    XMinCorr  = (XMin - XMin) + 1;
    XMaxCorr = (XMax - XMin) + 1;
    YMinCorr = (YMin - YMin) + 1;
    YMaxCorr = (YMax - YMin) + 1;

    PositionIm = zeros(YMaxCorr, XMaxCorr, 'logical');
    GroupsIm5DCellArray = cell(YMaxCorr, XMaxCorr);
    
    for i = 1:height(InfoTable)
        InfoTableNow = InfoTable(i,:);
        XNow = InfoTableNow.PosXInteger;
        XNow = (XNow - XMin) + 1;
        YNow = InfoTableNow.PosYInteger;
        YNow = (YNow - YMin) + 1;
        PositionIm(YNow, XNow) = 1;
    end
    
    for x = XMin:XMax
        for y = YMin:YMax
            InfoTableThisField = InfoTable(InfoTable.PosXInteger == x & InfoTable.PosYInteger == y, :);
            GroupsIm5DCellArray{y-YMin+1, x-XMin+1} = InfoTableThisField;
        end
    end
    
    %% Find Objects using connected components logic
    GroupsTable = bwlabeln(PositionIm);
    % it(PositionIm)
end

