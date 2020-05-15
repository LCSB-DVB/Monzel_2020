function [CroppedMosaics, StitchedIms] = f_stitching_opera(XYmosaicCells, XYmosaicContourCell, XPositions, YPositions, ResolutionXY, MaxPixelDrift, PreviewChannel, ShowProgress)
%This function stitches images acquired on the Opera and assumes that the
%sublayout is nonsparse and that each pair of neighboured images presents
%an overlap
%   author: Paul Antony 20160418
%   XYmosaicCells is a cell array containing 1 or more XYmosaicCell
%   outputed from f_load_mosaic_image3D. The order of the XYmosaicCell(s)
%   in XYmosaicCells is the channel order
%   XYmosaicContourCell is a cell array in the same format as XYmosaicCell
%   that contains a preprocessed image which is optimized for cross
%   correlation analysis
%   The variables XPositions, YPositions, and ResolutionXY
%   correspond to the outputs of the function f_load_mosaic_image3D
%   MaxPixelDrift is a positive integer number defining tha maximum
%   expected drift that shall be corrected by image stitching. Recommended
%   value: 10
%   PreviewChannel is an integer defining the channel to use for previews
%   ShowProgress 0 or 1 controls output of previews

% ShowProgress = 1;
% debug = 1;

Zoomfactor = 1;


ImSize = size(XYmosaicCells{1}{1,1});
Rows = 1:size(XYmosaicCells{1}, 1);
Columns = 1:size(XYmosaicCells{1}, 2);

%% Get coordinates of each mosaic element
[X, Y] = meshgrid(XPositions, YPositions);

%% Find pixel steps from image center to image center and derive overlap expressed in pixels 
%XnonZero = X; XnonZero(XnonZero == 0) = []; % Ignore the 0 center
XPositions = sort(XPositions, 'ascend');
XstepLength = XPositions(2)-XPositions(1); % m
StepsX = round((XstepLength*1E6) / (ResolutionXY)); % pixels = m*10E6/(um/pixel)
ExpectedOverlapX = ImSize(2) - StepsX;

%YnonZero = Y; YnonZero(YnonZero == 0) = []; % Ignore the 0 center
%StepsY = round(min(abs(YnonZero)) / ResolutionXY); % pixels = um/(um/pixel)
YPositions = sort(YPositions, 'ascend');
YstepLength = YPositions(2)-YPositions(1);
StepsY = round((YstepLength*1E6) / (ResolutionXY)); % pixels = m*10E6/(um/pixel)
ExpectedOverlapY = ImSize(1) - StepsY;

%% Prepare main output
CroppedMosaics = cell(1, size(XYmosaicCells, 2));
StitchedIms = cell(1, size(XYmosaicCells, 2));

%% Every image has four borders

%From the expected overlap imformation decide how the size and position the boxes for crosscorrelation analysis
BoxOffsetX = ExpectedOverlapX - MaxPixelDrift;
BoxOffsetY = ExpectedOverlapY - MaxPixelDrift;

%Prepare maximum projections for the stitching process
Idx = 0;
for cellcol = 1:size(XYmosaicCells{1}, 2)
    for cellrow = 1:size(XYmosaicCells{1}, 1)
        Idx = Idx+1;
        ImStack = XYmosaicCells{PreviewChannel}{cellrow,cellcol}; % Box preview channel is set here
        Im(:,:,Idx) = max(ImStack,[],4);
        ContourStack = XYmosaicContourCell{cellrow,cellcol};
        ContourIm(:,:,Idx) = max(ContourStack,[],4); %it(ContourIm(:,:,Idx))
    end
end
%vi(Im)

Boxes = {};

for i = 1:size(XYmosaicCells{1},1)*size(XYmosaicCells{1},2)
    ImageID = i;
    ThisIm = Im(:,:,i);
    
    for j = 1:4
        if i == 1 && j == 1
            BoxID = 1;
        else
            BoxID = BoxID + 1;
        end

        Boxes{BoxID,1} = ImageID;
        Boxes{BoxID,2} = BoxID;
        side = j;
        
        switch side
            case 1
                Boxes{BoxID,3} = 'Top';
                Box_X_Start = MaxPixelDrift;
                Box_X_End = size(Im, 2) - MaxPixelDrift;
                Box_Y_Start = MaxPixelDrift;
                Box_Y_End = BoxOffsetY + MaxPixelDrift;
                
            case 2
                Boxes{BoxID,3} = 'Right';
                Box_X_Start = size(Im, 2) - (BoxOffsetX + MaxPixelDrift);
                Box_X_End = size(Im, 2) - MaxPixelDrift;
                Box_Y_Start = MaxPixelDrift;
                Box_Y_End = size(Im, 1) - MaxPixelDrift;
                
            case 3
                Boxes{BoxID,3} = 'Bottom';
                Box_X_Start = MaxPixelDrift;
                Box_X_End = size(Im, 2) - MaxPixelDrift;
                Box_Y_Start = size(Im, 1) - (BoxOffsetY + MaxPixelDrift);
                Box_Y_End = size(Im, 1) - MaxPixelDrift;

            case 4
                Boxes{BoxID,3} = 'Left';
                Box_X_Start = MaxPixelDrift;
                Box_X_End = MaxPixelDrift + BoxOffsetX;
                Box_Y_Start = MaxPixelDrift;
                Box_Y_End = size(Im, 1) - MaxPixelDrift;             
        end
        
        Boxes{BoxID,4} = Box_X_Start;
        Boxes{BoxID,5} = Box_X_End;
        Boxes{BoxID,6} = Box_Y_Start;
        Boxes{BoxID,7} = Box_Y_End;
        


        BoxIm = ThisIm(Box_Y_Start:Box_Y_End, Box_X_Start:Box_X_End);
        
        

        if ShowProgress == 1
            imshow(Im(:,:,i), []); hold on; line([Box_X_Start, Box_X_End, Box_X_End, Box_X_Start, Box_X_Start], [Box_Y_End, Box_Y_End, Box_Y_Start, Box_Y_Start, Box_Y_End], 'LineWidth', 2, 'Color', [1 0 0])
            drawnow        
        end
%         if debug == 1
%             imtool(BoxIm, [])
%         end
        
        Boxes{BoxID,8} = sum(BoxIm(:)) / numel(BoxIm(:));
        Boxes{BoxID,9} = sum(BoxIm(:)) / numel(BoxIm(:)) > 20;
        [r c] = ind2sub(size(XYmosaicCells{1}),i);
        %%%%%%%%%%%%%%%%%%%%%%%%
        %Boxes{BoxID,10} = r; % y
        Boxes{BoxID,10} = 1 + (size(XYmosaicCells{1}, 1) - r); % y
        Boxes{BoxID,11} = c; % x
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        BoxCase = Boxes{BoxID,3};
        
        switch BoxCase
           
            case 'Top'
                
                if Boxes{BoxID,10} ~= min(Rows)
                    Boxes{BoxID,12} = 1;
                    Boxes{BoxID,13} = Rows(find(Rows == Boxes{BoxID,10})-1);
                    Boxes{BoxID,14} = Boxes{BoxID,11};    
                else
                    Boxes{BoxID,12} = 0;
                end
                
            case 'Right'
                
                if Boxes{BoxID,11} ~= max(Columns)
                    Boxes{BoxID,12} = 1;
                    Boxes{BoxID,13} = Boxes{BoxID,10};
                    Boxes{BoxID,14} = Columns(find(Columns == Boxes{BoxID,11})+1);
                else
                    Boxes{BoxID,12} = 0;
                end
                     
            case 'Bottom'
                                
                if Boxes{BoxID,10} ~= max(Rows)
                    Boxes{BoxID,12} = 1;
                    Boxes{BoxID,13} = Rows(find(Rows == Boxes{BoxID,10})+1);
                    Boxes{BoxID,14} = Boxes{BoxID,11};
                else
                    Boxes{BoxID,12} = 0;
                end
                
            case 'Left'
                
                if Boxes{BoxID,11} ~= min(Columns)
                    Boxes{BoxID,12} = 1;
                    Boxes{BoxID,13} = Boxes{BoxID,10};
                    Boxes{BoxID,14} = Columns(find(Columns == Boxes{BoxID,11})-1);
                else
                    Boxes{BoxID,12} = 0;
                end
                         
        end
        
        
    end
end

Boxes = cell2table(Boxes);
Boxes.Properties.VariableNames = {'ImageID', 'BoxID', 'Position', 'Box_X_Start', 'Box_X_End', 'Box_Y_Start', 'Box_Y_End', 'ContourVal', 'ContourOK', 'row', 'column' 'NeighbourAvailable', 'RowNeighbour', 'ColumnNeighbour'};
% if size(Boxes,2) == 14
%     Boxes.Properties.VariableNames = {'ImageID', 'BoxID', 'Position', 'Box_X_Start', 'Box_X_End', 'Box_Y_Start', 'Box_Y_End', 'ContourVal', 'ContourOK', 'row', 'column' 'NeighbourAvailable', 'RowNeighbour', 'ColumnNeighbour'};
% else if size(Boxes,2) == 12
%     Boxes.Properties.VariableNames = {'ImageID', 'BoxID', 'Position', 'Box_X_Start', 'Box_X_End', 'Box_Y_Start', 'Box_Y_End', 'ContourVal', 'ContourOK', 'row', 'column' 'NeighbourAvailable'};
% end
%% Set priorities and walk through the stitching process
BoxesRemaining = Boxes;
BoxesRemaining = sortrows(BoxesRemaining, {'row', 'column'});
SortStatusPrimer = array2table(zeros(height(BoxesRemaining), 1));
SortStatusPrimer.Properties.VariableNames = {'Stitched'};
NeighbourStatusPrimer = array2table(zeros(height(BoxesRemaining), 1));
NeighbourStatusPrimer.Properties.VariableNames = {'NeighbourStitched'};
BoxesRemaining = [BoxesRemaining, SortStatusPrimer, NeighbourStatusPrimer];

%% Initialization
MosaicCoordinatesXY = [];
StitchedIm = uint16(zeros(size(Rows,2)*(size(Im,1)), size(Columns,2)*(size(Im,2)), size(XYmosaicCells{1}{1},4)));
ProgressAll = sortrows(BoxesRemaining, 'ContourVal', 'descend');
StarterIm = sortrows(grpstats(ProgressAll(:,[1,8]), 'ImageID', 'min'), 'min_ContourVal', 'descend');
StarterIm = StarterIm.ImageID(1);
StarterRow = ProgressAll(ProgressAll.ImageID == StarterIm,'row');
StarterRow = StarterRow{1,1};
StarterColumn = ProgressAll(ProgressAll.ImageID == StarterIm,'column');
StarterColumn = StarterColumn{1,1};
ProgressAll(ProgressAll.ImageID == StarterIm, 'Stitched') = num2cell(ones(height(ProgressAll(ProgressAll.ImageID == StarterIm, 'Stitched')),1));
StarterTable = BoxesRemaining(BoxesRemaining.ImageID == StarterIm, :);
template_X_Start = (StarterTable.column(1) * size(XYmosaicCells{1}{1},2)) - (StarterTable.column(1) * ExpectedOverlapX);
template_Y_Start = (StarterTable.row(1) * size(XYmosaicCells{1}{1},1)) - (StarterTable.column(1) * ExpectedOverlapY);

for channel = 1:size(XYmosaicCells,2)
    StitchedIms{channel} = StitchedIm;
    StitchedIms{channel}(template_Y_Start:template_Y_Start+size(Im,1)-1, template_X_Start:template_X_Start+size(Im,2)-1,1:size(StitchedIms{channel},3)) = permute(XYmosaicCells{channel}{StarterIm},[1 2 4 3]);
end

%% Main stitching process
WhileStep = 0;

while height(ProgressAll(ProgressAll.Stitched == 0, :)) > 0
    
    try
        
    WhileStep = WhileStep + 1;
    ProgressMosaic = ProgressAll(ProgressAll.Stitched == 1, :);
    ProgressNext = ProgressMosaic(ProgressMosaic.NeighbourAvailable == 1 & ProgressMosaic.NeighbourStitched ~= 1, :); % Find image of mosaic where stitching continues
    ProgressNextBoxParent = ProgressAll(ProgressAll.row == ProgressNext.RowNeighbour{1} & ProgressAll.column == ProgressNext.ColumnNeighbour{1}, :);
    

    ExtensionDirection = ProgressNext.Position{1};

    switch ExtensionDirection
        case 'Top'
            BoxPosition = 'Bottom';
        case 'Right'
            BoxPosition = 'Left';
        case 'Bottom'
            BoxPosition = 'Top';
        case 'Left'
            BoxPosition = 'Right';
    end
        
    rowNow = ProgressNext.RowNeighbour{1};
    colNow = ProgressNext.ColumnNeighbour{1};
    box = Boxes(Boxes.row == rowNow & Boxes.column == colNow & strcmp(Boxes.Position, BoxPosition), :);    
    box_ParentIm = Im(:,:,ProgressNextBoxParent.ImageID(1));
    box_ParentContourIm = ContourIm(:,:,ProgressNextBoxParent.ImageID(1));
    [StackRow, StackCol] = ind2sub(size(XYmosaicCells{1}), ProgressNextBoxParent.ImageID(1));
    
    if WhileStep == 1
        templateID = StarterIm;
        template = Im(:,:,templateID);
        MosaicCoordinatesXY(StarterIm, 1) = template_X_Start;
        MosaicCoordinatesXY(StarterIm, 2) = template_Y_Start;
    else
        templateID = ProgressNext.ImageID(1);
        template = Im(:,:,templateID);
    end
    
    
    % Define sub-windows for xcorr2 which are searched for maxima
    boxSizeX = round((1/Zoomfactor) * (length([box.Box_X_Start:box.Box_X_End])));
    boxSizeY = round((1/Zoomfactor) * (length([box.Box_Y_Start:box.Box_Y_End])));
    templateSizeX = round((1/Zoomfactor) * (size(template,2)));
    templateSizeY = round((1/Zoomfactor) * (size(template,1)));
    
    crosscorrImSizeX = templateSizeX + boxSizeX - 1;
    crosscorrImSizeY = templateSizeY + boxSizeY - 1;
    
    switch box.Position{:}
        
        case 'Top'
            xcorr2_Max_Window_Rows = [boxSizeY + 0.5*MaxPixelDrift : boxSizeY + 2.5*MaxPixelDrift] 
            xcorr2_Max_Window_Columns = [round(crosscorrImSizeX/2 - MaxPixelDrift):round(crosscorrImSizeX/2 + MaxPixelDrift)]; 
        case 'Right'
            xcorr2_Max_Window_Rows = [round(crosscorrImSizeY/2 - MaxPixelDrift):round(crosscorrImSizeY/2 + MaxPixelDrift)];     
            xcorr2_Max_Window_Columns = [round(crosscorrImSizeX - 2.5*MaxPixelDrift - boxSizeX):round(crosscorrImSizeX - 0.5*MaxPixelDrift - boxSizeX)];         
        case 'Bottom'
             xcorr2_Max_Window_Rows = [round(crosscorrImSizeY - 2.5*MaxPixelDrift - boxSizeY):round(crosscorrImSizeY - 0.5*MaxPixelDrift - boxSizeY)];         
             xcorr2_Max_Window_Columns = [round(crosscorrImSizeX/2 - MaxPixelDrift):round(crosscorrImSizeX/2 + MaxPixelDrift)];
              
        case 'Left'
            xcorr2_Max_Window_Rows = [round(crosscorrImSizeY/2 - MaxPixelDrift):round(crosscorrImSizeY/2 + MaxPixelDrift)];     
            xcorr2_Max_Window_Columns = [round(0.5*MaxPixelDrift + boxSizeX):round(2.5*MaxPixelDrift + boxSizeX)];         
    
    end

    % Stitch and update progress tracking
    disp(['template_' num2str(templateID), '    box_', num2str(box.ImageID), '_', box.Position{:}])
    template_X_Start = MosaicCoordinatesXY(templateID, 1);
    template_Y_Start = MosaicCoordinatesXY(templateID, 2);
    [xStartThisImage, yStartThisImage] = f_BoxToTemplateStitching(box, box_ParentIm, box_ParentContourIm, template, template_X_Start, template_Y_Start, Zoomfactor, xcorr2_Max_Window_Rows, xcorr2_Max_Window_Columns);
    
    for channel = 1:size(XYmosaicCells,2)
        StitchedIms{channel}(yStartThisImage:yStartThisImage+size(box_ParentIm,1)-1, xStartThisImage:xStartThisImage+size(box_ParentIm,2)-1, 1:size(StitchedIms{channel},3)) = permute(XYmosaicCells{channel}{StackRow, StackCol},[1 2 4 3]);
    end
    
    MosaicCoordinatesXY(box.ImageID, 1) = xStartThisImage(1);
    MosaicCoordinatesXY(box.ImageID, 2) = yStartThisImage(1);    
    
    StitchedIm = StitchedIms{PreviewChannel};
    if ShowProgress == 1
        if length(size(StitchedIm)) == 2
            %imtool(imresize(StitchedIm,1/Zoomfactor), [0 100])
            imtool(imresize(StitchedIm,1/Zoomfactor), [])
        elseif length(size(StitchedIm)) == 3
            ShowIm = imresize(StitchedIm,1/Zoomfactor);
            %imtool(max(ShowIm,[],3), [0 100])
            imtool(max(ShowIm,[],3), [])
        end
    end
    
    if WhileStep == 3
    disp('check bug')
    end
    
    ProgressAll(ProgressAll.ImageID == ProgressNextBoxParent.ImageID(1), 'Stitched') = num2cell(ones(height(ProgressAll(ProgressAll.ImageID == ProgressNextBoxParent.ImageID(1), 'Stitched')),1));

    if WhileStep == 1
        RowStarterIDs = cellfun(@(x) x == StarterRow, ProgressAll.RowNeighbour, 'UniformOutput', false);
        EmptyRowStarterIDs = cellfun(@(x) isempty(x), RowStarterIDs);
        RowStarterIDs(EmptyRowStarterIDs) = {false};

        ColStarterIDs = cellfun(@(x) x == StarterColumn, ProgressAll.ColumnNeighbour, 'UniformOutput', false);
        EmptyColStarterIDs = cellfun(@(x) isempty(x), ColStarterIDs);
        ColStarterIDs(EmptyColStarterIDs) = {false};

        StarterNeighbourBoxes = ProgressAll(cell2mat(RowStarterIDs) & cell2mat(ColStarterIDs), :);
        ProgressAll(cell2mat(RowStarterIDs) & cell2mat(ColStarterIDs), 'NeighbourAvailable') = num2cell(zeros(height(StarterNeighbourBoxes), 1));
    end
    

    RowStitchedNeighbourIDs = cellfun(@(x) x == box.row, ProgressAll.RowNeighbour, 'UniformOutput', false);
    EmptyRowStitchedNeighbourIDs = cellfun(@(x) isempty(x), RowStitchedNeighbourIDs);
    RowStitchedNeighbourIDs(EmptyRowStitchedNeighbourIDs) = {false};

    ColStitchedNeighbourIDs = cellfun(@(x) x == box.column, ProgressAll.ColumnNeighbour, 'UniformOutput', false);
    EmptyColStitchedNeighbourIDs = cellfun(@(x) isempty(x), ColStitchedNeighbourIDs);
    ColStitchedNeighbourIDs(EmptyColStitchedNeighbourIDs) = {false};

    StitchedNeighbourBoxes = ProgressAll(cell2mat(RowStitchedNeighbourIDs) & cell2mat(ColStitchedNeighbourIDs), :);

    ProgressAll(cell2mat(RowStitchedNeighbourIDs) & cell2mat(ColStitchedNeighbourIDs), 'NeighbourStitched') = num2cell(ones(height(StitchedNeighbourBoxes), 1));
    
    catch % here the main stitching process ends
        break
    end
    
end

%% Crop and save stitched mosaic

for channel = 1:size(XYmosaicCells,2)
    CroppedMosaics{channel} = StitchedIms{channel}(find(sum(max(StitchedIm,[],3), 2) ~= 0), find(sum(max(StitchedIm,[],3), 1) ~= 0),:); 
end


end

