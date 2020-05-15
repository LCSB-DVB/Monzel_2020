function [Preview] = CreateLabelHelpPreview(Stencil, PreviewSavePath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

 %% Preview

%     SizeSingleIm = size(XYMosaicCells{1}{1,1});
%     SizeSingleIm = SizeSingleIm(1:2);
%     RescanGridSize = size(GroupsTable);
%     GreatPreview = zeros(SizeSingleIm(1)*RescanGridSize(1), SizeSingleIm(2)*RescanGridSize(2), 'uint16');
%     ImHeight = SizeSingleIm(1);
%     ImWidth = SizeSingleIm(2);
%     StartRCell = {};
%     StartCCell = {};
addpath('S:\HCS_Platform\Scripts_Repository\SilviaBolognin\OperettaBeginnings\AddTextToImage')
    Zoomfactor = 50;
    Preview = imresize(Stencil, Zoomfactor, 'nearest');
    Preview = Preview ./ max(Preview(:));
Groups = unique(Stencil);
Groups = Groups(Groups > 0)';

    for g = Groups
%         g
%         ThisGroup = Stencil == g;
%         
%         ThisGroupsRows = find(sum(ThisGroup,2) > 0);
%         ThisGroupsCols = find(sum(ThisGroup,1) > 0);
%         ThisGroup = ThisGroup(ThisGroupsRows, ThisGroupsCols);
        
%         StitchedGroupSize = size(ThisGroup);
        ZoneNow = Stencil == g;
        [R,C] = find(ZoneNow)
        StartR = min(R);
        StartRpixel = 1 + ((StartR-1)*Zoomfactor);
        StartC = min(C);
        StartCpixel = 1 + ((StartC-1)*Zoomfactor);
%         StartRPixel = ((StartR-1) * ImHeight) + 1;
%         %EndRPixel = StartRPixel + (3 * ImHeight) - 1;
%         EndRPixel = StartRPixel + StitchedGroupSize(1) - 1;
%         StartCPixel = ((StartC-1) * ImWidth) + 1;
%         %EndCPixel = StartCPixel + (3 * ImWidth) - 1;
%         EndCPixel = StartCPixel + StitchedGroupSize(2) - 1;
%         GreatPreview(StartRPixel:EndRPixel, StartCPixel:EndCPixel) = GroupPreviews{g};
%         StartRCell{g} = StartRPixel;
%         StartCCell{g} = StartCPixel;
%     it(Preview)
 
        %Preview = insertText(Preview, [StartCpixel, StartRpixel], num2str(g), 'FontSize', 5, 'BoxColor', 'red', 'TextColor', 'white');
        %Preview = AddTextToImage(Preview,'toto',Position,Color,Font,FontSize)
        Position = [StartRpixel, StartCpixel];
        if g == Groups(1)
            %Preview = cat(3, Preview, Preview, Preview);
            Preview = imoverlay(Preview, bwperim(Preview), [1 1 0]);%(3, Preview, Preview, Preview);
        end
        Preview = AddTextToImage(Preview, num2str(g), Position, [1 0 0], 'Arial',40);
        %Preview = insertText(Preview, [901, 101], num2str(g), 'FontSize', 5, 'BoxColor', 'red', 'TextColor', 'white');
%it(Preview)
 

    end


    %GreatPreviewResized = imresize(imadjust(GreatPreview), 1/Zoomfactor);


%     it(Preview)
%     for g = Groups
%         GreatPreviewResized = insertText(GreatPreviewResized, [round(StartCCell{g}/Zoomfactor), round(StartRCell{g}/Zoomfactor)], num2str(g), 'FontSize', 12, 'BoxColor', 'red', 'TextColor', 'white');
%     end
   
it(Preview)
end

