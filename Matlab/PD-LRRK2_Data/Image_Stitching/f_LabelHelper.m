function [ output_args ] = f_LabelHelper(Stencil, SavePath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    %% Preview of the whole slide

%     SizeSingleIm = size(XYMosaicCells{1}{1,1});
%     SizeSingleIm = SizeSingleIm(1:2);
%     RescanGridSize = size(GroupsTable);
%     GreatPreview = zeros(SizeSingleIm(1)*RescanGridSize(1), SizeSingleIm(2)*RescanGridSize(2), 'uint16');
%     ImHeight = SizeSingleIm(1);
%     ImWidth = SizeSingleIm(2);
%     StartRCell = {};
%     StartCCell = {};

    Groups = unique(Stencil(Stencil>0))';


% % 
% %     for g = Groups
% %         g
% %         %StitchedGroupSize = size(GroupPreviews{g});
% %         ZoneNow = Stencil == g;
% %         [R,C] = find(ZoneNow);
% %         StartR = min(R);
% %         StartC = min(C);
% %         StartRPixel = ((StartR-1) * ImHeight) + 1;
% %         %EndRPixel = StartRPixel + (3 * ImHeight) - 1;
% %         EndRPixel = StartRPixel + StitchedGroupSize(1) - 1;
% %         StartCPixel = ((StartC-1) * ImWidth) + 1;
% %         %EndCPixel = StartCPixel + (3 * ImWidth) - 1;
% %         EndCPixel = StartCPixel + StitchedGroupSize(2) - 1;
% %         GreatPreview(StartRPixel:EndRPixel, StartCPixel:EndCPixel) = GroupPreviews{g};
% %         StartRCell{g} = StartRPixel;
% %         StartCCell{g} = StartCPixel;
% %     end

% %     Zoomfactor = 50;
% %     GreatPreviewResized = imresize(imadjust(GreatPreview), 1/Zoomfactor);

Preview = Stencil > 0;
OrganoidLayout = figure
imagesc(Preview)


    for g = Groups
        %GreatPreviewResized = insertText(GreatPreviewResized, [round(StartCCell{g}/Zoomfactor), round(StartRCell{g}/Zoomfactor)], num2str(g), 'FontSize', 12, 'BoxColor', 'red', 'TextColor', 'white');
        
        ZoneNow = Stencil == g;
        [R,C] = find(ZoneNow);
        StartR = min(R);
        StartC = min(C);
        
        %GreatPreviewResized = insertText(Stencil, [StartC, StartR], num2str(g), 'FontSize', 12, 'BoxColor', 'red', 'TextColor', 'white');
        text(StartC, StartR, num2str(g))
        
        
    end
    
    saveas(OrganoidLayout, [SavePath, filesep, 'OrganoidLayout.fig'])
    
end

