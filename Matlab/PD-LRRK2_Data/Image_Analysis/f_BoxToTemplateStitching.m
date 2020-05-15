function [xStartThisImage, yStartThisImage] = f_BoxToTemplateStitching(box, box_Parent, box_ParentContourIm, template, template_X_Start, template_Y_Start, ZoomFactor, xcorr2_Max_Window_Rows, xcorr2_Max_Window_Cols)
%This function takes imformation from an image A subregion and a full image
%B within a mosaic, in order to stitch image A into the mosaic containing
%image B
%   box: table describing image A subregion
%   box_Side: description of the position of 'box' within image A. Possible
%   values are {'Top','Right','Bottom','Left'}
%   box_X_Start: column index of box start within image A
%   box_Y_Start: row index of box start within image A
%   box_Parent: image A
%   box_ParentContourIm: contour of image A
%   template: image B
%   template_X_Start: column index of box start within mosaic
%   template_Y_Start: row index of box start within mosaic
%   Zoomfactor: the Zoomfactor which is used to speed up cross correlation
%   xcorr2_Max_Window_Rows: The window rows within xcorr2 output to be searched for a maximum 
%   xcorr2_Max_Window_Cols: The window columns within xcorr2 output to be searched for a maximum 

normalize = 1;

template    = gpuArray(double(imresize(template, 1/ZoomFactor)));
boxIm       = box_ParentContourIm(box.Box_Y_Start:box.Box_Y_End, box.Box_X_Start:box.Box_X_End);
boxIm       = gpuArray(double(imresize(boxIm, 1/ZoomFactor)));
%imtool(gather(template), [])
%imtool(box_Parent, [])
%imtool(gather(boxIm), [])

if normalize == 0
    cc = xcorr2(boxIm,template);
elseif normalize == 1
    try
    cc = normxcorr2(boxIm,template);
    cc = flip(cc, 1);
    cc = flip(cc, 2);
    catch
    cc = xcorr2(boxIm,template);    
    end
end

c = gather(cc);
cbu = c;
c = zeros(size(c));
c(xcorr2_Max_Window_Rows, xcorr2_Max_Window_Cols) = cbu(xcorr2_Max_Window_Rows, xcorr2_Max_Window_Cols);
%c = cbu;
%imtool(cbu, [])
%imtool(c, [])
% figure
% imshow(cbu, []);
% hold on
% line([min(xcorr2_Max_Window_Cols), max(xcorr2_Max_Window_Cols), max(xcorr2_Max_Window_Cols), min(xcorr2_Max_Window_Cols), min(xcorr2_Max_Window_Cols)],...
%      [max(xcorr2_Max_Window_Rows), max(xcorr2_Max_Window_Rows), min(xcorr2_Max_Window_Rows), min(xcorr2_Max_Window_Rows), max(xcorr2_Max_Window_Rows)],'Color', [1 0 0])
% drawnow

% Top left corner of box within template
%[ypeak, xpeak] = find(c==max(c(:)));
%% Adjust to positive values
cAdjusted = zeros(size(c));
cAdjusted_Box = c + abs(min(c(:)));
cAdjusted(xcorr2_Max_Window_Rows, xcorr2_Max_Window_Cols) = cAdjusted_Box(xcorr2_Max_Window_Rows, xcorr2_Max_Window_Cols);
%imtool(c, [])
%imtool(cAdjusted, [])
[ypeak, xpeak] = find(cAdjusted==max(cAdjusted(:)));
GoodContrast = 1;

%% Take expected position in case of absence of acceptable contour information
%if cAdjusted(ypeak,xpeak) < mean(mean(cAdjusted(xcorr2_Max_Window_Rows, xcorr2_Max_Window_Cols))) + (6 * std(std(cAdjusted(xcorr2_Max_Window_Rows, xcorr2_Max_Window_Cols))))
%if cAdjusted(ypeak,xpeak) < mean(mean(cAdjusted_Box(xcorr2_Max_Window_Rows, xcorr2_Max_Window_Cols))) + (10 * std(std(cAdjusted_Box(xcorr2_Max_Window_Rows, xcorr2_Max_Window_Cols))))
%if cAdjusted(ypeak,xpeak) < 3 * mean(mean(cAdjusted(xcorr2_Max_Window_Rows, xcorr2_Max_Window_Cols)))
if cAdjusted(ypeak,xpeak) < (3 * mean(cAdjusted_Box(:)))
    xpeak = round(mean(xcorr2_Max_Window_Cols));
    ypeak = round(mean(xcorr2_Max_Window_Rows));
    GoodContrast = 0;
end

%% Visualize peak detection
show = 0;
if show == 1
    figure
    imshow(cbu, []);
    hold on
    plot(xpeak, ypeak, 'o', 'Color', [0 1 1])
    hold on
    if GoodContrast == 1
        line([min(xcorr2_Max_Window_Cols), max(xcorr2_Max_Window_Cols), max(xcorr2_Max_Window_Cols), min(xcorr2_Max_Window_Cols), min(xcorr2_Max_Window_Cols)],...  
         [max(xcorr2_Max_Window_Rows), max(xcorr2_Max_Window_Rows), min(xcorr2_Max_Window_Rows), min(xcorr2_Max_Window_Rows), max(xcorr2_Max_Window_Rows)],'Color', [0 1 0])
    elseif GoodContrast == 0
        line([min(xcorr2_Max_Window_Cols), max(xcorr2_Max_Window_Cols), max(xcorr2_Max_Window_Cols), min(xcorr2_Max_Window_Cols), min(xcorr2_Max_Window_Cols)],...  
         [max(xcorr2_Max_Window_Rows), max(xcorr2_Max_Window_Rows), min(xcorr2_Max_Window_Rows), min(xcorr2_Max_Window_Rows), max(xcorr2_Max_Window_Rows)],'Color', [1 0 0])
    end

    drawnow
end
%%


xOffset = (size(template, 2) - xpeak) * ZoomFactor;
yOffset = (size(template, 1) - ypeak) * ZoomFactor;

xStartThisImage = template_X_Start - box.Box_X_Start + xOffset;
yStartThisImage = template_Y_Start - box.Box_Y_Start + yOffset;            

%StitchedIm(yStartThisImage:yStartThisImage+size(box_ParentContourIm,1)-1, xStartThisImage:xStartThisImage+size(box_ParentContourIm,2)-1) = box_Parent;

end

