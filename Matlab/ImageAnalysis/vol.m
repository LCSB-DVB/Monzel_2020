function vol(Im, varargin)
%Spartanic viewer for 3D images featuring scaling of size and intensity, scrolling through planes, and pixel information
%   Author: Paul Antony 2016/11/24
%   Im: imge in 3D
%   The 2nd argument is optional and represents the minimum intensity to display (requires 3rd argument)
%   The 3rd argument is optional and represents the maximum intensity to display (requires 2nd argument)
%   The 4th argument is optional and can provide colormaps such as 'gray', 'jet', or 'hot'
%   Example: vol(ch1,0,200,'jet')


%% initialize display
ScaleFactor = 1000 / size(Im,1);
if ScaleFactor > 5
    ScaleFactor = floor(ScaleFactor);
end
fig_handle = figure('Position',[100, 100, ScaleFactor*size(Im,2), ScaleFactor*size(Im,1)])
set(fig_handle, 'WindowScrollWheelFcn', @wheel);

CurrentSlice = 1;
ha1 = axes('Units', 'norm', 'Position',[0,0,1,1]);		   
updatePlot()
set(findall(fig_handle, '-property', 'Units'), 'Units', 'Normalized')



    %% callback for mouse scrolling
    function wheel(hObject, callbackdata)
        if callbackdata.VerticalScrollCount < 0 & CurrentSlice < size(Im,3)
            CurrentSlice = CurrentSlice + 1;
            updatePlot()            
        elseif callbackdata.VerticalScrollCount > 0 & CurrentSlice > 1
            CurrentSlice = CurrentSlice - 1;
            updatePlot()           
        else
            return            
        end
    end

        
    function updatePlot()
        axes(ha1)
        if length(varargin) == 0
            imshow(Im(:, :, CurrentSlice),[])
            colormap(ha1, 'gray')
        elseif length(varargin) == 3
            imshow(Im(:, :, CurrentSlice),[varargin{1}, varargin{2}])
            colormap(ha1, varargin{3})
        else
            imshow(Im(:, :, CurrentSlice),[varargin{1}, varargin{2}])
        end
        impixelinfo
    end

end
