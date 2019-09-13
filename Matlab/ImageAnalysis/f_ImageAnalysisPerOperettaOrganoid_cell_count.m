function  [ObjectsThisOrganoid] = FINAL_f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, ch1, ch2, ch3, ch4, PreviewPath);
%% Check images
    % vol(ch1, 0, 3000) % Alexa 488 >>>  	TH
    % vol(ch2, 0, 5000) % Alexa 647 >>> 	MAP2
    % vol(ch3, 0, 3000) % HOECHST 33342 >>> Hoechst 
    % vol(ch4, 0, 3000) % TRITC >>> 		TUJ1

%% Initialize variables
    NucleiMask = [];
    THMask = [];
    MAP2Mask = [];
    Tuj1Mask = [];
    
%% Segment nuclei
    %vol(ch3, 0, 5000)
    ch3BlurSmall = imfilter(double(ch3), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch3BlurSmall)
    ch3BlurBig = imfilter(double(ch3), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch3BlurBig)
    ch3DoG = ch3BlurSmall - ch3BlurBig; %vol(ch3DoG, 0, 200, 'hot')
    NucleiMask = ch3DoG > 20; %75; %vol(NucleiMask)
    NucleiMask = bwareaopen(NucleiMask, 20);%vol(NucleiMask)
    ch3LP = imfilter(ch3, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch3LP, 0, 4000, 'hot')
    NucMaskHigh =  (ch3LP > 1800) .* NucleiMask; % %vol(NucMaskHigh, 0, 1) 800
    NucMaskAlive = NucleiMask & ~NucMaskHigh; % vol(NucMaskAlive)
      
%% TH segmentation (ch1) 
   %vol(ch1, 0, 3000)
    TH_FT = zeros(size(ch1), 'double');
parfor p=1:size(ch1, 3)
    TH_FT(:,:,p) = f_LPF_by_FFT(ch1(:,:,p), 'Butterworth', [7,1], 0);
end
THMask = TH_FT > 0.003; %modify this for different threshold
THMask = bwareaopen(THMask, 1000);
THMask =medfilt3(THMask);
%vol(THMask)   
%vol(TH_FT)    
 
%% TH skeleton3D EPFL
    disp('Start skel')
    tic
    skelTH = Skeleton3D(THMask);
    toc
    disp('Skel done')
%     vol(skelTH, 0, 1)
    [AdjacencyMatrixTH, nodeTH, linkTH] = Skel2Graph3D(skelTH,0);                       
    %imtool(AdjacencyMatrixTH, [])
    NodeTH = zeros(size(THMask), 'uint8');
    NodeIdxs = vertcat(nodeTH(:).idx);
    NodeTH(NodeIdxs) = 1;
%     vol(NodeTH)    
    if size(NodeIdxs, 1) == 0
        return
    end
    NodeTHPreview = uint8(skelTH) + NodeTH + uint8(THMask); 
    NodeTHPreview2D = max(NodeTHPreview, [], 3);
    %it(NodeTHPreview2D)
%   vol(NodeTHPreview, 0, 3, 'jet')    
    NodeDegreeVectorTH = sum(AdjacencyMatrixTH, 1);
    
    ZeroNodeExplanationNeeded = 0;
    if ZeroNodeExplanationNeeded
        ZeroNodes = find(NodeDegreeVectorTH == 0);
        ZeroNodesLinIdx = vertcat(nodeTH(ZeroNodes).idx);
        ZeroNodeMask = zeros(size(THMaskClipped), 'uint8');
        ZeroNodeMask(ZeroNodesLinIdx) = 1; %vol(ZeroNodeMask)
        NodePreviewZeroCase = uint8(skelTH) + NodeMaskTH + 10*uint8(ZeroNodeMask) + uint8(THMask);
    end  
    
%% TH Fragmentation    
    % Define structuring element for surface detection
    Conn6 = strel('sphere', 1); % 6 connectivity
    % Detect surface
    SurfaceTH = THMask & ~(imerode(THMask, Conn6));
%vol(SurfaceTH)
   
%% Tuj1 Segmentation (ch4)
        %vol(ch4, 1, 3000)
    Tuj1_FT = zeros(size(ch4), 'double');
    for p=1:size(ch4, 3)
        Tuj1_FT(:,:,p) = f_LPF_by_FFT(ch4(:,:,p), 'Butterworth', [7,1], 0);
    end
    Tuj1Mask = Tuj1_FT > 0.0015;%0.001;
    Tuj1Mask = bwareaopen(Tuj1Mask, 500);
    Tuj1Mask = medfilt3(Tuj1Mask);
        %vol(Tuj1Mask)    

%% MAP2 Segmentation (ch2)
       %vol(ch2, 1, 2000)
    MAP2_FT = zeros(size(ch2), 'double');
    for p=1:size(ch2, 3)
        MAP2_FT(:,:,p) = f_LPF_by_FFT(ch2(:,:,p), 'Butterworth', [7,1], 0);
    end
    MAP2Mask = MAP2_FT > 0.0015;%0.001;
    MAP2Mask = bwareaopen(MAP2Mask, 500);
    MAP2Mask = medfilt3(MAP2Mask);
        %vol(MAP2Mask)    
%% Perinuclear Volume    
    %vol(NucleiMask)
    NucleiMaskSingleCells = f_RemoveBigObjects (NucleiMask, 10000); 
    NucDil = imdilate(imdilate(NucleiMaskSingleCells, strel('disk', 4)), strel('sphere',1));
    NucPerim = logical(NucDil) & ~logical(NucleiMaskSingleCells);
    %vol(NucPerim)
    THMaskinNucPerim = THMask & NucPerim;% vol(THMaskinNucPerim)
    
%% Percent TH pos
        %split perinuc
    D = bwdist(NucleiMaskSingleCells);
        %vol(D, 0, 20, 'hot')
        %it(D(:,:,1))
    disp('start watershed')
    tic
    W = watershed(D);
    toc
    disp('watershed done')
        %vol(W)
    NucPerimStencil = uint16(W) .* uint16(imreconstruct(logical(imdilate(NucPerim, strel('disk', 1))), logical(NucleiMaskSingleCells))); % This line was causing the error 20171207 % Function imreconstruct expected MARKER and MASK to have the same class.
        %vol(NucPerimStencil)
        %vol(NucPerim)
        %vol(NucleiMaskSingleCells)
        %toto = imreconstruct(logical(NucPerim), logical(NucleiMaskSingleCells));
       
    PeriNucMask = logical(NucPerimStencil);
    PeriNucMask = bwareaopen(PeriNucMask, 500);
        %vol(PeriNucMask)
    PerinucLM = bwlabeln(PeriNucMask);%vol(PerinucLM); vol(uint16(PeriNucMask) .* uint16(THMask), 0,1); vol(THMask +2*PeriNucMask)
    PeriNucObjects = regionprops('table', PerinucLM, double(THMask), 'PixelValues');
    THproportions = rowfun(@(x) sum(x{:})/length(x{:}), PeriNucObjects, 'InputVariables', 'PixelValues');
    THPos = array2table(table2array(THproportions) > 0.01);
    THPos.Properties.VariableNames(end) = {'THpos'};
    PeriNucObjects = [PeriNucObjects, THproportions, THPos];
    PeriNucObjects.Properties.VariableNames(end-1) = {'THproportion'};
    PeriNucObjectsCompact = PeriNucObjects(:, {'THproportion','THpos'});
    THPercent = (sum(PeriNucObjectsCompact.THpos)/height(PeriNucObjectsCompact))*100;
     
%% Previews 
    
    %Scalebar  
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
        %it(BarMask)

    PreviewTH = imoverlay2(imadjust(max(ch1,[],3),[0 0.03]), bwperim(max(THMask,[],3)), [1 0 0]);
    PreviewTH = imoverlay2(PreviewTH, BarMask, [1 1 1]);
        %imtool(PreviewTH)
    PreviewHoechst = imoverlay2(imadjust(max(ch3,[],3),[0 0.04]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
        %imtool(PreviewHoechst)
    PreviewHoechstAlive = imoverlay2(imadjust(max(ch3,[],3),[0 0.04]), bwperim(max(NucMaskAlive,[],3)), [1 0 0]);
    PreviewHoechstAlive = imoverlay2(PreviewHoechstAlive, BarMask, [1 1 1]);
       % imtool(PreviewHoechstAlive)

    PreviewTuj1 = imoverlay2(imadjust(max(ch4, [], 3), [0 0.03]), bwperim(max(Tuj1Mask,[],3)), [0 0 1]);
    PreviewTuj1 = imoverlay2(PreviewTuj1, BarMask, [1 1 1]);
       % imtool(PreviewTuj1)
    
    PreviewMAP2 = imoverlay2(imadjust(max(ch2, [], 3), [0 0.025]), bwperim(max(MAP2Mask,[],3)), [0 0 1]);
    PreviewMAP2 = imoverlay2(PreviewMAP2, BarMask, [1 1 1]);
        %imtool(PreviewMAP2)
  
   PreviewPeriNucMask = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(PeriNucMask,[],3)), [1 0 0]);
   PreviewPeriNucMask = imoverlay2(PreviewPeriNucMask, BarMask, [1 1 1]);
  
    IdentityString = ['Preview', '_', Label{1}];
    imwrite(PreviewTH, [PreviewPath, filesep, IdentityString, '_', 'TH', '.png'])
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])
    imwrite(PreviewTuj1, [PreviewPath, filesep, IdentityString, '_', 'Tuj1', '.png'])
    imwrite(PreviewMAP2, [PreviewPath, filesep, IdentityString, '_', 'MAP2', '.png'])
    imwrite(PreviewPeriNucMask, [PreviewPath, filesep, IdentityString, '_', 'PeriNucMask', '.png'])
    imwrite(PreviewHoechstAlive, [PreviewPath, filesep, IdentityString, '_', 'NucAlive', '.png'])
    
%% Feature extraction
    
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = {Label{1}};
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.Tuj1MaskSum = sum(Tuj1Mask(:));
    ObjectsThisOrganoid.MAP2MaskSum = sum(MAP2Mask(:));
    ObjectsThisOrganoid.THMaskSum = sum(THMask(:));
    ObjectsThisOrganoid.Tuj1ByNuc = sum(Tuj1Mask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.MAP2ByNuc = sum(MAP2Mask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.THByNuc = sum(THMask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.NucMaskHigh = sum(NucMaskHigh(:));
    ObjectsThisOrganoid.NucMaskAlive = sum(NucMaskAlive(:));
    ObjectsThisOrganoid.THFragmentation = sum(SurfaceTH(:)) / sum(THMask(:));
    ObjectsThisOrganoid.SkelTH = sum(skelTH(:));
    ObjectsThisOrganoid.Nodes = size(nodeTH, 2);
    ObjectsThisOrganoid.Links = size(linkTH, 2);
    ObjectsThisOrganoid.THPercent = THPercent;
    
end

