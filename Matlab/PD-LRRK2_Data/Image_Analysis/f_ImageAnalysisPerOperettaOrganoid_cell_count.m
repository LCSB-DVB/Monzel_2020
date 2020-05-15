function  [ObjectsThisOrganoid] = FINAL_f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, ch1, ch2, ch3, PreviewPath);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % vol(ch1, 0, 500) % Alexa 488 >>> TH
    % vol(ch2, 0, 10000) % HOECHST 33342 >>> Hoechst imtool(max(ch2, [], 3))
    % vol(ch3, 0, 1000) % TRITC >>> FOXA2

    %% Initialize variables
    NucleiMask = [];
    THMask = [];
    FOXA2Mask = [];
    
    %% Segment nuclei
    %vol(ch2, 0, 3000)
    ch2BlurSmall = imfilter(double(ch2), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch2BlurSmall)
    ch2BlurBig = imfilter(double(ch2), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch2BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch2DoG = ch2BlurSmall - ch2BlurBig; %vol(ch2DoG, 0, 200, 'hot')
    NucleiMask = ch2DoG > 75; %vol(NucleiMask)
    NucleiMask = bwareaopen(NucleiMask, 20);%vol(NucleiMask)
    ch2LP = imfilter(ch2, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch2LP, 0, 10000, 'hot')
    NucMaskHigh =  (ch2LP > 4000) .* NucleiMask; %vol(NucMaskHigh, 0, 1) %%before3500
    NucMaskAlive = NucleiMask & ~NucMaskHigh; % vol(NucMaskAlive)
     
    %% TH (ch1)
    ch1MedFilt = []; 
    SizeZ = size(ch1, 3);
    parfor p = 1:SizeZ
        ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p));
    end
    %vol(ch1MedFilt, 0, 2000, 'hot')
    THMask = ch1MedFilt > 400; % vol(THMask)
    THDoG = imfilter(ch1, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch1, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(THDoG, 0, 300, 'hot')
    THDoGMask = THDoG > 150;
    %vol(THDoGMask)
    THMask = THMask & THDoGMask;
    %vol(THMask, 0, 1)
    %it(max(THMask, [], 3))
    %THMask = THMask & ~NucleiMask; %%%%%% IF THIS IS ON WILL CREATE ERROR IN THE %TH
    THMask = bwareaopen(THMask, 200);
    
    %% TH skeleton3D EPFL
    disp('Start skel')
    tic
    skelTH = Skeleton3D(THMask);
    toc
    disp('Skel done')
    %vol(skelTH, 0, 1)
    [AdjacencyMatrixTH, nodeTH, linkTH] = Skel2Graph3D(skelTH,0);                       
    %imtool(AdjacencyMatrixTH, [])
    NodeTH = zeros(size(THMask), 'uint8');
    NodeIdxs = vertcat(nodeTH(:).idx);
    NodeTH(NodeIdxs) = 1;
    %vol(NodeTH)    
    if size(NodeIdxs, 1) == 0
        return
    end
    NodeTHPreview = uint8(skelTH) + NodeTH + uint8(THMask); 
    NodeTHPreview2D = max(NodeTHPreview, [], 3);
    %it(NodeTHPreview2D)
    %vol(NodeTHPreview, 0, 3, 'jet')    
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
    
    %% FOXA2 (ch3)
    
    %vol(ch3, 0, 1000)
    ch3MedFilt = [];
    parfor p = 1:size(ch3, 3)
        ch3MedFilt(:,:,p) = medfilt2(ch3(:,:,p));
    end
    %vol(ch3MedFilt, 0, 1500, 'hot')    
    FOXA2Mask = ch3MedFilt > 350; %FOXA2Mask = ImFOXA2MedFilt > 150; BEFORE 800
%     FOXA2Mask = FOXA2Mask & NucleiMask;
    FOXA2Mask = bwareaopen(FOXA2Mask, 50); %vol(FOXA2Mask) 
    
%     PreviewFOXA2 = imoverlay2(imadjust(max(ch3, [], 3), [0 0.02]), bwperim(max(FOXA2Mask,[],3)), [0 0 1]);
%     PreviewFOXA2 = imoverlay2(PreviewFOXA2, BarMask, [1 1 1]);
    %imtool(PreviewFOXA2)

    %% Perinuclear Volume (to detect amount of cells positive for TH-FOXA2)
    
    %vol(NucleiMask)
    NucleiMaskSingleCells = f_RemoveBigObjects (NucleiMask, 10000); 
    NucDil = imdilate(imdilate(NucleiMaskSingleCells, strel('disk', 4)), strel('sphere',1));
    NucPerim = logical(NucDil) & ~logical(NucleiMaskSingleCells);
    %vol(NucPerim)
    THMaskinNucPerim = THMask & NucPerim;% vol(THMaskinNucPerim)
    FOXA2Dil = imdilate(imdilate(FOXA2Mask, strel('disk', 4)), strel('sphere',1));
    FOXA2Perim = FOXA2Dil & ~NucleiMask;
    Duplex = FOXA2Perim & THMask;
    %vol(Duplex)vol(Duplex+2*THMask)
    
    %% FOXA2 Mask in non-TH
    FOXA2MaskNonTH = FOXA2Perim & ~THMask;
    %vol(FOXA2MaskNonTH,0,1)        
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
    PeriNucObjects.Properties.VariableNames(end-1) = {'THproportion'};%{'PixelValues', 'THproportion', 'THpos'};
    PeriNucObjectsCompact = PeriNucObjects(:, {'THproportion','THpos'});
    THPercent = (sum(PeriNucObjectsCompact.THpos)/height(PeriNucObjectsCompact))*100;
     
    %% Previews 
    
    % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)

    PreviewTH = imoverlay2(imadjust(max(ch1,[],3),[0 0.05]), bwperim(max(THMask,[],3)), [1 0 0]);
    PreviewTH = imoverlay2(PreviewTH, BarMask, [1 1 1]);
    %imtool(PreviewTH)
    
    PreviewHoechst = imoverlay2(imadjust(max(ch2,[],3),[0 0.075]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
    % imtool(PreviewHoechst)

    PreviewFOXA2 = imoverlay2(imadjust(max(ch3, [], 3), [0 0.02]), bwperim(max(FOXA2Mask,[],3)), [0 0 1]);
    PreviewFOXA2 = imoverlay2(PreviewFOXA2, BarMask, [1 1 1]);
    %imtool(PreviewFOXA2)
    
    PreviewPeriNucMask = imoverlay2(imadjust(max(ch2,[],3),[0 0.075]), bwperim(max(PeriNucMask,[],3)), [1 0 0]);
    PreviewPeriNucMask = imoverlay2(PreviewPeriNucMask, BarMask, [1 1 1]);
    
    PreviewFOXA2MaskNonTH = imoverlay2(imadjust(max(ch3, [], 3), [0 0.02]), bwperim(max(FOXA2MaskNonTH,[],3)), [0 0 1]);
    PreviewFOXA2MaskNonTH = imoverlay2(PreviewFOXA2MaskNonTH, BarMask, [1 1 1]);
    %imtool(PreviewFOXA2MaskNonTH)
    
    PreviewFOXA2MaskTH = imoverlay2(imadjust(max(ch3, [], 3), [0 0.02]), bwperim(max(Duplex,[],3)), [0 0 1]);
    PreviewFOXA2MaskTH = imoverlay2(PreviewFOXA2MaskTH, BarMask, [1 1 1]);
    %imtool(PreviewFOXA2MaskTH)
   
    IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
    imwrite(PreviewTH, [PreviewPath, filesep, IdentityString, '_', 'TH', '.png'])
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])
    imwrite(PreviewFOXA2, [PreviewPath, filesep, IdentityString, '_', 'FOXA2', '.png'])
    imwrite(PreviewPeriNucMask, [PreviewPath, filesep, IdentityString, '_', 'PeriNucMask', '.png'])
    imwrite(PreviewFOXA2MaskNonTH, [PreviewPath, filesep, IdentityString, '_', 'FOXA2MaskNonTH', '.png'])
    imwrite(PreviewFOXA2MaskTH, [PreviewPath, filesep, IdentityString, '_', 'FOXA2MaskTH', '.png'])
    
    
    %% Feature extraction
    
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = {Label.Idx};
    ObjectsThisOrganoid.AreaName = {Label.AreaName};
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.FOXA2MaskSum = sum(FOXA2Mask(:));
    ObjectsThisOrganoid.NonTH_FOXA2 = sum(FOXA2MaskNonTH(:));
    ObjectsThisOrganoid.TH_FOXA2 = sum(Duplex(:));
    ObjectsThisOrganoid.THMaskSum = sum(THMask(:));
    ObjectsThisOrganoid.FOXA2ByNuc = sum(FOXA2Mask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.THByNuc = sum(THMask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.NucMaskHigh = sum(NucMaskHigh(:));
    ObjectsThisOrganoid.NonTH_FOXA2byNucPerim = sum(FOXA2MaskNonTH(:)) / sum(NucPerim(:));
    ObjectsThisOrganoid.TH_FOXA2byNucPerim = sum(Duplex(:)) / sum(NucPerim(:));
    ObjectsThisOrganoid.THFragmentation = sum(SurfaceTH(:)) / sum(THMask(:));
    ObjectsThisOrganoid.SkelTH = sum(skelTH(:));
    ObjectsThisOrganoid.Nodes = size(nodeTH, 2);
    ObjectsThisOrganoid.Links = size(linkTH, 2);
    ObjectsThisOrganoid.THPercent = THPercent;
    

end

