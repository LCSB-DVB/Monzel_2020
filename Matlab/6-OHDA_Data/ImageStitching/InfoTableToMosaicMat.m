function [ output_args ] = InfoTableToMosaicMat( InfoTable, SavePath, SlideLayout, SetupMode, SaveRootPath )

    %% Run mode control
    if SetupMode ==0;
        RunMode = 0;
    else
        RunMode = 1;
    end

    %% Common part
        channelID = 1; % channel to use for overview
        ChannelNames = unique(InfoTable.Channel);
        Channels = length(unique(InfoTable.Channel));
        Planes = unique(InfoTable.Plane)';
        Timepoints = unique(InfoTable.Timepoint)' + 1;
        [GroupsTable, GroupsIm5DCellArray] = FindGroups(InfoTable); % it(GroupsTable)

    %% Setup mode

    if SetupMode == 1

        Preview = CreateLabelHelpPreview(GroupsTable, PreviewSavePath);
        imwrite(Preview, [PreviewSavePath, filesep, 'layout.png'])
        Message = ['The plate Layout has been saved at ', [PreviewSavePath, filesep, 'layout.png'], '. Please save a text file without header, using tab separation, where the first column is the index number as shown in the preview and the second is the area name. Save the text file as SlideLayout_Date.txt in your working directory and set the variable SetupMode to 0 >>> Run'];
        h = msgbox(Message);

    else    
    %% Analysis mode

        % Load annotations
        Layout = readtable(SlideLayout)
        Layout.Properties.VariableNames = {'Idx', 'AreaName'};

        % Load images and organize in an XYC array
        Groups = unique(GroupsTable(GroupsTable > 0))';
        GroupPreviews = {};
        ObjectsAll = {};

        for g = Groups

            XYMosaicCells = {};
            GroupZone = GroupsTable == g;
            [GroupIdxRowVec, GroupIdxColVec] = find(GroupZone); % linear indexes
            Elements = sum(GroupZone(:));
            InfoTablesThisGroup = {};
            for e = 1:Elements % Fields of a given organoid
                for c = 1:Channels
                InfoTableThisField = GroupsIm5DCellArray{GroupIdxRowVec(e), GroupIdxColVec(e)};
                InfoTablesThisGroup{e} = InfoTableThisField;
                InfoTableThisChannel = InfoTableThisField(strcmp(InfoTableThisField.Channel, ChannelNames{c}), :);
                    clear Im4D
                    for t = Timepoints
                        for p = Planes
                            InfoTableThisChannelThisPlane = InfoTableThisChannel(InfoTableThisChannel.Plane == p, :);
                            ImPathThisPlane = InfoTableThisChannelThisPlane.Path{:};   
                            Im4D(:,:,t,p) = imread(ImPathThisPlane); % it(Im4D(:,:,t,p))
                        end
                    end
                   XYMosaicCells{c}{GroupIdxRowVec(e), GroupIdxColVec(e)} = Im4D; % Operetta counterpart of XYmosaicCells for Opera
                end
            end

            InfoTableThisGroup = vertcat(InfoTablesThisGroup{:});

            %% Remove empty cells
            XYMosaicCells = cellfun(@(x) GroupClipper(x),  XYMosaicCells, 'UniformOutput', false);

            %% Stitch
            XYmosaicContourCell = cellfun(@(x) stdfilt(x, ones(3)), XYMosaicCells{1}, 'UniformOutput', false);
            XPositions = unique(InfoTableThisGroup.PositionX); % m
            YPositions = unique(InfoTableThisGroup.PositionY); % m
            ResolutionXY = 675 / 1360; % um per pixel
            MaxPixelDrift = 30;
            PreviewChannel = 1;
            ShowProgress = 0;
            [CroppedMosaic, StitchedIm] = f_stitching_operetta(XYMosaicCells, XYmosaicContourCell, XPositions, YPositions, ResolutionXY, MaxPixelDrift, PreviewChannel, ShowProgress);
            GroupPreviews{g} = max(CroppedMosaic{channelID},[],3); %it(GroupPreviews{g})

            %% Save Channels
            TH488 = CroppedMosaic{1};
            MAP2647 = CroppedMosaic{2};
            Hoechst = CroppedMosaic{3};
            TUJ1568 = CroppedMosaic{4};
            Label = Layout(g,:);
            IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
            save ([SavePath, filesep, IdentityString, 'TH488.mat'], 'TH488')
            save ([SavePath, filesep, IdentityString, 'MAP2647.mat'], 'MAP2647')
            save ([SavePath, filesep, IdentityString, 'Hoechst.mat'], 'Hoechst')
            save ([SavePath, filesep, IdentityString, 'TUJ1568.mat'], 'TUJ1568')

        %% Preview of the whole slide

        SizeSingleIm = size(XYMosaicCells{1}{1,1});
        SizeSingleIm = SizeSingleIm(1:2);
        RescanGridSize = size(GroupsTable);
        GreatPreview = zeros(SizeSingleIm(1)*RescanGridSize(1), SizeSingleIm(2)*RescanGridSize(2), 'uint16');
        ImHeight = SizeSingleIm(1);
        ImWidth = SizeSingleIm(2);
        StartRCell = {};
        StartCCell = {};

        for g = Groups
            g
            StitchedGroupSize = size(GroupPreviews{g});
            ZoneNow = GroupsTable == g;
            [R,C] = find(ZoneNow)
            StartR = min(R);
            StartC = min(C);
            StartRPixel = ((StartR-1) * ImHeight) + 1;
            %EndRPixel = StartRPixel + (3 * ImHeight) - 1;
            EndRPixel = StartRPixel + StitchedGroupSize(1) - 1;
            StartCPixel = ((StartC-1) * ImWidth) + 1;
            %EndCPixel = StartCPixel + (3 * ImWidth) - 1;
            EndCPixel = StartCPixel + StitchedGroupSize(2) - 1;
            GreatPreview(StartRPixel:EndRPixel, StartCPixel:EndCPixel) = GroupPreviews{g};
            StartRCell{g} = StartRPixel;
            StartCCell{g} = StartCPixel;
        end

        Zoomfactor = 50;
        %GreatPreviewResized = imresize(imadjust(GreatPreview), 1/Zoomfactor);
        GreatPreviewResized = imresize(imadjust(GreatPreview, [0 0.02], [0 1]), 1/Zoomfactor);

        for g = Groups
            GreatPreviewResized = insertText(GreatPreviewResized, [round(StartCCell{g}/Zoomfactor), round(StartRCell{g}/Zoomfactor)], num2str(g), 'FontSize', 12, 'BoxColor', 'red', 'TextColor', 'white');
        end

        %imtool(GreatPreview)
        %imtool(GreatPreviewResized)
        imwrite(GreatPreviewResized, [SavePath, filesep, 'GreatPreview.png'])

        % it(GreatPreviewResized)
        % save([SavePath, filesep, 'WorkspaceIncludingObjects.mat'])

    end







end

