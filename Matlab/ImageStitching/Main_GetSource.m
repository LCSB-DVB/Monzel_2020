%% Authors: Paul Antony and Anna Monzel, 2019-08-20
%% This script stitches the original operetta images using data from "ImagesOperetta" and "Metadata"
%% No image analysis, no further processing
clear
clc
SetupMode = 0; % 1 for creating numeric organoid labels OR 0 for linking the final analysis to human labels

%% Parallel pool control
    delete(gcp('nocreate'))
    myCluster = parcluster;
    Workers = myCluster.NumWorkers;
    % parpool(28) %for HURRICANE
    parpool(Workers) % for MEGATRON
%%
RootSourcePath = 'S:\Operetta\OperettaDB_LCSB';
SaveRootPath = 'S:\HCS_Platform\Path-where-you-want-to-save-the-data'
folders = dir([RootSourcePath, filesep, '*tox*Neurons*']) %find all images that contain "tox" and "neurons"
folders = struct2table(folders)
FullPath = rowfun(@(a,b) strcat(a,filesep,b), folders, 'InputVariables', {'folder','name'})
FullPath = table2cell(FullPath)
folders.FullPath = FullPath
%writetable(folders, [SaveRootPath, filesep, 'Objects.xlsx'])
filesAll = cell(height(folders), 1)
SlideLayout = cell(height(folders), 1)
InfoTable = {} ;
SavePath = {};
for i = 1:height(folders) 
    csvThis = folders{i,'FullPath'}{:};
    files = dirrec(csvThis, '.csv'); %csv file containing the Metadata "metadata.csv" in folder "Metadata"
    txtfiles = dirrec(csvThis, '.txt'); %txt file containing the slide IDs "SlidePreviewID.txt" in folder "Metadata"
    filesAll(i) = files(cellfun(@(x) ~isempty(x), strfind(files, 'metadata')));
    InfoTable{i} = readtable(filesAll{i});
    SavePath{i} = [SaveRootPath, filesep, folders{i,'name'}{:}];
    mkdir(SavePath{i});
    SlideLayout(i) = txtfiles(cellfun(@(x) ~isempty(x), strfind(txtfiles, 'SlidePreviewID')));
    InfoTableToMosaicMat(InfoTable{i}, SavePath{i}, SlideLayout{i}, SetupMode, SaveRootPath)
end


