
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
SaveRootPath = 'S:\PathWhereYouWantToSaveTheData'
folders = dir([RootSourcePath, filesep, '*LS*BILL*']) %
folders = struct2table(folders)
folders = folders(strcmp(folders.name,'LS_20171003_BILLWT-MutP4-5-6'),:) % | ... 
% %     strcmp(folders.name,'LS_20180625_DAN2018_11cd') | ...%
% %     strcmp(folders.name,'LS_20180621_DAN2018_11e') | ...
% %     strcmp(folders.name,'LS_20180621_DAN2018_11f') | ...
% %     strcmp(folders.name,'LS_20180608_DAN2018_12') | ...
% %     strcmp(folders.name,'LS_20180608_DAN2018_13') | ...
% %     strcmp(folders.name,'LS_20180612_DAN2018_14') | ...
% %     strcmp(folders.name,'LS_20180621_DAN2018_15') | ...
% %     strcmp(folders.name,'LS_20180621_DAN2018_15b') | ...
% %     strcmp(folders.name,'LS_20180620_DAN2018_16') |...
% %     strcmp(folders.name,'LS_20180620_DAN2018_17') |...
% %     strcmp(folders.name,'LS_20180621_DAN2018_18') |...
% %     strcmp(folders.name,'LS_20180620_DAN2018_19') |...
% %     strcmp(folders.name,'LS_20180625_DAN2018_20') ,:)

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
 %   ifmultiple = find(~cellfun(@isempty,strfind(files,'723f79bb-00ab-4abd-b099-a5ab13f4ffa2')))
 %   files = files(ifmultiple(:))
    txtfiles = dirrec(csvThis, '.txt'); %txt file containing the slide IDs "SlidePreviewID.txt" in folder "Metadata"
    filesAll(i) = files(cellfun(@(x) ~isempty(x), strfind(files, 'metadata')));
    InfoTable{i} = readtable(filesAll{i});
    SavePath{i} = [SaveRootPath, filesep, folders{i,'name'}{:}];
    mkdir(SavePath{i});
    SlideLayout(i) = txtfiles(cellfun(@(x) ~isempty(x), strfind(txtfiles, 'SlidePreviewID')));
    InfoTableToMosaicMat(InfoTable{i}, SavePath{i}, SlideLayout{i}, SetupMode, SaveRootPath)
end


