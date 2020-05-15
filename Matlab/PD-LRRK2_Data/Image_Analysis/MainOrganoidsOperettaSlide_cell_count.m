clear
clc

SavePath = 'S:\HCS_Platform\Path_where_you_want_to_save_your_data';
mkdir(SavePath)
PreviewSavePath = [SavePath, filesep, 'Previews'];
mkdir(PreviewSavePath)
%channelID = 1; % channel to use for overview

%% Parallel pool control
delete(gcp('nocreate'))
myCluster = parcluster;
Workers = myCluster.NumWorkers;
parpool(Workers) % for MEGATRON

%% Load the data 
% The stitched images are available as .mat files. These images are
% stitched tile scans and include all different planes. To open the images after loading,
% use the vol.m function, which lies in the folder of this script i.e. its workspace. vol(chHoechst, 0, 3000). This
% allows scrolling through the different planes in the stitched image.
    
   Files = dir('S:\HCS_Platform\Path_where_your_data_is_stored_as_stitched_images\*.mat')
   Files = struct2table(Files) 
   Groups = regexp(Files.name,'(.*_.*_.*_.*_.*_\d{1,})', 'tokens') %regular expression of the files
   Groups = cellfun(@(x) x{:}{:}, Groups, 'UniformOutput', false)
   Files.Groups = Groups;
   Groups= unique(Groups);
   GroupPreviews = {};
   ObjectsAll = {};

    for g = 1:numel(Groups)
        
    FileThis = Files(strcmp(Files.Groups, Groups(g)), :)
    PathHoechst = cellfun(@(x) ~isempty(x), strfind(FileThis.name, 'Hoechst'));
    chHoechst = load([FileThis{PathHoechst, 'folder'}{:}, filesep, FileThis{PathHoechst, 'name'}{:}]); chHoechst = chHoechst.Hoechst;
    PathTH = cellfun(@(x) ~isempty(x), strfind(FileThis.name, 'TH'));
    chTH = load([FileThis{PathTH, 'folder'}{:}, filesep, FileThis{PathTH, 'name'}{:}]); chTH = chTH.TH488;
    PathFOXA2 = cellfun(@(x) ~isempty(x), strfind(FileThis.name, 'FOXA2'));
    chFOXA2 = load([FileThis{PathFOXA2, 'folder'}{:}, filesep, FileThis{PathFOXA2, 'name'}{:}]); chFOXA2 = chFOXA2.FOXA2568;
       
%% Image analysis 
        Label = Groups(g);
        try
            ObjectsThisOrganoid = f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, chTH, chHoechst, chFOXA2, PreviewSavePath);
            catch
            Errors{g} = 'Image analysis failed';
            continue % next group g
        end
            ObjectsAll{g} = ObjectsThisOrganoid; %use this if you analyze all sections
    
    Objects = vertcat(ObjectsAll{:});
    save([SavePath, filesep, 'Objects.mat'], 'Objects');
    writetable(Objects, [SavePath, filesep, 'Objects.csv'])
    writetable(Objects, [SavePath, filesep, 'Objects.xlsx'])

    
end
