clear
clc

RootSourcePath = 'S:\Operetta\OperettaDB_LCSB';
SaveRootPath = 'S:\HCS_Platform\Data\AnnaMonzel\3D\TOX_Project_Publication'
folders = dir([RootSourcePath, filesep, '*tox*Neurons*'])
folders = struct2table(folders)
FullPath = rowfun(@(a,b) strcat(a,filesep,b), folders, 'InputVariables', {'folder','name'})
FullPath = table2cell(FullPath)
folders.FullPath = FullPath

filesAll = cell(height(folders), 1)
InfoTable = {} ;
SavePath = ();
for i = 1:height(folders)
    csvThis = folders{i,'FullPath'}{:};
    files = dirrec(csvThis, '.csv');
    filesAll(i) = files(cellfun(@(x) ~isempty(x), strfind(files, 'metadata')));
    InfoTable{i} = readtable(filesAll{i});
    SavePath{i} = [SaveRootPath, filesep, folders{i,'name'}{:}];
    mkdir(SavePath{i});
    
end

InfoTable = readtable(filesAll(i);

SavePath = S:\HCS_Platform\Data\AnnaMonzel\3D\TOX_Project_Publication  folder(i)

