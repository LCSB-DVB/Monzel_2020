clc
clear

%% Load data
dataCell{1}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180508_1716BDay10_FOXA2_TH_Tuj1Batch3\Objects.mat');
dataCell{2}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180509_1716Untr_FOXA2_TH_Tuj1Batch3\Objects.mat');
dataCell{3}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180510_1716A8BDay10_FOXA2_TH_Tuj1Batch3\Objects.mat');
dataCell{4}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180516_1716CDay10_FOXA2_TH_Tuj1Batch3\Objects.mat');
dataCell{5}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180426_1716A8Untr10days_PARP_TH_Tuj1Batch2_A\Objects.mat');
dataCell{6}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180503_1716C10Day_FOXA2_TH_Tuj1Batch2\Objects.mat');
dataCell{7}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180504_1716Veh10Day_FOXA2_TH_Tuj1Batch2\Objects.mat');
dataCell{8}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180428_1716A8Veh10days_FOXA2_TH_Tuj1Batch2\Objects.mat');
dataCell{9}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180514_1716A8VehDay10_FOXA2_TH_Tuj1Batch3\Objects.mat');
dataCell{10}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180502_1716A8Untr_FOXA2_TH_Tuj1Batch2\Objects.mat');
dataCell{11}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180501_1716A8C10days_FOXA2_TH_Tuj1Batch2\Objects.mat');
dataCell{12}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180515_1716A8Untre_FOXA2_TH_Tuj1Batch3\Objects.mat');
dataCell{13}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180426_1716A8Untr10days_PARP_TH_Tuj1Batch2_B\Objects.mat');
dataCell{14}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180430_1716A8B10days_FOXA2_TH_Tuj1Batch2_C\Objects.mat');
dataCell{15}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180504_1716B10Day_FOXA2_TH_Tuj1Batch2\Objects.mat');
dataCell{16}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180415_1716A8_CDay10_PARP_TH_Tuj1\Objects.mat');
dataCell{17}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180420_1716C10days_PARP_TH_Tuj1\Objects.mat');
dataCell{18}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180414_1716A8_BDay10_PARP_TH_Tuj1\Objects.mat');
dataCell{19}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180418_1716Veh10days_PARP_TH_Tuj1\Objects.mat');
dataCell{20}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180506_1716A8C10Day_FOXA2_TH_Tuj1Batch3\Objects.mat');
dataCell{21}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180412_1716A8Veh10days_PARP_TH_Tuj1\Objects.mat');
dataCell{22}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180610_1716VehDay10_FOXA2_TH_Tuj1Batch3\Objects.mat');
dataCell{23}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180411_1716Untr_PARP488_TH568_Tuj1647\Objects.mat');
dataCell{24}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180429_1716A8B10days_FOXA2_TH_Tuj1Batch2_A\Objects.mat');
dataCell{25}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180430_1716A8B10days_FOXA2_TH_Tuj1Batch2_B\Objects.mat');
dataCell{26}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180615_1716VehDay10_FOXA2_TH_Tuj1Batch2_2ndAquisi_A\Objects.mat');
dataCell{27}= load('S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\SB_20180418_1716WT_BDay10_PARP_TH_Tuj1_2nd_A\Objects.mat');

for d = 1:size(dataCell, 2)
    dataThis = dataCell{d};
    dataThistakenApart = dataThis.Objects;
    GoodRows = sum(cellfun(@(x) isempty(x), table2cell(dataThistakenApart)), 2) == 0;
    dataThistakenApart = dataThistakenApart(GoodRows, :);
    dataThistakenApart2 = cellfun(@(x) cellUnpacker(x), table2cell(dataThistakenApart), 'UniformOutput', false);
    dataThistakenApart2 = cell2table(dataThistakenApart2);
    dataThistakenApart2.Properties.VariableNames = dataThistakenApart.Properties.VariableNames;
    data{d} = dataThistakenApart2;
end
 
Objects=vertcat(data{:});
Objects.AreaName

%%Remove undesired rows: Wrong or other tiem points initially included

Objects.AreaName = strrep(Objects.AreaName, '10days', '');
% Objects.AreaName = strrep(Objects.AreaName, 'B', 'PB1mM Rectas 25uM');
% Objects.AreaName = strrep(Objects.AreaName, 'C', 'PB1mM Rectas 2.5uM');

AreaName = unique(Objects.AreaName);
AreaNameSorted = {'1716Untr', '1716Veh','1716B', '1716C', '1716A8Untr', '1716A8Veh','1716A8B', '1716A8C'};

% Objects = Objects (1:3 5:end, :)
Objects = Objects(Objects.NucMaskSum > 5000, :);

Features = Objects.Properties.VariableNames

%% grpstat idea

AreaNamesAll = unique(Objects.AreaName);
ObjectsForGrpStats = Objects(ismember(Objects.AreaName, AreaName), :);
ObjectsForHeatMap = ObjectsForGrpStats;

%% HeatMap

GroupStats = grpstats(ObjectsForHeatMap, 'AreaName', 'mean');
ConditionsNow = strrep({'1716Untr', '1716Veh','1716B', '1716C', '1716A8Untr', '1716A8Veh','1716A8B', '1716A8C'}, '_', '\_');
ReferenceNow = {'1716Untr'};
HeatMapDataNow = GroupStats(:,4:end)

% Remove NaN and Inf by removing the corresponding features
NaNmask = cellfun(@(x) isnan(x), table2cell(HeatMapDataNow));
Infmask = cellfun(@(x) isinf(x), table2cell(HeatMapDataNow));
ToRemove = NaNmask | Infmask;
ToKeep = min(~ToRemove);
HeatMapDataNow = HeatMapDataNow(:, ToKeep);

% Transform data to percent of reference condition
ReferenceRow = find(strcmp(HeatMapDataNow.Properties.RowNames, ReferenceNow));
HeatMapDataArray = table2array(HeatMapDataNow);
for c = 1:size(HeatMapDataArray,2)
    HeatMapDataArray(:,c) = (HeatMapDataArray(:,c) ./ HeatMapDataArray(ReferenceRow,c)) * 100;
end

FeaturesThisHeatmap = strrep(HeatMapDataNow.Properties.VariableNames, '_', '\_');
clustergram(HeatMapDataArray', 'RowLabels', FeaturesThisHeatmap, 'ColumnLabels', HeatMapDataNow.Row, 'DisplayRange', max(HeatMapDataArray(:)), 'LogTrans', true, 'Symmetric', false,'ColorMap', hot)
clustergram(HeatMapDataArray', 'RowLabels', FeaturesThisHeatmap, 'ColumnLabels', HeatMapDataNow.Row, 'DisplayRange', max(HeatMapDataArray(:)), 'LogTrans', true, 'Symmetric', false,'ColorMap', redbluecmap)

%% Create boxplots
AreaNameSorted = {'1716Untr', '1716Veh','1716B', '1716C', '1716A8Untr', '1716A8Veh','1716A8B', '1716A8C'};
s_CreateManuscriptPlotsOrganoids(ObjectsForGrpStats, 'S:\HCS_Platform\Data\SilviaBolognin\Operetta\Organoids_Ibo_DJ1\Graph', AreaNameSorted);
