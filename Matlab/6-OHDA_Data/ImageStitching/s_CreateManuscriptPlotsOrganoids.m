function [] = s_CreateManuscriptPlotsOrganoids(Data, SavePath, AreaNameSorted);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   Example: s_CreateManuscriptPlots(ObjectsForGrpStats, 'S:\HCS_Platform\Data\SilviaBolognin\General_Script_Data_analysis_high_level');

% AreaNames = unique(Data.AreaName)
AreaNameSorted = {'C5Untr', 'C5Veh','C5PB1mM Rec2.5uM', 'C5PB1mM Rec25uM', 'C5mutUntr', 'C5mutVeh','C5mutPB1mM Rec2.5uM', 'C5mutPB1mM Rec25uM'};
Features = Data.Properties.VariableNames(3:end);

StatsMean = grpstats(Data, {'AreaName'}, {'mean'});
StatsMean = StatsMean(:, 4:end);
StatsMean = StatsMean(ismember(StatsMean.Row, AreaNameSorted),:)

StatsStd = grpstats(Data, {'AreaName'}, {'std'});
StatsStd = StatsStd(:, 4:end);
StatsStd = StatsStd(ismember(StatsStd.Row, AreaNameSorted),:)

%% Boxplots
BoxPlotsDynamic = figure('Position', [486.2000 404.6000 252.8000 261.6000], 'Units', 'pixels')
Features = Features(:, 3:end);
for f = 1:size(Features,2)
    FeatureNow = Features{f}
    boxplot(eval(['Data.' FeatureNow]), Data.AreaName, 'GroupOrder', AreaNameSorted, 'Notch', 'on')
    set(gca, 'XTicklabelRotation', 45, 'XTick', [1 2 3 4 5 6 7 8], 'XLim', [0.5, 9])
    ylabel(FeatureNow)
    pause(1)
    saveas(BoxPlotsDynamic, [SavePath, filesep, FeatureNow, '.fig'])
    saveas(BoxPlotsDynamic, [SavePath, filesep, FeatureNow, '.pdf'])
end
