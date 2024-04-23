clear all; close all; clc;
rng('default')

fileName = {'Seek_Mixed_final_122.csv';...
    'Share_Mixed_final_125.csv';...
    'SeekMixed_rep_red_Int.csv';...
    'ShareMixed_rep_red_Int.csv'};

%% Colors
% Blue
colorLineI = [36, 36, 255]/256;
colorCII = [217, 217, 255]/256;
% Green
colorLineU = [49, 151, 49]/256;
colorCIU = [217, 236, 217]/256;
% Red
colorLineV = [255, 30, 30]/256;
colorCIV = [255, 217, 217]/256;
% Full
colorLineF= [255, 143, 30] / 256; % Orang
colorCIF = [255, 236, 80] / 256;

%%%%%%%%%%%%%%%%%%%%%
%%% Main Analysis %%%
%%%%%%%%%%%%%%%%%%%%%

for iter = 1:length(fileName)
    %% Load files
    data = importfile(fileName{iter});

    %% Extract subjects
    list_ID = unique(data.ID);
    nSubject = length(list_ID);
    seek_mc = [];
    seek_un = [];
    seek_inst = [];
    All_un = []; All_mc = []; All_inst = []; All_unXmc = []; All_unXinst = []; All_mcXinst = []; All_unXmcXinst = [];
    R2 = [];
    seek_mc0 = [];
    seek_mc1 = [];
    seek_un0 = [];
    seek_un1 = [];
    
    for iSubject = 1:nSubject
        dataSubject = data(data.ID == list_ID(iSubject), :);
        dataSubject.uncertainty = zscore(dataSubject.uncertainty);
        dataSubject.mc = zscore(dataSubject.mc);
        dataSubject.Instr = zscore(dataSubject.Instr);
        dataSubject.seek = zscore(dataSubject.seek);
        
        %% Base model
        mdl = fitlm(dataSubject,'seek ~ mc + uncertainty + Instr');
        BIC(iSubject, 1) = mdl.ModelCriterion.BIC;
        seek_un(iSubject, 1) = mdl.Coefficients.Estimate(2);
        seek_mc(iSubject, 1) = mdl.Coefficients.Estimate(3);
        seek_inst(iSubject, 1) = mdl.Coefficients.Estimate(4);
                        
        %% Interactions model
        mdlAll = fitlm(dataSubject,'seek ~ mc*uncertainty*Instr');
        All_un(iSubject, 1) = mdlAll.Coefficients.Estimate(2);
        All_mc(iSubject, 1) = mdlAll.Coefficients.Estimate(3);
        All_inst(iSubject, 1) = mdlAll.Coefficients.Estimate(4);
        All_unXmc(iSubject, 1) = mdlAll.Coefficients.Estimate(5);
        All_unXinst(iSubject, 1) = mdlAll.Coefficients.Estimate(6);
        All_mcXinst(iSubject, 1) = mdlAll.Coefficients.Estimate(7);
        All_unXmcXinst(iSubject, 1) = mdlAll.Coefficients.Estimate(8);
              
    end
    
    %% Remove subjects
    ind2delete = find(seek_un == 0 & seek_mc == 0 & seek_inst == 0);
    seek_un(ind2delete) = []; seek_mc(ind2delete) = []; seek_inst(ind2delete) = [];
    All_un(ind2delete) = []; All_mc(ind2delete) = []; All_inst(ind2delete) = []; All_unXmc(ind2delete) = [];
    All_unXinst(ind2delete) = []; All_mcXinst(ind2delete) = []; All_unXmcXinst(ind2delete) = [];
    list_ID(ind2delete) = [];
         
    %%%%%%%%%%%%%%%%
    %% Clustering %%
    %%%%%%%%%%%%%%%%
    
    % Example data
    data = [seek_un, seek_mc, seek_inst];
    
    % Initialize silhouette scores array
    evaluationCH = evalclusters(data,'kmeans','CalinskiHarabasz','KList',1:7);
    CH(iter, 1:7) = evaluationCH.CriterionValues;
    
    [idx, ~] = kmeans(data, 3, 'Replicates',10);
    if iter == 1
        code = [2, 3, 1];
    elseif iter == 2
        code = [3, 1, 2];
    elseif iter == 3
        code = [2, 1, 3];
    else
        code = [2, 1, 3];
    end
    mapping = containers.Map(code, 1:numel(code));
    new_idx = arrayfun(@(x) mapping(x), idx);
    p3(1, 1 + 3*(iter - 1):3*iter) = mean(data(idx==code(1), :));
    p3(2, 1 + 3*(iter - 1):3*iter) = mean(data(idx==code(2), :));
    p3(3, 1 + 3*(iter - 1):3*iter) = mean(data(idx==code(3), :));
    p3(4, 1 + 3*(iter - 1):3*iter) = std(data(idx==code(1), :))./sqrt(length(idx));
    p3(5, 1 + 3*(iter - 1):3*iter) = std(data(idx==code(2), :))./sqrt(length(idx));
    p3(6, 1 + 3*(iter - 1):3*iter) = std(data(idx==code(3), :))./sqrt(length(idx));
    [~, p3(7, 1 + 3*(iter - 1):3*iter)] = ttest(data(idx==code(1), :));
    [~, p3(8, 1 + 3*(iter - 1):3*iter)] = ttest(data(idx==code(2), :));
    [~, p3(9, 1 + 3*(iter - 1):3*iter)] = ttest(data(idx==code(3), :));
       
    [~, tests(iter, 1)] = ttest(All_un);
    [~, tests(iter, 2)] = ttest(All_mc);
    [~, tests(iter, 3)] = ttest(All_inst);
    [~, tests(iter, 4)] = ttest(All_unXmc);
    [~, tests(iter, 5)] = ttest(All_unXinst);
    [~, tests(iter, 6)] = ttest(All_mcXinst);
    [~, tests(iter, 7)] = ttest(All_unXmcXinst);
    
    % Assuming `iter` is the iteration number
    filename = ['matrixBF' num2str(iter) '.csv'];
    
    % Assuming `matrixBF` is defined as described
    matrixBF = [All_un, All_mc, All_inst, All_unXmc, All_unXinst, All_mcXinst, All_unXmcXinst];
    
    % Create a table with the matrix and column names
    dataTable = array2table(matrixBF, 'VariableNames', {'All_un', 'All_mc', 'All_inst', 'All_unXmc', 'All_unXinst', 'All_mcXinst', 'All_unXmcXinst'});
    
    % Write the table to a CSV file
    writetable(dataTable, filename);
    
    meanBeta(iter, 1) = mean(All_un);  
    meanBeta(iter, 2) = mean(All_mc);
    meanBeta(iter, 3) = mean(All_inst);
    meanBeta(iter, 4) = mean(All_unXmc);
    meanBeta(iter, 5) = mean(All_unXinst);
    meanBeta(iter, 6) = mean(All_mcXinst);
    meanBeta(iter, 7) = mean(All_unXmcXinst);
    
    %% Test
    [~, ps(iter, 1:3)] = ttest([seek_un, seek_mc, seek_inst]); 
   
    %%%%%%%%%%%%%%%
    %% Save data %%
    %%%%%%%%%%%%%%%
    data2save = table(list_ID, seek_mc, seek_un, seek_inst, new_idx);
    writetable(data2save, ['DataRobust' num2str(iter) '.csv']);
end

%%%%%%%%%%%%%%%%%
%% Across data %%
%%%%%%%%%%%%%%%%%

for iExp = 1:2
    if iExp == 1
        dataSeek = readtable('DataRobust1.csv');
        dataShare = readtable('DataRobust2.csv');
    else
        dataSeek = readtable('DataRobust3.csv');
        dataShare = readtable('DataRobust4.csv');
    end
    
    % Extract relevant columns from dataSeek
    list_ID_seek = dataSeek.list_ID;
    
    % Extract relevant columns from dataShare
    list_ID_share = dataShare.list_ID;
    
    % Find common list_ID values
    common_list_ID = intersect(list_ID_seek, list_ID_share);
    
    % Initialize arrays to store values for correlation calculation
    correlation_values = [];
    
    all_mc_seek = [];
    all_mc_share = [];
    all_un_seek = [];
    all_un_share = [];
    all_inst_seek = [];
    all_inst_share = [];
    all_cluster_seek = [];
    all_cluster_share = [];
    
    % Loop through common_list_ID and calculate correlation for each
    for i = 1:length(common_list_ID)
        id = common_list_ID{i};
        
        % Extract relevant rows for the current list_ID from both tables
        seek_mc_seek_current = dataSeek.seek_mc(ismember(list_ID_seek,id));
        seek_mc_share_current = dataShare.seek_mc(ismember(list_ID_share,id));
        seek_un_seek_current = dataSeek.seek_un(ismember(list_ID_seek,id));
        seek_un_share_current = dataShare.seek_un(ismember(list_ID_share,id));
        seek_inst_seek_current = dataSeek.seek_inst(ismember(list_ID_seek,id));
        seek_inst_share_current = dataShare.seek_inst(ismember(list_ID_share,id));
        seek_cluster_seek_current = dataSeek.new_idx(ismember(list_ID_seek,id));
        seek_cluster_share_current = dataShare.new_idx(ismember(list_ID_share,id));
        
        %
        all_mc_seek = [all_mc_seek; seek_mc_seek_current];
        all_mc_share = [all_mc_share; seek_mc_share_current];
        all_un_seek = [all_un_seek; seek_un_seek_current];
        all_un_share = [all_un_share; seek_un_share_current];
        all_inst_seek = [all_inst_seek; seek_inst_seek_current];
        all_inst_share = [all_inst_share; seek_inst_share_current];
        all_cluster_seek = [all_cluster_seek; seek_cluster_seek_current];
        all_cluster_share = [all_cluster_share; seek_cluster_share_current];
        
    end
    
    % Calculate correlation for the current list_ID
    if iExp == 1
        figure
        correlation_value2 = fitlm(all_mc_share, all_mc_seek, 'RobustOpts', 'huber'); 
        plotCorr2(all_mc_share, all_mc_seek, '', '', colorLineV, colorCIV)
        correlation_value1 = fitlm(all_un_share, all_un_seek,  'RobustOpts', 'huber');
        figure
        plotCorr2(all_un_share, all_un_seek, '', '', colorLineU, colorCIU)
        correlation_value3 = fitlm(all_inst_share, all_inst_seek,  'RobustOpts', 'huber');
        figure
        plotCorr2(all_inst_share, all_inst_seek, '', '', colorLineI, colorCII)
        allcluster1 = mean(all_cluster_seek == all_cluster_share);
    else
        figure
        correlation_value5 = fitlm(all_mc_share, all_mc_seek, 'RobustOpts', 'huber');
        plotCorr2(all_mc_share, all_mc_seek, '', '', colorLineV, colorCIV)
        figure
        correlation_value4 = fitlm(all_un_share, all_un_seek,  'RobustOpts', 'huber');
        plotCorr2(all_un_share, all_un_seek, '', '', colorLineU, colorCIU)
        figure
        correlation_value6 = fitlm(all_inst_share, all_inst_seek,  'RobustOpts', 'huber');
        plotCorr2(all_inst_share, all_inst_seek, '', '', colorLineI, colorCII)
        allcluster2 = mean(all_cluster_seek == all_cluster_share);
    end
end

%%%%%%%%%%%%%%%%%
%% Supp figure %%
%%%%%%%%%%%%%%%%%
figure
for iter = 1:length(fileName)
    % Load files
    data = importfile(fileName{iter});
    
    % Extract subjects
    list_ID = unique(data.ID);
    nSubject = length(list_ID);
    seek_abs_mc = [];
    seek_sign_mc = [];
    seek_mc = [];
    seek_un = [];
    seek_inst = [];
  
    All_un = []; All_mc = []; All_inst = []; All_unXmc = []; All_unXinst = []; All_mcXinst = []; All_unXmcXinst = [];
    R2 = [];
    time = [];
    seek_mc0 = [];
    seek_mc1 = [];
    
    data.uncertainty = zscore(data.uncertainty);
    data.mc = zscore(data.mc);
    data.Instr = zscore(data.Instr);
    data.seek = zscore(data.seek);
    
    subplot(4, 3, 1 + 3*(iter-1))
    dataModel = data;
    mdlAll1 = fitlme(dataModel,'seek ~  mc  + Instr  + (Instr|ID) + (mc|ID)');
    fittedV = fitted(mdlAll1);
    dataModel.seek = dataModel.seek - fittedV;
    mdlAll2 = fitlme(dataModel,'seek ~ uncertainty + (uncertainty|ID)');
    b0 = mdlAll2.Coefficients.Estimate(1); b1 = mdlAll2.Coefficients.Estimate(2);
    scatter(dataModel.seek, fitted(mdlAll2),1,colorCIU)
    mdlAll2.Coefficients.pValue(2);
    
    h1 = lsline
    set(h1, 'LineWidth', 1, 'Color', colorLineU)
    axis equal
    axis([-2.5, 2.5, -2.5,2.5])
    if iter == 1
        title('Uncertainty')
    end
    
    if iter == 1 | iter == 3
        xlabel({'Information Seeking'; '(Observed, Residulas)'}, 'FontSize', 8)
        ylabel({'Information Seeking'; '(Predicted, Residulas)'}, 'FontSize', 8)
    end
    if iter == 2 | iter == 4
        xlabel({'Information Sharing' ;'(Observed, Residulas)'}, 'FontSize', 8)
        ylabel({'Information Sharing'; '(Predicted, Residulas)'}, 'FontSize', 8)
        
    end
    axis square
    subplot(4, 3, 2 + 3*(iter-1))
    dataModel = data;
    mdlAll1 = fitlme(dataModel,'seek ~  uncertainty  + Instr  + (Instr|ID) + (uncertainty|ID)');
    fittedV = fitted(mdlAll1);
    dataModel.seek = dataModel.seek - fittedV;
    mdlAll2 = fitlme(dataModel,'seek ~ mc + (mc|ID)');
    b0 = mdlAll2.Coefficients.Estimate(1); b1 = mdlAll2.Coefficients.Estimate(2);
    scatter(dataModel.seek, fitted(mdlAll2),1,colorCIV)
    h1 = lsline
    set(h1, 'LineWidth', 1, 'Color', colorLineV)
    [~,p]=corr(dataModel.seek, fitted(mdlAll2));
    mdlAll2.Coefficients.pValue(2);
    
    axis([-2.5, 2.5, -2.5,2.5])
    axis square
    
    if iter == 1
        title('Valence')
    end
    if iter == 1 | iter == 3
        xlabel({'Information Seeking' ;'(Observed, Residulas)'}, 'FontSize', 8)
    end
    if iter == 2 | iter == 4
        xlabel({'Information Sharing'; '(Observed, Residulas)'}, 'FontSize', 8)
        
    end
    axis square
     
    subplot(4, 3, 3 + 3*(iter-1))
    dataModel = data;
    mdlAll1 = fitlme(dataModel,'seek ~ mc + uncertainty  + (mc|ID) + (uncertainty|ID)');
    fittedV = fitted(mdlAll1);
    dataModel.seek = dataModel.seek - fittedV;
    mdlAll2 = fitlme(dataModel,'seek ~ Instr  + (Instr|ID)');
    b0 = mdlAll2.Coefficients.Estimate(1); b1 = mdlAll2.Coefficients.Estimate(2);
    scatter(dataModel.seek, fitted(mdlAll2),1,colorCII)
    h1 = lsline
    set(h1, 'LineWidth', 1, 'Color', colorLineI)
    mdlAll2.Coefficients.pValue(2);
    if iter == 1
        title('Instrumentality ')
    end
    
    if iter == 1 | iter == 3
        xlabel({'Information Seeking';'(Observed, Residulas)'}, 'FontSize', 8)
    end
    if iter == 2 | iter == 4
        xlabel({'Information Sharing'; '(Observed, Residulas)'}, 'FontSize', 8)
        
    end
    axis square
    axis equal
    axis([-2.5, 2.5, -2.5,2.5])
    axis square  
end

%%%%%%%%%%%%%%%
%% Functions %%
%%%%%%%%%%%%%%%
function SeekMixedfinal122 = importfile(filename, startRow, endRow)

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

formatSpec = '%f%f%f%f%f%f%f%f%C%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.

%% Create output variable
SeekMixedfinal122 = table(dataArray{1:end-1}, 'VariableNames', {'n_trial','uncertainty','seek','Market_change_1up','mc','Instr','Task1Seek','check','ID','Old1'});
end

%% plotCorr2
function plotCorr2(X,Y, xLabel, yLabel, colorLine, colorCI)

[r, p] = corr(X, Y);
axis([-1,1,-1,1])
xlabel(xLabel, 'FontSize', 12);
ylabel(yLabel, 'FontSize', 12);
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
hold on
mdl = fitlm(X, Y, 'RobustOpts', 'huber');
xl = xlim;
[~, yci] = predict(mdl,[xl(1):.5:xl(2)]');
ciUpper = yci(:,2);
ciLower = yci(:,1);
fill([[xl(1):.5:xl(2)], fliplr([xl(1):.5:xl(2)])], [ciUpper; flip(ciLower)]', colorCI,'EdgeColor','none');
line([-1, 1], [predict(mdl,-1),predict(mdl,1)],'LineWidth', 2.5, 'Color', colorLine)
scatter(X, Y ,25, 'filled' ,'MarkerEdgeColor',colorLine,...
    'MarkerFaceColor',colorLine,...
    'LineWidth',0.5,...
    'MarkerFaceAlpha',1,...
    'MarkerEdgeAlpha',1)

if mdl.Coefficients.pValue(2) < 0.001
    text(-0.8, 0.8, '***', 'FontSize', 25)
elseif mdl.Coefficients.pValue(2) > 0.001 & mdl.Coefficients.pValue(2) < 0.01
    text(-0.8, 0.8, '**', 'FontSize', 25)
elseif mdl.Coefficients.pValue(2) > 0.01 & mdl.Coefficients.pValue(2) < 0.05
    text(-0.8, 0.8, '*', 'FontSize', 25)
end
axis square;
end