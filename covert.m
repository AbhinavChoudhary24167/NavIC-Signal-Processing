% convert_all_mat_to_unique_csv.m
% Converts multiple .mat files to CSV with unique filenames for each variable

% List of input .mat files
matFiles = {'in_data2.mat', 'in_data3.mat', 'in_data4.mat', 'in_data5.mat'};

% Loop through each .mat file
for k = 1:length(matFiles)
    matFile = matFiles{k};
    
    % Load the .mat file
    fprintf('Processing "%s"...\n', matFile);
    data = load(matFile);
    
    % Extract variable names
    varNames = fieldnames(data);
    
    % Export each variable to CSV with a unique name
    for i = 1:length(varNames)
        varName = varNames{i};
        signal = data.(varName);
        
        % Generate a unique file name: e.g., in_data1_varname.csv
        outputFileName = sprintf('%s_%s.csv', matFile(1:end-4), varName);
        
        % Save to CSV
        writematrix(signal, outputFileName);
        fprintf('  -> Saved "%s"\n', outputFileName);
    end
end

fprintf('✅ All .mat files converted to uniquely named .csv files.\n');
