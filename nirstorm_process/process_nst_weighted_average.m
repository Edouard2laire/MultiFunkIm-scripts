function varargout = process_nst_weighted_average( varargin )
% Comput the weight average of fNIRS channels based on their sensitivity to
% specific regions and export it as tsv file. 
% Authors: Edouard Delaire (2024) 

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Weighted averaged based on sensitivity';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'Multifunkim', 'fNIRS'};
    sProcess.Index       = 7002;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;


    % === Process description

    SelectOptions = {...
        '', ...                            % Filename
        '', ...                            % FileFormat
        'open', ...                        % Dialog type: {open,save}
        'Select weight...', ...     % Window title
        'ImportData', ...                  % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                      % Selection mode: {single,multiple}
        'files', ...                        % Selection mode: {files,dirs,files_and_dirs}
        {{'.tsv'}, '*.tsv'}, ... % Available file formats
        'ResultsIn'};                         % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,AnatIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}
    
    % Option definition
    sProcess.options.weight_file.Comment = 'Sensitivity file:';
    sProcess.options.weight_file.Type    = 'filename';
    sProcess.options.weight_file.Value   = SelectOptions;


    SelectOptions = {...
        '', ...                            % Filename
        '', ...                            % FileFormat
        'save', ...                        % Dialog type: {open,save}
        'Select output folder...', ...     % Window title
        'ExportData', ...                  % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                      % Selection mode: {single,multiple}
        'files', ...                        % Selection mode: {files,dirs,files_and_dirs}
        {{'.tsv'}, '*.tsv'}, ... % Available file formats
        'MriOut'};                         % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,AnatIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}
    
    % Option definition
    sProcess.options.outputdir.Comment = 'Output file:';
    sProcess.options.outputdir.Type    = 'filename';
    sProcess.options.outputdir.Value   = SelectOptions;

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    OutputFiles = {};
        
    weight_matrix = readtable(sProcess.options.weight_file.Value{1},'FileType','text');
    
    sData     = in_bst_data(sInputs.FileName);
    sChannels = in_bst_channel(sInputs.ChannelFile);


    [sDataNIRS, sChannels_data] = process_nst_mbll('filter_bad_channels',sData.F', sChannels, sData.ChannelFlag);

    hb_type = {'HbO', 'HbR', 'HbT'};


    varTypes = [ repmat("double", 1 , size(weight_matrix,2))];
    varNames = [{'Time'}, weight_matrix.Properties.VariableNames(2:end)];
    sz = [size(sDataNIRS,1), length(varNames)];
    
    T_Hbo = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    T_Hbo.Time = sData.Time';

    T_HbR = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    T_HbR.Time = sData.Time';

    T_HbT = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    T_HbT.Time = sData.Time';

    for iRegion = 2:size(weight_matrix,2)
        region_weight = weight_matrix{:, iRegion};
        channels_name = weight_matrix{region_weight > 0 , 1};
        
        if isempty(channels_name)
            % No channel for this region, skip
            continue;
        end

        weight_channels = region_weight (region_weight > 0) ./ sum(region_weight);

        idx_chann = channel_find( sChannels_data.Channel, cellfun(@(x)[ x hb_type{1}], channels_name , 'UniformOutput', false));
        [C, ia,ib] = intersect( cellfun(@(x)[ x hb_type{1}], channels_name , 'UniformOutput', false)  , {sChannels_data.Channel(idx_chann).Name});

        if ~isempty(ia)
            T_Hbo{:,iRegion} = mean( weight_channels(ia)' .* sDataNIRS(:,idx_chann),2);
        end

        idx_chann = channel_find( sChannels_data.Channel, cellfun(@(x)[ x hb_type{2}], channels_name , 'UniformOutput', false));
        [C, ia,ib] = intersect( cellfun(@(x)[ x hb_type{1}], channels_name , 'UniformOutput', false)  , {sChannels_data.Channel(idx_chann).Name});
        if ~isempty(ia)
            T_HbR{:,iRegion} = mean( weight_channels(ia)' .* sDataNIRS(:,idx_chann),2);
        end

        idx_chann = channel_find( sChannels_data.Channel, cellfun(@(x)[ x hb_type{3}], channels_name , 'UniformOutput', false));        
        [C, ia,ib] = intersect( cellfun(@(x)[ x hb_type{1}], channels_name , 'UniformOutput', false)  , {sChannels_data.Channel(idx_chann).Name});
        if ~isempty(ia)
            T_HbT{:,iRegion} = mean( weight_channels(ia)' .* sDataNIRS(:,idx_chann),2);
        end
    end
    
    
    fileName = sProcess.options.outputdir.Value{1};
    writetable(T_Hbo, strrep(fileName, '.tsv','_HbO.tsv'), 'FileType','text');
    writetable(T_HbR, strrep(fileName, '.tsv','_HbR.tsv'), 'FileType','text');
    writetable(T_HbT, strrep(fileName, '.tsv','_HbT.tsv'), 'FileType','text');

end
