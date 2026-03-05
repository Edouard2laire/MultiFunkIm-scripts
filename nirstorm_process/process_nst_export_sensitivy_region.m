function varargout = process_nst_export_sensitivy_region( varargin )
% Export the sensitivity of the montage to each region of a definied atlas.
% For each channel, report the average sensitivity to eah region of the
% atlas and export it as tsv file. 
% Authors: Edouard Delaire (2024) 

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Export nirs sensitivity to tsv';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'Multifunkim', 'fNIRS'};
    sProcess.Index       = 7001;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;


    % === Process description
    % === CLUSTERS
    sProcess.options.scouts.Comment = '';
    sProcess.options.scouts.Type    = 'scout';
    sProcess.options.scouts.Value   = {};


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

    sSubject = bst_get('Subject', sInputs.SubjectName);
    sStudy   = bst_get('Study', sInputs.iStudy);
    
    sForward      = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName, 1);
    ChannelMat    = in_bst_channel(sInputs(1).ChannelFile);
    if ndims(sForward.Gain) == 3
        sForward = process_nst_import_head_model('convert_head_model', ChannelMat, sForward, 0);
    end

    if ~strcmp(sSubject.Surface(sSubject.iCortex).FileName, sForward.SurfaceFile)
        bst_error('Headmodel and default cortical surface are not the same');
        return
    end

    sCortex     = in_tess_bst(sForward.SurfaceFile);
    iNIRS       = channel_find(ChannelMat.Channel, 'NIRS');
    sChannel    = ChannelMat.Channel(iNIRS);
    assert(length(sChannel) == size(sForward.Gain,1), 'Headmodel size dont match')

    groups       = {sChannel.Group};
    unique_group = unique(groups);
    nChannel     = length(sChannel) / length(unique_group);

    % average accross multiple wavelength
    gain_matrix  = zeros(nChannel, size(sForward.Gain, 2));
    for iGroup = 1:length(unique_group)
        gain_matrix = gain_matrix + sForward.Gain(strcmp(groups, unique_group{iGroup}), :);
    end
    gain_matrix = gain_matrix ./ length(unique_group);
    max_gain    = max(max(gain_matrix));


    % ROI selection
    ROI     = sProcess.options.scouts.Value;
    iAtlas  = find(strcmp( {sCortex.Atlas.Name},ROI{1}));
    iRois   = cellfun(@(x)find(strcmp( {sCortex.Atlas(iAtlas).Scouts.Label},x)),   ROI{2});
    
    % Threshold : remove all sensitivity lower than -2db
    threshold_value = -2; % in db

    varTypes = ["string",  repmat("double", 1 , length(iRois))];
    varNames = [{'Channel'}, {sCortex.Atlas(iAtlas).Scouts.Label}];
    sz = [nChannel, length(varNames)];
    
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    
    iRow = 1;
    for iPair = 1:nChannel

        gain_channel = squeeze(gain_matrix(iPair,:)); 
        gain_channel(gain_channel <  10^(threshold_value)*max_gain) = 0;
        
        assert(any(gain_channel < 0), 'Found channel with negative gain.')

        for iCluster = 1:length(iRois)
            sROI = sCortex.Atlas(iAtlas).Scouts(iRois(iCluster));
            vertex = sROI.Vertices;

            T{iRow, sROI.Label} = sum(gain_channel(vertex));
        end

        T{iRow,'Channel'} = sForward.pair_names(iPair);
        iRow = iRow + 1;
    end
    
    
    fileName = sProcess.options.outputdir.Value{1};
    writetable(T,fileName, 'FileType','text');

end
