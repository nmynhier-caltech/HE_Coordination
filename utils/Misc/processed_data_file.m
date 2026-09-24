function filePath = processed_data_file(savepath, area)
%PROCESSED_DATA_FILE Resolve a read-only pseudonymized public handoff file.

    [~, sessionFolder] = fileparts(savepath);
    tokens = regexp(sessionFolder, ...
        '^(session_\d{2})_(M1|SPL|all)_results$', 'tokens', 'once');
    if isempty(tokens)
        error('ProcessedDataFile:InvalidSessionFolder', ...
            'Expected session_XX_<area>_results, received: %s', savepath);
    end
    if ~strcmp(tokens{2}, area)
        error('ProcessedDataFile:AreaMismatch', ...
            'Folder area "%s" does not match requested area "%s".', ...
            tokens{2}, area);
    end
    publicDataRoot = configure_public_data_root();
    filePath = fullfile(publicDataRoot, ...
        sprintf('%s_%s_Processed_data.mat', tokens{1}, area));
end
