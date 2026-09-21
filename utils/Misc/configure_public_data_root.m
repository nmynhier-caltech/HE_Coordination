function publicDataRoot = configure_public_data_root(varargin)
%CONFIGURE_PUBLIC_DATA_ROOT Store the read-only public-data location.
%   CONFIGURE_PUBLIC_DATA_ROOT(ROOT) configures PROCESSED_DATA_FILE for a
%   pipeline run. CONFIGURE_PUBLIC_DATA_ROOT() returns the configured root.

    persistent configuredRoot

    if nargin > 0
        candidate = varargin{1};
        if ~(ischar(candidate) || isstring(candidate)) || isempty(candidate)
            error('ConfigurePublicDataRoot:InvalidPath', ...
                'publicDataRoot must be a nonempty text path.');
        end
        configuredRoot = char(candidate);
    end

    if isempty(configuredRoot)
        error('ConfigurePublicDataRoot:NotConfigured', ...
            ['No public data root is configured. Run Analysis_Macro with ' ...
             'publicDataRoot set to the extracted CaltechDATA release.']);
    end
    publicDataRoot = configuredRoot;
end
