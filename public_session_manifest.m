function sessions = public_session_manifest()
%PUBLIC_SESSION_MANIFEST Pseudonymized sessions required by the paper.
% The public analysis deliberately contains no acquisition dates or raw names.

sessionIDs = arrayfun(@(n) sprintf('session_%02d', n), 1:11, ...
    'UniformOutput', false);
areas = {'M1', 'SPL', 'all'};
sessions = cell(numel(sessionIDs) * numel(areas), 2);
row = 0;
for s = 1:numel(sessionIDs)
    for a = 1:numel(areas)
        row = row + 1;
        sessions(row, :) = {sessionIDs{s}, areas{a}};
    end
end
end
