function createMnxRefReactionDbFile(mnxref_dir)
%% createMnxRefReactionDbFile(mnxref_dir)
% Create a reaction template file from the MNXref database.
% Input:
%       char mnxref_dir:            directory where the reac_prop.tsv file
%                                   is located
% Creates a file:
% * reaction_MNXref_balanced.lst
% Format:
% Reaction ID: TAB general reaction equation

if isfile(fullfile(mnxref_dir, 'reac_prop.tsv'))
    % Columns (copied from MetaNetX):
    % 1) The identifier of the reaction in the MNXref namespace
    % 2) Equation of the reaction in the MNXref namespace (compartmentalized and undirected)
    % 3) The original best resource from where this reaction comes from
    % 4) The EC(s) associated to this reaction
    % 5) Is the equation balanced with respect to elemental composition and charge
    % 6) Is this a transport reaction
    reac_prop_table = readtable(fullfile(mnxref_dir, 'reac_prop.tsv'),...
        'Delimiter', '\t', ...
        'FileType', 'text', ...
        'CommentStyle', '#', ...
        'ReadVariableNames', false ...
        );
else
    error('File not found: %s', fullfile(mnxref_dir, 'reac_prop.tsv'))
end

% only consider reactions that are balanced:
keep_idx = ismember(reac_prop_table.(5), 'B');
reaction_ids = reac_prop_table.(1)(keep_idx);
reaction_formulas = reac_prop_table.(2)(keep_idx);
clear reac_prop_table

% change format of reaction equations
reaction_formulas = strrep(reaction_formulas, '@MNXD1', '[0]');
reaction_formulas = strrep(reaction_formulas, '@MNXD2', '[1]');
reaction_formulas = strrep(reaction_formulas, '@MNXD3', '[0]');
reaction_formulas = strrep(reaction_formulas, '@MNXD4', '[0]');
reaction_formulas = strrep(reaction_formulas, '@MNXD5', '[0]');
reaction_formulas = strrep(reaction_formulas, '@MNXDX', '[0]');
reaction_formulas = strrep(reaction_formulas, 'MNXM01', 'MNXM1');
reaction_formulas = strrep(reaction_formulas, '=', '<=>');

% remove all reactions that involve BIOMASS or EMPTY
keep_idx = ~contains(reaction_formulas, 'BIOMASS') & ...
    ~contains(reaction_ids, 'EMPTY');
reaction_ids = reaction_ids(keep_idx);
reaction_formulas = reaction_formulas(keep_idx);

% add colon to reaction IDs
reaction_ids = strcat(reaction_ids, ':');

% write file
writetable(...
    cell2table([reaction_ids reaction_formulas]),...
    fullfile(mnxref_dir, 'reaction_MNXref_balanced.lst'), ...
    'Delimiter', '\t', ...
    'FileType', 'text', ...
    'WriteVariableNames', false ...
    );
    
end