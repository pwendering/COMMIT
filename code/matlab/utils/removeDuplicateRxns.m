function changed_model = removeDuplicateRxns(model, dbCheck, database, verbose)
% Check for identity of reactions that share the same identifier and remove
% the identical reactions from the model
% Input:
%           struct model:           metabolic model
%           logical dbCheck:        whether each reaction should be checked
%                                   against a reaction with the identical ID
%                                   in the provided database model
%           struct database:        database model to be checked against
%                                   (ensure mass balanced reactions!)
%           logical verbose:        true if function should print progess
%                                   and warnings on stdout
% Output:
%           struct changed_model:   updated model
fid1 = fopen('/stud/wendering/gene_rules_identical_reactions.txt', 'w');
fid2 = fopen('/stud/wendering/equations_lower_similarity.txt', 'w');
fprintf(fid2, 'ID\tEquation\tgprJaccard\n');
fprintf(fid1, 'IDs\tRules\n');
if nargin == 1
    dbCheck = false;
elseif ~islogical(dbCheck)
    error('Argument dbCheck is not of logical type')
elseif nargin == 2
    error('You did not give a database to search against')
elseif nargin == 3
    verbose = 0;
end

% check mass balancing for the model
[~, imbalanced] = checkMassChargeBalance(model);
imbalanced = ~cellfun('isempty', imbalanced);

%% Find reactions with equal stoichiometries
if verbose
    fprintf('\nFind duplicate reactions by stoichiometry\n')
    fprintf('\n\t> Calculating the Cosine similarity matrix... ')
end
n_rxns = numel(model.rxns);

stoich_mat = full(model.S);
similarity_matrix = zeros(n_rxns);

tic;
parfor i=1:n_rxns-1
    A = stoich_mat(:,i);
    row_tmp = zeros(1, n_rxns);
    for j=i+1:n_rxns
        B = stoich_mat(:,j);
        row_tmp(j) = cosineSimilarity(A, B);
    end
    similarity_matrix(i,:) = row_tmp;
end
if verbose
    fprintf('done. (%.2fs)\n', toc)
end

% identity thresholds
epsilon = 10E-6;
t_identity = 1-epsilon;
t_identity_lower = 0.95;

% counters
count_identity_lower = 0;
count_identity = 0;
count_diff_rev = 0;
count_opposite = 0;
count_added = 0;

% array that contains indices of reactions that should be removed
remove = [];

% declare the vectors to stors gpr Jaccard distances for identical and
% almost identicalreactions
gpr_scores_id = [];
gpr_scores_idl = [];

% find exchange reactions
sink_rxns = (sum(model.S==-1,1)==1) & (sum(model.S~=0) == 1);
uptake_rxns = (sum(model.S==1,1)==1) & (sum(model.S~=0) == 1);
exchange_rxns = sink_rxns + uptake_rxns;
clear sink_rxns uptake_rxns

% remove eventually pre-existing number inidices from reactions
model.rxns = regexprep(model.rxns, '_[0-9]*$', '');

if verbose
    fprintf('\n\t> Searching for identical reactions...\n')
end
for  i=1:n_rxns
    
    % variables for the current reaction
    current_id = model.rxns{i};
    
    if ~ismember(i, remove)
        
        % row of the similarity matrix for the current reaction
        current_row = similarity_matrix(i,:);
        
        % identical reactions
        identical = (current_row >= t_identity) | (current_row <= -t_identity);
        % exclude reactions that should be removed
        identical(remove) = 0;
        % exclude all preceding reactions
        identical = setdiff(find(identical), 1:i);
        % lower identity threshold
        identical_lower = (current_row >= t_identity_lower) | (current_row <= -t_identity_lower);
        % exclude reactions that should be removed
        identical_lower(remove) = 0;
        % exclude identical, exchange and preceding reactions
        identical_lower(identical) = 0;
        identical_lower(logical(exchange_rxns)) = 0;
        identical_lower = setdiff(find(identical_lower), 1:i);
        
        % reversibility of the current reaction
        rev_i = model.lb(i) & model.ub(i);
        
        %% Completely identical reactions
        if ~isempty(identical)
            count_identity = count_identity + numel(identical);
            fprintf(fid1, '%s\t%s\n', char(current_id), model.rules{i});
            for j=1:numel(identical)
                
                % index of the current duplicate candidate
                idx_cmp = identical(j);
                id_cmp = model.rxns{idx_cmp};
                
                % determine reversibility of the current duplicate candidate
                rev_j = model.lb(idx_cmp) & model.ub(idx_cmp);
                
                % take the MNXref id in case the other reaction in not in
                % MNXref namespace
                rxn_id = unique(regexp([current_id, id_cmp], 'MNXR\d+_\w', 'match'));
                if numel(rxn_id) == 2
                    % different IDs in MNXref
                    rxn_id = current_id;
                else
                    rxn_id = char(rxn_id);
                end
                
                if rev_i && rev_j
                    % both reversible ==> remove the second reaction
                    remove = vertcat(remove, idx_cmp);
                    % rename the first ID
                    model.rxns{i} = rxn_id;
                elseif rev_i
                    % only the current reaction is reversible
                    count_diff_rev = count_diff_rev + 1;
                    remove = vertcat(remove, i);
                    model.rxns{idx_cmp} = rxn_id;
                elseif rev_j
                    % only the other reaction is reversible
                    count_diff_rev = count_diff_rev + 1;
                    remove = vertcat(remove, idx_cmp);
                    model.rxns{i} = rxn_id;
                elseif current_row(idx_cmp)<0
                    % both unidirectional and opposite directions
                    id_cmp_rev = strcat(rxn_id, '_r');
                    if ~ismember(id_cmp_rev, model.rxns)
                        % if the reversible reaction is not contained in
                        % the model yet, it will be added
                        model.rxns{idx_cmp} = id_cmp_rev;
                        count_opposite = count_opposite + 1;
                    else
                        % if already contained, remove the current
                        % duplicate candidate
                        remove = vertcat(remove, idx_cmp);
                        model.rxns{i} = rxn_id;
                    end
                elseif current_row(idx_cmp)>0
                    % both unidirectional and same direction
                    remove = vertcat(remove, idx_cmp);
                    model.rxns{i} = rxn_id;
                else
                    if verbose
                        warning('No reaction found to delete: %s\n', char(current_id))
                    end
                    continue
                end
                
                
                % Compare and merge GPR rules
                gpr_scores_id = vertcat(gpr_scores_id, compareGPR(model, i, model, idx_cmp));
                try
                    merged_rule = mergeGPR(model.rules(i), model.rules(idx_cmp));
                catch EM
                    if verbose
                        fprintf('%s\n', EM.message)
                    end
                end
                if remove(end)==i
                    % if reaction i will be removed, merge GPR rules in the
                    % other reaction
                    model.rules(idx_cmp) = merged_rule;
                    % break the for loop here because later this reaction
                    % will be compared to the succeeding duplicates
                    fprintf(fid1, '%s\t%s\n', model.rxns{idx_cmp}, model.rules{identical(j)});
                    fprintf(fid1, 'Merged\t%s\n', merged_rule{:});
                    break
                else
                    % otherwise if the other reaction will be removed
                    model.rules(i) = merged_rule;
                end
                
                fprintf(fid1, '%s\t%s\n', model.rxns{idx_cmp}, model.rules{identical(j)});
                fprintf(fid1, 'Merged\t%s\n', merged_rule{:});
            end
            clear id_cmp idx_cmp rxn_id
            fprintf(fid1, '\n');
            
        end
        
        %% Reaction that have lower similarity but are almost identical
        if ~isempty(identical_lower)
            count_identity_lower = count_identity_lower + numel(identical_lower);
            
            mass_balance_current = 0;
            if dbCheck && ismember(current_id, database.rxns)
                % Is the database reaction with the identical ID identical
                % to the current reaction?
                mass_balance_current = checkIdenticalReaction(model, i,...
                    database, current_id, {}, true);
            end
            
            if ~mass_balance_current
                % if the reaction is not identical to the database
                % reaction, check mass balance by the formula
                mass_balance_current = ~imbalanced(i);
            end
            
            id_orig = model.rxns(i);
            model.rxns(i) = strcat(model.rxns(i), '_R1');
            fprintf(fid2, '%s\t%s\n', model.rxns{i},...
                char(printRxnFormula(model, model.rxns(i), 0)));
            model.rxns(i) = id_orig;
            
            
            for j=1:numel(identical_lower)
                
                % index and ID of the current duplicate candidate
                idx_cmp = identical_lower(j);
                id_cmp = model.rxns{idx_cmp};
                
                % determine reversibility of the current duplicate candidate
                rev_j = model.lb(idx_cmp) & model.ub(idx_cmp);
                
                mass_balance_cmp = 0;
                if dbCheck && ismember(id_cmp, database.rxns)
                    % Is the database reaction with the identical ID identical
                    % to the current reaction?
                    mass_balance_cmp = checkIdenticalReaction(model, idx_cmp,...
                        database, current_id, {}, true);
                end
                
                if ~mass_balance_cmp
                    % is the reaction mass-balanced?
                    mass_balance_cmp = ~imbalanced(idx_cmp);
                end
                
                model.rxns(idx_cmp) = strcat(model.rxns(idx_cmp), '_R2');
                fprintf(fid2, '%s\t%s\n', model.rxns{idx_cmp},...
                    char(printRxnFormula(model, model.rxns(idx_cmp), 0)));
                model.rxns{idx_cmp} = id_cmp;
                
                % metabolite its for both reactions
                mets_i = model.mets(any(model.S(:, i), 2));
                mets_j = model.mets(any(model.S(:, idx_cmp), 2));
                
                ids_equal = isequal(current_id, id_cmp);
                ids_mnxref = all(contains([current_id, id_cmp], 'MNXR'));
                mets_sym_diff = setxor(mets_i, mets_j);
                
                if ~ids_equal && ids_mnxref && ~isempty(mets_sym_diff)
                    % IDs are not identical but both in MNXref namespace
                    % AND they do not have the same metabolite set
                    % (very likely different reactions that just happen
                    % to be very similar) ==> do not consider this reaction
                    continue
                end
                
                % reactions have the same set of metabolites
                if mass_balance_current && ~mass_balance_cmp
                    % The current is mass-balanced, while the other one
                    % is not
                    remove = vertcat(remove, idx_cmp);
                elseif mass_balance_cmp && ~mass_balance_current
                    % The current reaction is not mass-balanced, while the
                    % second one is
                    remove = vertcat(remove, i);
                elseif ~mass_balance_current && ~mass_balance_cmp && any(ismember([current_id, id_cmp], database.rxns))
                    % Both reactions are not mass-balanced and
                    % at least one of the IDs is contained in the database
                    if ids_equal || (~ids_mnxref && contains(current_id, 'MNXR')) || ids_mnxref
                        % IDs are equal or both are in MNXref
                        % namespace
                        new_id = current_id;
                        new_name = model.rxnNames{i};
                        model.rxns{i} = strcat(current_id, '_delete');
                        remove = vertcat(remove, idx_cmp);
                    else
                        % only one is in MNXref namespace
                        new_id = id_cmp;
                        new_name = model.rxnNames{idx_cmp};
                        model.rxns{idx_cmp} = strcat(id_cmp, '_delete');
                        remove = vertcat(remove, i);
                    end
                    
                    new_equation = char(printRxnFormula(database, new_id, 0));
                    
                    
                    % remove the second reaction from the model first and
                    % change the ID so it will not be updated
                    
                    if ismember(new_id, model.rxns)
                        % there is another reaction with the same ID
                        % ==> check mass balance
                        prev_identical = checkIdenticalReaction(model, new_id,...
                            database, new_id, {}, true);
                    else
                        prev_identical = 0;
                    end
                    
                    % add the reaction from the database
                    if ~prev_identical
                        n_tmp = numel(model.rxns);
                        model = addReaction(model, new_id,...
                            'reactionName', new_name,...
                            'reactionFormula', new_equation,...
                            'checkDuplicate', false,...
                            'printLevel', 0);
                        
                        if n_tmp > numel(model.rxns)
                            % if a new reaction is added (not only updated),
                            % remove the current one
                            if remove(end)==i
                                remove = vertcat(remove, idx_cmp);
                            else
                                remove = vertcat(remove, i);
                            end
                        end
                        
                        clear n_tmp new_name new_id new_equation
                        % this field is added by the 'addReaction' function
                        if isfield(model, 'grRules')
                            model = rmfield(model, 'grRules');
                        end
                        
                        fprintf(fid2, '%s\t%s\n', strcat(current_id, '_db'),...
                            char(printRxnFormula(model, current_id, 0)));
                        
                        if verbose
                            fprintf('\t> Added reaction %s from the database\n', char(current_id))
                        end
                        count_added = count_added + 1;
                    elseif verbose
                        fprintf('\t> Both removed, an identical is already contained in the model: %s.\n', current_id)
                    end
                else
                    % either none of them are mass-balanced and not
                    % contained in the database or both are mass-balanced
                    % ==> decide by presence of protons and reversibility
                    
                    % get metabolites that are only contained in the compared
                    mets_diff_i = setdiff(mets_i, mets_j);
                    mets_diff_j = setdiff(mets_j, mets_i);
                    
                    if numel(mets_sym_diff)==1 && ~isempty(regexp(mets_diff_i, 'MNXM1[.]'))
                        % the first reaction contains H+ whereas the second one does
                        % not ==> choose the first reaction
                        remove = vertcat(remove, idx_cmp);
                    elseif numel(mets_sym_diff)==1 && ~isempty(regexp(mets_diff_j, 'MNXM1[.]'))
                        % the second reaction contains H+ whereas the first one does
                        % not ==> choose the first reaction
                        remove = vertcat(remove, i);
                    elseif rev_i && rev_j
                        % almost identical both reversible
                        remove = vertcat(remove, idx_cmp);
                    elseif rev_i
                        % almost identical and only the current reaction is reversible
                        count_diff_rev = count_diff_rev + 1;
                        remove = vertcat(remove, i);
                        % break the for loop here because later this reaction
                        % will be compared to the succeeding duplicates
                        break
                    elseif rev_j
                        % almost identical and only the other reaction is reversible
                        count_diff_rev = count_diff_rev + 1;
                        remove = vertcat(remove, idx_cmp);
                    elseif current_row(idx_cmp)<0
                        % almost identical, both unidirectional with opposite
                        % directions
                        id_cmp_rev = strcat(id_cmp, '_r');
                        if ~ismember(id_cmp_rev, model.rxns)
                            % if the reversible reaction is not contained in
                            % the model yet, it will be added
                            model.rxns{idx_cmp} = id_cmp_rev;
                            count_opposite = count_opposite + 1;
                        else
                            % if already contained, remove the current
                            % duplicate candidate
                            remove = vertcat(remove, idx_cmp);
                        end
                    elseif current_row(idx_cmp)>0
                        % both unidirectional and same direction
                        remove = vertcat(remove, idx_cmp);
                    elseif verbose
                        warning('Two almost identical reactions but none could be chosen: %s', char(current_id))
                    end
                end
                
                % Compare and merge GPR rules
                gpr_scores_idl = vertcat(gpr_scores_idl, compareGPR(model, i, model, idx_cmp));
                try
                    merged_rule = mergeGPR(model.rules(i), model.rules(idx_cmp));
                catch EM
                    if verbose
                        fprintf('%s\n', EM.message)
                    end
                end
                fprintf(fid2, '\t\t%.2f\n', gpr_scores_idl(end));
                
                if remove(end)==i
                    % if reaction i will be removed, merge GPR rules in the
                    % other reaction
                    model.rules(idx_cmp) = merged_rule;
                    fprintf(fid2, 'Chose R2\n');
                    % break the for loop here because later this reaction
                    % will be compared to the succeeding duplicates
                    break
                else
                    % otherwise if the other reaction will be removed
                    model.rules(i) = merged_rule;
                    fprintf(fid2, 'Chose R1\n');
                end
            end
            fprintf(fid2, '\n');
        end
    end
end

% change the idenitfier of the duplicate reactions to they can be
% discriminated from the reactions with identical IDs (The for loop allows
% to delete multiple reactions with the same ID, otherwise only the first
% of them would be deleted)
model.rxns(remove) = strcat(model.rxns(remove), '_delete');
remove_ids = model.rxns(contains(model.rxns, '_delete'));
for i=1:numel(remove_ids)
    model = removeRxns(model, remove_ids(i));
end

if verbose
    fprintf('\n\t\t> %d completely identical reactions found.', count_identity)
    fprintf('\n\t\t> %d almost identical reactions found.', count_identity_lower)
    fprintf('\n\t\t> %d (almost) identical reactions had different reversibility', count_diff_rev)
    fprintf('\n\t\t> %d (almost) identical reactions had opposite directions', count_opposite)
    fprintf('\n\t\t> Average for GPR rule Jaccard index for identical reactions: %.2f', mean(gpr_scores_id))
    fprintf('\n\t\t> Average for GPR rule Jaccard index for almost identical reactions: %.2f', mean(gpr_scores_idl))
    fprintf('\n\t\t> %d reactions were replaced by reactions from the database', count_added)
    fprintf('\n\t==> Removed %d reactions.\n', n_rxns-numel(model.rxns))
end
fclose(fid1);


%% Find reactions with equal identifiers
fprintf(fid2, '---------------------------------------------------------------------------\n');
fprintf(fid2, 'From identical IDs\n\n');
if verbose
    fprintf('\nFind duplicate reactions by ID\n')
end

% remove indices from reaction IDs if present
model.rxns = regexprep(model.rxns, '_[0-9]*$', '');

% Check mass balancing
[~, imbalanced] = checkMassChargeBalance(model);
imbalanced = ~cellfun('isempty', imbalanced);

n_rxns = numel(model.rxns);
remove = [];
count_imbalanced = 0;

for i=1:n_rxns
    
    % current ID
    current_id = model.rxns(i);
    % indices of the repeated occurrences of the same reaction
    idx_cmp = find(strcmp(model.rxns, current_id));
    idx_cmp = setdiff(idx_cmp, union(remove, i));
    
    % index to add to non-identical mass-balanced reaction with the same ID
    add_index = 2;
    
    % reversibility of the current reaction
    rev_i = model.lb(i) & model.ub(i);
    
    % check the identity of the current reaction with the (mass-balanced)
    % database reaction
    mass_balance_current = 0;
    if dbCheck && ismember(current_id, database.rxns)
        mass_balance_current = checkIdenticalReaction(model, i, database, current_id,...
            {}, true);
    end
    
    if ~mass_balance_current
        mass_balance_current = ~imbalanced(i);
    end
    
    % go through the duplicates and check for stoichiometric identity
    for j = idx_cmp
        
        % metabolite its for both reactions
        mets_i = model.mets(any(model.S(:, i), 2));
        mets_j = model.mets(any(model.S(:, j), 2));
        % get metabolites that are only contained in the compared
        mets_diff_j = setdiff(mets_j, mets_i);
        mets_diff_i = setdiff(mets_i, mets_j);
        mets_sym_diff = setxor(mets_i, mets_j);
        
        % reversibility of the compared reaction
        rev_j = model.lb(j) & model.ub(j);
        
        if isempty(mets_sym_diff)
            % all metabolite IDs equal ==> check for idenitity
            [mass_balance_cmp, opposite] = checkIdenticalReaction(model, i,...
                model, j, {}, true);
        elseif ~isempty(regexp(mets_diff_i, 'MNXM1[.]'))
            % the first reaction contains H+ whereas the second one does
            % not ==> choose the first reaction
            remove = vertcat(remove, j);
            continue
        elseif ~isempty(regexp(mets_diff_j, 'MNXM1[.]'))
            % the second reaction contains H+ whereas the first one does
            % not ==> choose the first reaction
            remove = vertcat(remove, i);
            continue
        elseif (numel(mets_diff_i) == 1) && (numel(mets_diff_j) == 1)
            % one metabolite is different between the reactions
            % exclude the symmetric difference from the identity check
            blackList = vertcat(mets_diff_i, mets_diff_j);
            [mass_balance_cmp, opposite] = checkIdenticalReaction(model, i,...
                model, j, blackList, true);
        else
            % more than one metabolite difference
            mass_balance_cmp = 0;
        end
        
        if ~mass_balance_cmp
            % if two reactions are not identical but have the same identifier,
            % chose the one that is identical with the database
            if dbCheck && ismember(current_id, database.rxns)
                % check identity of the compared reaction with the
                % mass-balanced database
                mass_balance_cmp = checkIdenticalReaction(model, j,...
                    database, current_id, {}, true);
            end
            if ~mass_balance_cmp
                mass_balance_cmp = ~imbalanced(j);
            end
            
            if mass_balance_current && mass_balance_cmp
                % both reactions are mass-balanced but not identical
                % assume that both are correct and add an index to the
                % compared reaction
                model.rxns(j) = strcat(model.rxns(j),...
                    '_', num2str(add_index));
                add_index = add_index + 1;
            elseif mass_balance_current
                % if only the current reaction is mass-balanced,
                % remove the compared reaction
                remove = vertcat(remove, j);
            elseif mass_balance_cmp
                % if only the compared reaction is mass-balanced,
                % remove the current reaction
                remove = vertcat(remove, i);
            elseif ismember(current_id, database.rxns)
                % none of them is mass-balanced, so poor quality is
                % assumed and the reaction from the database is added
                remove = vertcat(remove, j);
                n_tmp = numel(model.rxns);
                
                id_cmp_orig = model.rxns(j);
                model.rxns(j) = strcat(id_cmp_orig,...
                    '_R2');
                fprintf(fid2, '%s\t%s\n', strcat(model.rxns{i}, '_R1'),...
                    char(printRxnFormula(model, model.rxns(i), 0)));
                fprintf(fid2, '%s\t%s\n', model.rxns{j},...
                    char(printRxnFormula(model, model.rxns(j), 0)));
                
                model.rxns(j) = id_cmp_orig;
                
                model.rxns(i) = strcat(model.rxns(i), '_delete');
                
                
                % set the direction based on the consensus of the two
                % reactions
                new_lb = -1000;
                new_ub = 1000;
                
                %                 % cosine similarity to determine whether the database
                %                 % reaction has the same direction (coefficients)
                %                 t_dir = 0.9;
                %                 db_rxn = find(strcmp(database.rxns, current_id));
                %                 shared_mets = intersect(database.mets(any(database.S(:, db_rxn), 2)),...
                %                     union(mets_i, mets_j));
                %                 coeff_db = database.S(cell2mat(cellfun(@(x)find(contains(database.mets, x)),...
                %                     shared_mets,'UniformOutput', false)), db_rxn);
                %                 coeff_i = model.S(cell2mat(cellfun(@(x)find(contains(model.mets, x)),...
                %                     shared_mets,'UniformOutput', false)),i);
                %                 coeff_j = model.S(cell2mat(cellfun(@(x)find(contains(model.mets, x)),...
                %                     shared_mets,'UniformOutput', false)),j);
                %                 similarity_i = cosineSimilarity(coeff_db, coeff_i);
                %                 similarity_j = cosineSimilarity(coeff_db, coeff_j);
                %
                %                 if (similarity_i <= -t_dir) && (similarity_j <= -t_dir)
                %                     % both stoichiometries opposite to database
                %                     new_lb = -1000 - max(model.lb(i), model.lb(j));
                %                     new_ub = 1000 - min(model.ub(i), model.ub(j));
                %                 elseif (similarity_i <= -t_dir) && (similarity_j >= t_dir)
                %                     % stoichiometry of first reaction opposite to db
                %                     new_lb = max(-1000 - model.lb(i), model.lb(j));
                %                     new_ub = min(1000 - model.ub(i), model.ub(j));
                %                 elseif (similarity_i >= t_dir) && (similarity_j <= -t_dir)
                %                     % stoichiometry of second reaction opposite to db
                %                     new_lb = max(model.lb(i), -1000 - model.lb(j));
                %                     new_ub = min(model.ub(i), 1000 - model.ub(j));
                %                 elseif (similarity_i >= t_dir) && (similarity_j >= t_dir)
                %                     new_lb = max(model.lb(i), model.lb(j));
                %                     new_ub = min(model.ub(i), model.ub(j));
                %                 else
                %                     warning('Reaction has a similarity of less than %.2f with the database reaction', t_dir)
                %                     new_lb = -1000; new_ub = 1000;
                %                 end
                %
                %                 if new_lb==0 && new_ub==0
                %                     new_lb = -1000; new_ub = 1000;
                %                 end
                
                
                model = addReaction(model, char(current_id),...
                    'reactionName', model.rxnNames{i},...
                    'reactionFormula', char(printRxnFormula(database, current_id, 0)),...
                    'lowerBound', new_lb,...
                    'upperBound', new_ub,...
                    'printLevel', 0);
                
                fprintf(fid2, '%s\t%s\n', strcat(char(current_id), '_db'),...
                    char(printRxnFormula(model, current_id, 0)));
                
                
                if n_tmp > numel(model.rxns)
                    % if a new reaction is added (not only updated),
                    % remove the current one
                    remove = vertcat(remove, i);
                end
                clear n_tmp new_formula new_lb new_ub db_idx similarity_i similarity_j
                % this field is added by the 'addReaction' function
                if isfield(model, 'grRules')
                    model = rmfield(model, 'grRules');
                end
                count_imbalanced = count_imbalanced + 1;
                if verbose
                    fprintf('\t==> Added reaction %s from the database\n', char(current_id))
                end
            else
                warning('Two imbalanced, non-identical reactions have the same ID %s (not is database).',...
                    char(current_id))
            end
        else
            % the two reactions are identical
            if rev_i && rev_j
                % identical and reversible
                remove = vertcat(remove, j);
            elseif rev_i
                % identical and only the current reaction is reversible
                remove = vertcat(remove, i);
            elseif rev_j
                % identical and only the other reaction is reversible
                remove = vertcat(remove, j);
            elseif opposite
                % identical, both unidirectional with opposite
                % directions
                id_cmp_rev = strcat(current_id, '_r');
                if ~ismember(id_cmp_rev, model.rxns)
                    % if the reversible reaction is not contained in
                    % the model yet, it will be added
                    model.rxns(j) = id_cmp_rev;
                else
                    % if already contained, remove the current
                    % duplicate candidate
                    remove = vertcat(remove, j);
                end
            else
                % both unidirectional and same direction
                remove = vertcat(remove, j);
            end
            
        end
        
        % Compare and merge GPR rules
        score = compareGPR(model, i, model, j);
        try
            merged_rule = mergeGPR(model.rules(i), model.rules(j));
        catch EM
            if verbose
                fprintf('%s\n', EM.message)
            end
        end
        
        if ~isempty(remove)
            if remove(end)==i
                % if reaction i will be removed, merge GPR rules in the
                % other reaction
                model.rules(j) = merged_rule;
                fprintf(fid2, 'Chose R2\n');
                % break the for loop here because later this reaction
                % will be compared to the succeeding duplicates
                break
            else
                % otherwise if the compared reaction will be removed
                model.rules(i) = merged_rule;
                fprintf(fid2, 'Chose R1\n');
            end
        end
        fprintf(fid2, '\t\t%.2f\n', score);
        fprintf(fid2, '\n');
    end
end
fclose(fid2);
% remove the marked reactions from the model
model.rxns(remove) = strcat(model.rxns(remove), '_delete');
remove_ids = model.rxns(contains(model.rxns, '_delete'));
for i=1:numel(remove_ids)
    model = removeRxns(model, remove_ids(i));
end

if verbose
    fprintf('\n\t==> Removed %d reactions.\n', n_rxns-numel(model.rxns))
    fprintf('\n\t> %d were non-identical and not mass-balanced\n', count_imbalanced)
    fprintf('time: %.2fs\n', toc)
end
model.S = sparse(model.S);
changed_model = model;
end
