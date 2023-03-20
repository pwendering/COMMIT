function changed_model = removeDuplicateRxns(model, database, verbose)
%% changed_model = removeDuplicateRxns(model, database, verbose)
% Check for identity of reactions that share the same identifier and remove
% the identical reactions from the model.
% Input:
%           struct model:           metabolic model
%           struct database:        database model to be checked against
%                                   (ensure mass balanced reactions!)
%           logical verbose:        true if function should print progess
%                                   and warnings on stdout
% Output:
%           struct changed_model:   updated model

warning('off','all')

%% Check user input
if nargin == 0
    error('USAGE: removeDuplicateRxns(model, database, verbose(optional)')
elseif nargin == 1
    error('You did not give enough arguments.')
elseif nargin == 2
    verbose = 0;
end

fid1 = fopen(fullfile('data','models','consensus','duplicate_reactions.txt'), 'w');
fprintf(fid1, 'ID\tEquation\n');

%% Initialize class
if verbose; fprintf('\n\t> Initializing the MergeModelsClass object... '); end

tic;
M = mergeModelsClass(model, database);
if verbose; fprintf('done. (%.2fs)\n', toc); end
n_rxns = M.getRxnNum(model);

% find biomass reaction index
% * by objective vector
% * by matching keyword 'biomass'
BIOMASS = union(find(M.model.c),...
    find(contains(M.model.rxns, 'biomass', 'IgnoreCase', true)));
fprintf('Detected the following biomass (-related) reaction(s): %s\n',...
    strjoin(M.model.rxns(BIOMASS), ', '))

% remove eventually pre-existing number inidices from reactions
M.model = M.removeRxnIndices;

if verbose
    fprintf('\n\t> Searching for duplicate reactions...\n')
end

%% Start search
for  i=1:n_rxns

    % create reaction 1
    R1 = Reaction(...
        M.getProperty(i, 'rxns'),...
        i,...
        M.getRxnMets(i),...
        M.getProperty(i, 'lb'),...
        M.getProperty(i, 'ub'),...
        M.checkMassBalance(i)...
        );
    
    if ~ismember(i, M.remove)
        
        %% completely identical reactions
        % row of the similarity matrix for the current reaction
        identical = M.findIdentical(i);
        % exclude reactions that should be removed
        identical = M.filterArray(identical, M.remove);
        % exclude biomass reaction (ATPase / NGAM reaction matches with
        % high identity)
        identical = M.filterArray(identical, BIOMASS);
        % exclude all preceding reactions
        identical = M.filterArray(identical, [1:i]);
        identical = find(identical);
        identical = reshape(identical, 1, numel(identical));
        M = M.increaseCounter('count_identity', numel(identical));
        
        %% reactions with a lower identity score
        % lower identity threshold
        identical_lower = M.findIdenticalLower(i);
        % exclude reactions that should be removed
        identical_lower = M.filterArray(identical_lower, M.remove);
        % exclude exchange reactions
        identical_lower = M.filterArray(identical_lower, M.findExcReactions);
        % exclude biomass reaction (ATPase / NGAM reaction matches with
        % high identity)
        identical_lower = M.filterArray(identical_lower, BIOMASS);
        % exclude preceding reactions
        identical_lower = M.filterArray(identical_lower, 1:i);
        identical_lower = find(identical_lower);
        identical_lower = reshape(identical_lower, 1, numel(identical_lower));
        M = M.increaseCounter('count_identity_lower', numel(identical_lower));
        
        %% reactions with equal IDs
        equal_ids = M.findRxnIndex(R1.id);
        if ~isequal(equal_ids, R1.index)
            equal_ids = M.filterArray(equal_ids, M.remove);
            equal_ids = M.filterArray(equal_ids, [1:i]);
            equal_ids = M.filterArray(equal_ids, identical);
            equal_ids = M.filterArray(equal_ids, identical_lower);
            equal_ids = reshape(equal_ids, 1, numel(equal_ids));
            M = M.increaseCounter('count_equal_ids', numel(equal_ids));
        else
            equal_ids = [];
        end
        
        
        candidates = unique([identical, identical_lower, equal_ids]);
        if isempty(candidates) || isequal(candidates, 0)
            continue
        else
            % convert to column vector (for iteration)
            candidates = reshape(candidates, 1, numel(candidates));
        end
        
        % print ID and equation to file
        M = M.changeRxnID(R1.index, strcat(R1.id, '_R1'));
        fprintf(fid1, '%s\t%s\n', char(strcat(R1.id, '_R1')), char(printRxnFormula(M.model, strcat(R1.id, '_R1'), 0)));
        M = M.changeRxnID(R1.index, R1.id);
        
        %% Start comparison
        for j = candidates
            
            % create reaction 2
            R2 = Reaction(...
                M.getProperty(j, 'rxns'),...
                j,...
                M.getRxnMets(j),...
                M.getProperty(j, 'lb'),...
                M.getProperty(j, 'ub'),...
                M.checkMassBalance(j)...
                );
            
            % print ID and equation to file
            M = M.changeRxnID(R2.index, strcat(R2.id, '_R2'));
            fprintf(fid1, '%s\t%s\n', char(strcat(R2.id, '_R2')), char(printRxnFormula(M.model, strcat(R2.id, '_R2'), 0)));
            M = M.changeRxnID(R2.index, R2.id);
            
            RP = ReactionPair(R1, R2, M.getSimilarity(i,j));
            
            if ~RP.compareIDs && ~RP.compareMets && numel(unique(RP.findMNXRxn))==2
                % IDs are not identical but both in MNXref namespace
                % AND they do not have the same metabolite set
                % (very likely different reactions that just happen
                % to be very similar) ==> do not consider this reaction
                continue
            end
            
            % take the MNXref id in case the other reaction is not in
            % MNXref namespace; if both IDs are not in MNXref namespace,
            % take the first ID
            rxn_id = RP.selectMNXref;
            if isempty(rxn_id)
                rxn_id = R1.id;
            end
            M = M.changeRxnID(i, rxn_id);
            M = M.changeRxnID(j, rxn_id);
            
            % initialize the variable that should indicate of reactions
            % that have one different metabolite have opposite directions
            opposite_dir = 0;
            
            if ~RP.identity_lower && numel(RP.getMetDiff) <= 2
                if RP.compareMets
                    % all metabolites are shared
                    RP = RP.setSimilarity(M.t_identity_lower);
                elseif numel(RP.getMetDiff) == 2
                    % one metabolite is different but not H+ AND the
                    % reactions have equal identifiers if in MNXref namespace
                    % (one per reaction)
                    blackList = RP.getMetDiff;
                    [tmp_identity, opposite_dir] = M.compareRxns(RP, blackList);
                    RP = RP.setSimilarity(tmp_identity);
                elseif contains(RP.getMetDiff, 'MNXM1[')
                    % one metabolite is different in this is H+
                    RP = RP.setSimilarity(M.t_identity_lower);
                end
                
            end
            
            if RP.identity_lower
                
                if RP.decideByMassBalance
                    decision = RP.decideByMassBalance;
                    M = M.appendRemove(decision);
                elseif RP.decideByProtons
                    decision = RP.decideByProtons;
                    M = M.appendRemove(decision);
                end
                
            end
            
            if RP.identity_lower && ~any(ismember([i j], M.remove))
                
                if ~R1.massBalance && ~R2.massBalance && (M.checkRxnDB(R1.id) || M.checkRxnDB(R2.id))
                    % none of them is mass-balanced but at least one is
                    % contained in the database
                    
                    if isequal(R1.id, RP.selectMNXref)
                        % IDs are equal or both are in MNXref
                        % namespace
                        new_name = M.getProperty(R1.index, 'rxnNames');
                        M = M.appendRemove(R2.index);
                        M = M.changeRxnID(R1.index, strcat(R1.id, '_delete'));
                    else
                        new_name = M.getProperty(R2.index, 'rxnNames');
                        M = M.appendRemove(R1.index);
                        M = M.changeRxnID(R2.index, strcat(R2.id, '_delete'));
                    end
                    
                    new_id = char(RP.selectMNXref);
                    new_equation = char(printRxnFormula(database, new_id, 0));
                    new_name = char(new_name);
                    
                    % Check whether there is another reaction with the same ID
                    tmp_idx = setdiff(M.findRxnIndex(new_id), [i j]);
                    if ~isempty(tmp_idx) && M.checkMassBalance(tmp_idx(1))
                        M = M.appendRemove(R1.index);
                        M = M.appendRemove(R2.index);
                        if verbose
                            fprintf('\t> Both removed, an identical reaction is already contained in the model: %s.\n', new_id)
                        end
                    else
                        M.model = addReaction(M.model, new_id,...
                            'reactionName', new_name,...
                            'reactionFormula', new_equation,...
                            'checkDuplicate', false,...
                            'printLevel', 0);
                        
                        % print ID and equation to file
                        fprintf(fid1, '%s\t%s\n', strcat(new_id, '_db'), new_equation);
                        
                        
                        % this field is added by the 'addReaction' function
                        if isfield(M.model, 'grRules')
                            M.model = rmfield(M.model, 'grRules');
                        end
                        
                        if verbose
                            fprintf('\t> Added reaction %s from the database\n', new_id)
                        end
                        
                        M = M.increaseCounter('count_added');
                    end
                    
                    clear new_name new_id new_equation tmp_idx
                end
            end
            
            % check if any decision has been made
            if ~any(ismember([i j], M.remove)) && (RP.identity_lower || RP.identity)
                % either identical or no decision has been made yet
                
                if RP.decideByReversibility
                    decision = RP.decideByReversibility;
                    M = M.appendRemove(decision);
                    if RP.compareReversibility
                        M = M.increaseCounter('count_diff_rev');
                    end
                elseif ~RP.compareDirections || opposite_dir
                    % reactions have opposite directions
                    Rev = RP.R2.createRevID;
                    if M.checkExistRxn(Rev.id)
                        % if the reaction already exists, remove the second
                        % reaction
                        M = M.appendRemove(R2.index);
                    else
                        M = M.changeRxnID(R2.index, Rev.id);
                        M = M.increaseCounter('count_opposite');
                    end
                else
                    % both unidirectional and same direction
                    M = M.appendRemove(R2.index);
                end
            end
            
            % Compare and merge GPR rules
            score = M.compareRxnGPR(i, j);
            if RP.identity
                M = M.appendGPRScoreIdentical(score);
            else
                M = M.appendGPRScoreLower(score);
            end
            
            try
                new_rule = M.mergeRxnGPR(i, j);
                M = M.changeRxnGPR(i, new_rule);
                M = M.changeRxnGPR(j, new_rule);
            catch EM
                if verbose
                    fprintf('%s\n', EM.message)
                end
            end
            
            fprintf(fid1, 'GPR Jaccard index:\t%.2f\n', score);
            fprintf(fid1, 'Cosine similarity:\t%.2f\n', M.getSimilarity(i,j));
            if ismember(j, M.remove)
                fprintf(fid1, 'Chose R1\n');
            elseif ismember(i, M.remove)
                fprintf(fid1, 'Chose R2\n');
            else
                fprintf(fid1, 'No decision\n');
            end
                
            if ismember(i, M.remove)
                % break the for loop here because later this reaction
                % will be compared to the succeeding duplicates
                break
            end
            fprintf(fid1, '\n');
        end
        fprintf(fid1, '\n');
    end
end

% change the idenitfier of the duplicate reactions to they can be
% discriminated from the reactions with identical IDs (The for loop allows
% to delete multiple reactions with the same ID, otherwise only the first
% of them would be deleted) model = M.model;
remove = M.remove;
model.rxns(remove) = strcat(model.rxns(remove), '_delete');
remove_ids = model.rxns(contains(model.rxns, '_delete'));

% 'EC' is not dealt with in removeRxns function 
if isfield(model, 'EC')
    model.EC(remove) = [];
end

for i=1:numel(remove_ids)
    model = removeRxns(model, remove_ids(i));
%     model.rules(remove(i)) = [];
end

if verbose
    fprintf('\n\t\t> %d completely identical reactions found.', M.count_identity)
    fprintf('\n\t\t> %d almost identical reactions found.', M.count_identity_lower)
    fprintf('\n\t\t> %d (almost) identical reactions had different reversibility', M.count_diff_rev)
    fprintf('\n\t\t> %d (almost) identical reactions had opposite directions', M.count_opposite)
    fprintf('\n\t\t> Average for GPR rule Jaccard index for identical reactions: %.2f', mean(M.gpr_scores_identical))
    fprintf('\n\t\t> Average for GPR rule Jaccard index for almost identical reactions: %.2f', mean(M.gpr_scores_identical_lower))
    fprintf('\n\t\t> %d reactions were replaced by reactions from the database', M.count_added)
    fprintf('\n\t==> Removed %d reactions.\n', n_rxns-numel(model.rxns))
end

model.S = sparse(model.S);
changed_model = model;
end
