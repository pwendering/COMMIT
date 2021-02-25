classdef mergeModelsClass
    properties
        model                   % merged stoichiometric model
        database                % database used for stoichiometry check
        similarity_matrix       % cosine similarity matrix
        mass_balance            % binary vector; 1 of mass-balanced
        remove                  % indices of reactions to be removed
        count_identity          % conter for identical reactions
        count_identity_lower    % coutner for 'almost' identical reactions
        count_equal_ids         % counter for reactions that have equal IDs
        count_diff_rev          % counter for different reversibility
        count_opposite          % counter for opposite directions
        count_added             % counter for reactions added from database
        gpr_scores_identical    % gpr Jaccard index for identical reactions
        gpr_scores_identical_lower % gpr Jaccard index for 'almost' identical reactions
    end
    
    properties (Constant)
        epsilon = 10E-6;
        t_identity = 1;
        t_identity_lower = 0.9;
    end
    
    methods
        function M = mergeModelsClass(model, database)
            % constructor method
            if nargin == 0
                M.model = struct;
                M.database = struct;
                M.similarity_matrix = 0;
                M.mass_balance = 0;
            elseif nargin == 1
                M.model = model;
                M.database = struct;
                M.similarity_matrix = M.getSimilarityMatrix;
                M.mass_balance = M.getMassBalance;
            else
                M.model = model;
                M.database = database;
                M.similarity_matrix = M.getSimilarityMatrix;
                M.mass_balance = M.getMassBalance;
            end
            M.remove = [];
            M.count_identity = 0;
            M.count_identity_lower = 0;
            M.count_equal_ids = 0;
            M.count_added = 0;
            M.count_diff_rev = 0;
            M.count_opposite = 0;
            M.gpr_scores_identical = [];
            M.gpr_scores_identical_lower = [];
        end
        function mb = getMassBalance(M, model)
            if nargin == 1
                model = M.model;
            end
            if M.checkField('metFormulas', model)
                [~, imbalanced] = checkMassChargeBalance(model);
                % 1 if mass-balanced, 0 otherwise
                mb = cellfun('isempty', imbalanced);
            else
                mb = 0;
            end
        end
        function n = getRxnNum(M, model)
            if nargin == 1
                model = M.model;
            end
            if M.checkField('rxns', model)
                n = numel(model.rxns);
            else
                n = 0;
            end
        end
        function similarity_matrix = getSimilarityMatrix(M, model)
            if nargin == 1
                model = M.model;
            end
            
            if M.checkField('S', model)
                % Cosine similarity matrix
                n = numel(model.rxns);
                stoich_mat = full(model.S);
                similarity_matrix = zeros(n);
                
                parfor i=1:n-1
                    A = stoich_mat(:,i);
                    row_tmp = zeros(1,n);
                    for j=i+1:n
                        B = stoich_mat(:,j);
                        row_tmp(j) = cosineSimilarity(A,B);
                    end
                    similarity_matrix(i,:) = row_tmp;
                end
            else
                similarity_matrix = 0;
            end
        end
        function M = setDatabase(M, database)
            M.database = database;
        end 
        function identical = findIdentical(M, index, matrix)
            % find completely identical reactions in the cosine similarity
            % matrix compared to the current one (index)
            if nargin == 2
                matrix = M.similarity_matrix;
            end
            row = M.getRow(matrix, index);
            e = mergeModelsClass.epsilon;
            t = mergeModelsClass.t_identity - e;
            identical = (row >= t) | ( row <= -t);
        end
        function identical_lower = findIdenticalLower(M, index, matrix)
            % find 'almost' identical reactions in the stoichiometric matrix
            % compared to the current one (index)
            if nargin == 2
                matrix = M.similarity_matrix;
            end
            row = M.getRow(matrix, index);
            t_identical = mergeModelsClass.t_identity;
            t_lower = mergeModelsClass.t_identity_lower;
            identical_lower = (row >= t_lower) & (row < t_identical) | ( row <= -t_lower) & (row > -t_identical);
        end
        function check = checkField(M, fieldname, model)
            if nargin == 2
                model = M.model;
            end
            if isfield(model, fieldname)
                check = 1;
            else
                warning('Model property has no field ''%s''.', fieldname)
                check = 0;
            end
        end
        function exc_rxns = findExcReactions(M, model)
            if nargin == 1
                model = M.model;
            end
            
            if M.checkField('S', model)
                stoich_mat = model.S;
            else
                stoich_mat = 0;
            end

            sink_rxns = (sum(stoich_mat==-1,1)==1) & (sum(stoich_mat~=0) == 1);
            uptake_rxns = (sum(stoich_mat==1,1)==1) & (sum(stoich_mat~=0) == 1);
            exc_rxns = sink_rxns + uptake_rxns;
        end
        function model = removeRxnIndices(M, model)
            if nargin == 1
                if M.checkField('rxns')
                    model = M.model;
                else
                    return
                end
            end
            model.rxns = regexprep(model.rxns, '_[0-9]*$', '');
        end
        function n = getFieldRange(M, property, model)
            if nargin == 2
                model = M.model;
            end
            
            if M.checkField(property, model)
                n = numel(model.(property));
            else
                n = 0;
            end
            
        end
        function p = getProperty(M, rxnIndex, property, model)
            if nargin == 3
                model = M.model;
            end
            
            if M.checkField(property, model)
                if rxnIndex <= M.getFieldRange(property, model)
                    p = model.(property)(rxnIndex);
                else
                    warning('Index %d is out of range for field ''%s''.',...
                    rxnIndex, property)
                end
            end
        end
        function s = getSimilarity(M, index_1, index_2, similarity_matrix)
            if nargin == 3
                similarity_matrix = M.similarity_matrix;
            end
            
            if ~isequal(similarity_matrix, 0)
                row = M.getRow(similarity_matrix, index_1);
                s = row(index_2);
            else
                s = 0;
            end
        end
        function mets = getRxnMets(M, rxnIndex, model)
            if nargin == 2
                model = M.model;
            end
            
            if M.checkField('mets', model) && M.checkField('S', model)
                if rxnIndex <= size(model.S, 2)
                    mets = model.mets(logical(model.S(:, rxnIndex)));
                else
                    warning('Index %d is out of range for field ''S''.',...
                    rxnIndex)
                    mets = {};
                end
            else
                mets = {};
            end
        end
        function mb = checkMassBalance(M, rxnIndex, model, database)
            % check the mass balance with the database or the COBRA
            % function 'checkMassChargeBalance'
            mb = 0;
            
            if nargin == 2
                model = M.model;
                database = M.database;
            elseif nargin == 3
                database = M.database;
            end
            
            % first, check with the mass balance array from the formulas
            if M.mass_balance ~= 0
                mb = M.mass_balance(rxn_index);
            end
            
            % if not mass-balanced, a reason for this could also be a
            % missing formula, so the identity with the according reaction
            % in the database is checked
            if ~mb && M.checkField('S', database)
                id = M.getProperty(rxnIndex, 'rxns');
                mb = checkIdenticalReaction(model, rxnIndex, database, id,...
                    {}, true);
            end
        end
        function M = increaseCounter(M, counter, value)
            if nargin == 2
                value = 1;
            end
            M.(counter) = M.(counter) + value;
        end
        function M = appendRemove(M, rxnIndex)
            M.remove = unique([M.remove, rxnIndex], 'stable');
        end
        function M = changeRxnID(M, rxnIndex, rxnID)
            if ischar(rxnID)
                rxnID = cellstr(rxnID);
            end
            
            if M.checkField('rxns')
                M.model.rxns(rxnIndex) = rxnID;
            end
        end
        function r = checkExistRxn(M, rxnID, model)
            if nargin == 2
                model = M.model;
            end
            
            if M.checkField('rxns', model)
                r = ismember(rxnID, model.rxns);
            end
        end
        function M = changeRxnGPR(M, rxnIndex, gpr)
            if ischar(gpr)
                gpr = cellstr(gpr);
            end
            if M.checkField('rules')
                M.model.rules(rxnIndex) = gpr;
            end
        end
        function score = compareRxnGPR(M, index_1, index_2)
            score = compareGPR(M.model, index_1, M.model, index_2);
        end
        function gpr = mergeRxnGPR(M, index_1, index_2)
            rule_1 = M.getProperty(index_1, 'rules');
            rule_2 = M.getProperty(index_2, 'rules');
            gpr = mergeGPR(rule_1, rule_2);
        end
        function M = appendGPRScoreIdentical(M, score)
            M.gpr_scores_identical = vertcat(M.gpr_scores_identical, score);
        end
        function M = appendGPRScoreLower(M, score)
            M.gpr_scores_identical_lower = vertcat(M.gpr_scores_identical_lower, score);
        end
        function checkDB = checkRxnDB(M, rxnID)
            checkDB = ismember(rxnID, M.database.rxns);
        end
        function index = findRxnIndex(M, rxnID)
            index = find(strcmp(M.model.rxns, rxnID));
        end
        function [similarity, opposite] = compareRxns(M, reactionPair, blackList)
            % column indices
            index_1 = reactionPair.R1.index;
            index_2 = reactionPair.R2.index;
            % metabolite indeices for the metabolites to be ignored
            mets_diff = contains(M.model.mets, blackList);
            % reaction columns
            col_1 = M.model.S(:, index_1);
            col_2 = M.model.S(:, index_2);
            % set blackList entries to 0
            col_1(mets_diff) = 0;
            col_2(mets_diff) = 0;
            % calculate the cosine similarity
            similarity = cosineSimilarity(col_1, col_2);
            opposite = similarity < 0;
        end
    end
    methods (Static)
        function similarity = cosineSimilarity(A, B)
            % Calculate the cosine similarity
            A = full(A);
            B = full(B);
            if isvector(A) && isvector(B)
                similarity = dot( A, B ) / ( norm(A) * norm(B) );
            else
                error('At least one of the inputs is not a vector')
            end
        end
        function row = getRow(matrix, index)
            row = matrix(index, :);
        end
        function array = filterArray(array, filter)
            % set entries of array to zero that have a 1 entry in filter or
            % are contained as indices in filter
            if all(ismember(unique(array), [0 1]))
                % binary array
                if all(ismember(unique(filter), [0 1])) && numel(filter) == numel(array) && numel(array) > 1
                    % filter is binary array of equal length
                    array(logical(filter)) = 0;
                elseif max(filter) <= numel(array)
                    array(filter) = 0;
                else
                    return
                end
            else
                % array contains indices
                if max(filter) <= max(array)
                    % filter is not a binary array and no index in filter
                    % is larger than the maximum in array
                    array = setdiff(array, filter);
                else 
                    return
                end
            end
        end
        

        

    end
end















