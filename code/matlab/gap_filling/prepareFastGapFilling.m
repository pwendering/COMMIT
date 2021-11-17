function dbModel = prepareFastGapFilling(databaseFile, compartments, blackList)
% Parses a database file into a stoichiometric matrix (works for KEGG,
% ModelSEED and MNXref). The file has to be in a form:
% 
% 'Reaction: M1[.] + M2[.] <=> M3[.] + M4[.]'
% 
% Input:
%           char databaseFile:          path to the correctly formatted
%                                       database file
%           cellstr compartments:       compartment ids of the model
%           cellstr blackList:          is used to exclude reactions that 
%                                       contain metabolites in this list
% Output:   struct dbModel:             model structure that contains the
%                                       stoichiometric matrix of the
%                                       universal database in 'S', the
%                                       corresponding metabolites and
%                                       reaction in 'mets' and 'rxns', 
%                                       the binary information of transport in a
%                                       logical array 'transport' and upper
%                                       and lower bounds on the reactions
%                                       in 'ub' and 'lb'

% read the database file
fprintf('\nLoading the database file.\n')
DB = importdata(databaseFile);

% remove blackList reactions and metabolites
blackList = strcat(blackList, '[');
DB = DB(~contains(DB, blackList));

% Extract all metabolites and reactions from the database
if any(cell2mat(regexp(DB, 'C[0-9]{5}', 'once')))
    fprintf('\nRecognized KEGG namespace.\n')
    mets = columnVector(regexp(DB, 'C[0-9]{5}', 'match'));
elseif any(cell2mat(regexp(DB, 'cpd[0-9]{5}', 'once')))
    fprintf('\nRecognized ModelSEED namespace.\n')
    mets = columnVector(regexp(DB, 'cpd[0-9]{5}', 'match'));
elseif any(cell2mat(regexp(DB, 'MNXM[0-9]*', 'once')))
    fprintf('\nRecognized MNXRef namespace.\n')
    mets = columnVector(regexp(DB, 'MNXM[0-9]*', 'match'));
else
    error('No namespace recognized.')
end

% The total number of metabolites is later used for the number of non-zero
% entries in the sparse matrix
n_mets = numel([mets{:}]);

% Construct the arrays that contain the reaction ans metabolite identifiers
% of the database model:
mets = unique([mets{:}]);
rxns = strtok(DB, ':');

dbRxns = rxns;%repmat({''}, numel(rxns)*numel(compartments), 1);
dbMets = mets;%repmat({''}, numel(mets)*numel(compartments), 1);

% If loop_max == number of compartments, all reactions will be inserted in
% all compartments. For simplicity, metabolic reactions are only inserted
% into the cytosol, as it is the only compartment in many bacterial GEMs.
loop_max = 1;
ncomps = numel(compartments);

% Add metabolites and reactions in all compartments
for i=1:ncomps %loopmax
    % reactions
    if contains(compartments{i}, {'e', 'e0'})
        continue
    end
    start_idx = (i-1)*numel(rxns) + 1;
    end_idx = i*numel(rxns);
    dbRxns(start_idx:end_idx) = strcat(rxns, '_', compartments{i});
end

% metabolites have to be in all compartments
for i=1:ncomps
    % metabolites
    start_idx = (i-1)*numel(mets) + 1;
    end_idx = i*numel(mets);
    dbMets(start_idx:end_idx) = strcat(mets, '[', compartments{i}, ']');
end
clear start_idx end_idx rxns mets

% Initialize structures for transport reactions:
% set the pattern for transport reactions
if any(contains(DB, '[1]'))
    pattern_tr = '[1]';
else
    pattern_tr = '[e]';
end

% pre-allocate space for transport reaction matrix
n_tr = sum(contains(DB, pattern_tr))*0.5*ncomps*(ncomps-1);
% clear DB

% matrix that will contain all transport reactions (will be added to the
% matrix of metabolic reactions)
transportRxns = single(zeros(numel(dbMets), n_tr));
% array that collects all reaction identifiers of transport reactions
rxns_tr = repmat({''}, n_tr, 1);
% All transport reactions are made reversible
lb_tr = repmat(-1000,n_tr, 1);
ub_tr = repmat(1000, n_tr, 1);
% Specifies the current column of the transportRxns matrix
rev_col = 0;

% pre-allocate space for the universal stoichiometric matrix
dbModel.S = spalloc(numel(dbMets), numel(dbRxns), n_mets);
% All metabolic reactions are made reversible, will be changed later if
% a reaction is unidirectional
dbModel.lb = repmat(-1000,numel(dbRxns), 1);
dbModel.ub = repmat(1000, numel(dbRxns), 1);

% pairs of compartments that are connected, transport reactions will be
% inserted between all pairs; the second identifier must be unique because
% it is added to the reaction id of the respective transport reaction
if any(contains(compartments, 'p'))
    connectedComps = { {'e', 'p'}, {'p', 'c'} };
else
    connectedComps = { {'c', 'e'}, {'c', 'm'}, {'c', 'r'}, {'c', 'x'}};
end

% The file containing the database is now read line by line
fprintf('\nStarting to add reactions to the universal matrix...\n')
% fid = fopen(databaseFile, 'r');
% dbLine=fgetl(fid);
i = 0;

DB = cellstr(DB);
for dbLine = DB'
% while ischar(dbLine)%for i=1:numel(dbRxns)/loop_max %numel(compartments)
    i = i + 1;
    if mod(i,1000)==0
        fprintf('\nProcessed %d out of %d reactions\n', i, numel(dbRxns))
    end
    % parse the current line
    [s, p, d] = parseFormula(char(dbLine));
    
    % the reaction is a transport reaction if it contains at least one
    % '[1]'
    is_transport = sum([s.tr, p.tr]);
    
    % if coefficients are larger than ~20, the reaction will likely be
    % lumped or lead to incorrectly high production of a compound
    if any(s.coeff>20) || any(p.coeff>20)
        s.id = {};
        p.id = {};
    elseif ~is_transport && ~isempty(intersect(s.id, p.id))
        % check if the products and substrates are identical, if yes it results
        % from lumping different generic compartments and the reaction has to
        % be excluded
        if all(isequal(s.id, p.id))
            s.id = {};
            p.id = {};
        % if not all metabolites are equal but one except water or hydrogen appears 
        % on both sides, exclude the reaction
        elseif ~all(ismember(intersect(s.id, p.id), {'MNXM2', 'MNXM1'}))
            s.id = {};
            p.id = {};
        end
    end

    % if metabolites have been parsed correctly and it is not a transport reaction
    % but it can still be a boundary reaction that would have either no
    % products or no substrates
    if (~isempty(s.id) || ~isempty(p.id)) &&  ~is_transport
        % check directionality and change boundaries accordingly
        if d == 1
            dbModel.lb(i) = 0;
        elseif d == -1
            dbModel.ub(i) = 0;
        end
        
        % get database indices of all metabolites involved in the reaction 
        % (as many a compartments, respectively)
        if ~isempty(s.id) && ~isempty(p.id)
%             s.id = strcat(s.id, '[', compartments(loop_max), ']');
%             p.id = strcat(p.id, '[', compartments(loop_max), ']');
        % Boundary reactions:
        elseif isempty(s.id)
            p.id = strcat(p.id, '[e]');
            dbRxns{i} = regexprep(dbRxns{i}, '_c', '_e');
        else
            s.id = strcat(s.id, '[e]');
            dbRxns{i} = regexprep(dbRxns{i}, '_c', '_e');
        end
        
        % get the row indices of the involved metabolites
        idx_s = cellfun(@(x)find(contains(dbMets, x)), s.id, 'UniformOutput', false);
        idx_p = cellfun(@(x)find(contains(dbMets, x)), p.id, 'UniformOutput', false);
        
        % if all metabolites indices could be obtained from the database,
        % add the corresponding coefficients to the matrix
        if ~isempty(idx_s) || ~isempty(idx_p)
            if ~isempty(idx_s) && all(~cellfun('isempty',idx_s))
                idx_s = cell2mat(idx_s);
            else 
                idx_s = cell.empty(0,0);
            end

            if ~isempty(idx_p) && all(~cellfun('isempty',idx_p))
                idx_p = cell2mat(idx_p);
            else
                idx_p = cell.empty(0,0);
            end
            
            % in case there are multiple compartments for metabolic
            % reactions, they will be inserted here:
            for j=1:ncomps
                % Create column for reaction and insert it into S
                col_idx = i + (j-1)*numel(dbRxns)/ncomps;
                row_idx_p = j:ncomps:(ncomps*numel(p.id));
                row_idx_s = j:ncomps:(ncomps*numel(s.id));
                new_col = zeros(numel(dbMets), 1);
                if ~isempty(idx_p)
                    if isempty(idx_s)
                        new_col(idx_p) = p.coeff';
                    else
                        new_col(idx_p(row_idx_p)) = p.coeff';
                    end
                end
                if ~isempty(idx_s)
                    new_col(idx_s(row_idx_s)) = new_col(idx_s(row_idx_s)) - s.coeff';
                end
                
                dbModel.S(:, col_idx) = new_col;
            end

        else
            warning('%s: Reaction could not be parsed correctly, not all identifiers in from the equation could be found in the database',...
                dbRxns{i})
        end
        
    % Transport reaction:
    elseif is_transport
        % Insert transport reactions between all pairs specified earlier
        for j=1:numel(connectedComps)
            if all(ismember(connectedComps{j}, compartments))
                
                % increase the column index
                rev_col = rev_col + 1;
                
                % add the second compartment id to the reaction identifier 
                % (this way, the reaction identifiers of multiple equal transport
                % reactions can be distinguished)
                rxn_id = strtok(dbRxns{i}, '_');
                rxn_id = strcat(rxn_id, '_', connectedComps{j}{2});
                rxns_tr{rev_col} = rxn_id;
                
                % Here, the compartment identifiers for the transported
                % metabolites are set.
                if ~isempty(s.id) && ~isempty(p.id)
                    id_s = s.id;
                    id_p = p.id;
                    id_s(s.tr==2) = strcat(s.id(s.tr==2), '[', connectedComps{j}{2}, ']');
                    id_s(s.tr~=2) = strcat(s.id(s.tr~=2), '[', connectedComps{j}{1}, ']');
                    id_p(p.tr==2) = strcat(p.id(p.tr==2), '[', connectedComps{j}{2}, ']');
                    id_p(p.tr~=2) = strcat(p.id(p.tr~=2), '[', connectedComps{j}{1}, ']');
                    
                    % find the indices of the metabolites (+comp.) in database
                    idx_s = cellfun(@(x)find(contains(dbMets, x)),id_s, 'UniformOutput', false);
                    idx_p = cellfun(@(x)find(contains(dbMets, x)),id_p, 'UniformOutput', false);
                    clear id_s id_p
                    
                    % write the entries to the database
                    if all(~cellfun('isempty',idx_s)) && all(~cellfun('isempty',idx_p))
                        transportRxns(cell2mat(idx_s), rev_col) = -s.coeff;
                        transportRxns(cell2mat(idx_p), rev_col) = p.coeff;
                    else
                        warning('%s: Reaction could not be parsed correctly, not all identifiers in from the equation could be found in the database',...
                            dbRxns{i})
                    end
                    
                else
                    warning('%s: Exchange or unbalanced reaction, will not be included in universal model', dbRxns{i})
                end

            end
        end
    end
    % next line
%     dbLine = fgetl(fid);
end

% fclose(fid);
% Add transport reactions for all metabolites between all compartments and
% concatenate the transport data structures with the ones for metabolic
% reactions
dbModel.transport = vertcat(zeros(size(dbModel.S,2), 1), ones(n_tr, 1));
dbModel.S = [dbModel.S, double(transportRxns)];
dbModel.lb = vertcat(dbModel.lb, lb_tr);
dbModel.ub = vertcat(dbModel.ub, ub_tr);
dbModel.rxns = vertcat(dbRxns, rxns_tr);
dbModel.mets = columnVector(dbMets);

% delete all-zero columns from S and change the other fields accordingly
non_zero_idx = any(dbModel.S, 1);
dbModel.S = dbModel.S(:,non_zero_idx);
dbModel.lb = dbModel.lb(non_zero_idx);
dbModel.ub = dbModel.ub(non_zero_idx);
dbModel.rxns = dbModel.rxns(non_zero_idx);
dbModel.transport = logical(dbModel.transport(non_zero_idx));

% delete all-zero rows from S and change the other fields accordingly
non_zero_idx = any(dbModel.S, 2);
dbModel.S = dbModel.S(non_zero_idx,:);
dbModel.mets = dbModel.mets(non_zero_idx);

% Since all default biomass reactions have been removed from the universal
% database, the metabolite 'BIOMASS' has to be added to be able to add the
% biomass reaction from the model lateron
dbModel.S(end+1, :) = 0;
dbModel.mets(end+1) = {'BIOMASS'};

end
