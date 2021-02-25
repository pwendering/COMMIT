function [ec, notes] = correctEC(ecList, ecTranslation)
% Correct a list of Enzyme Commission (EC) numbers using a translation table 
% or the KEGG REST API.
% 
% Input:
%               cell ecList:        List of EC numbers
%               table ecTranslation:Table that contains all EC numbers with
%                                   their respective corrections
% Output:
%               cell ec:            List of corrected EC numbers
%               cell notes:         Information about tranfer

% Check input
if ~iscell(ecList)
    error('The EC input list in not of type cell')
end

if nargin > 1
    if ~istable(ecTranslation)
        error('The input for ecTranslation is not of type table')
    end
else
    ecTranslation = [];
end
% Initialize output variables
ec = cell(numel(ecList),1);
notes = cell(numel(ecList),1);

% Regex pattern for EC number (EC numbers containing letters are not
% considered)
if isempty(ecTranslation)
    pattern = '[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+';
else
    pattern = '[0-9A-Z]+\.[0-9A-Z]+\.[0-9A-Z]+\.[0-9A-Z]+';
end

dlm = '|';
        
for i=1:numel(ecList)
    query = regexp(ecList{i}, pattern, 'match');
    if isempty(query)
        ec{i} = ecList{i};
        notes{i} = '';
        continue
    end
    query = query(~cellfun('isempty', query));
    % Some reactions can have multiple EC numbers assigned to them
    entry = cell(1,numel(query));
    entry(:) = {''};
    note = cell(1,numel(query));
    note(:) = {''};
    if isempty(ecTranslation)
        for j=1:numel(query)
            % Multiple lines of response possible that are received as a
            % string. I split the string by newlines to get a cell array
            response = strtrim(urlread(fullfile('http://rest.kegg.jp/find/enzyme/',...
                query{j})));
            response = strsplit(response, '\n');
            % Only this line is selected that contains exactly the desired
            % EC number
            response = regexp(response, strcat('ec:',query{j},'\s.*'), 'match');
            response = response(~cellfun('isempty', response));
            response = [response{:}];
            if ~isempty(response) && ~isempty(regexpi(response{:}, 'Transferred'))
                % delimiter is '\t', select the field that contains the transfer information
                response = regexpi(strsplit(response{:}, '\t'),...
                    strcat('Transferred to .*'),'match');
                response = response(~cellfun('isempty', response));
                response = [response{:}];
                % Extract the new EC number(s)
                response = regexp(response{:}, pattern, 'match');
                if numel(response) > 1
                    entry{j} = char(strjoin(string(response), dlm));
                else
                    entry{j} = response{:};
                end
                note{j} = 'transferred';
            else
                entry{j} = query{j};
                note{j} = '';
            end
        end
    else
        for j=1:numel(query)
            p = regexprep(query{j},'\.', '\\.');
            idx = regexp(ecTranslation{:,1}, strcat(p,'$'));
            idx = ~cellfun('isempty', idx);
            
            if sum(idx) == 0
                entry{j} = query{j};
            else
                entry(j) = table2cell(ecTranslation(idx,2));
                
                if ~isequal(entry{j},query{j})
                    note{j} = 'transferred';
                end
            end
        end
    end
    
    % If the reaction has multiple EC assignments, join them by a delimiter
    if numel(query) > 1
        entry = entry(~cellfun('isempty', entry));
        
        % Check is all EC numbers are equal now
        if all(cellfun(@(x)isequal(x,entry{1}),entry))
            ec{i} = entry{1};
        else
            ec{i} = strjoin(entry(:), dlm);
        end
        
        if contains(note, 'transferred')
            notes{i} = strjoin(note(:), dlm);
        else
            notes{i} = note{:};
        end
        
    else
        ec{i} = entry{:};
        notes{i} = note{:};
    end
    
end
end





