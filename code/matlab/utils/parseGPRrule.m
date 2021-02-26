function components = parseGPRrule(rule, splitAND)
%% split a GPR rule into the first-level components
% Input:
%   cellstr/char rule:          gene-protein-reaction (GPR) rule
%   logical splitAND:           whether or not the GPR rule should be split
%                               by AND operator(s)
% Output:
%   cellstr components:         main components of the GPR rule

if nargin < 2
    splitAND = 0;
end

if iscellstr(rule)
    rule = char(rule);
end

components = {};
match = {'start'};
next_operator = '|';
while ~isempty(match)
    % find the next component
    pos_complex = regexp(rule, '\(\D.+?\D\)', 'once');
    pos_single = regexp(rule, 'x\(\d+\)',  'once');
    if isempty(pos_complex) || (pos_single < pos_complex)
        match = regexp(rule, 'x\(\d+\)', 'match', 'once');
    else
        match = regexp(rule, '\(\D.+?\D\)', 'match', 'once');
    end
    
    if ~isempty(match)
        % add it to the array of components
        if isequal(next_operator, '&') && ~splitAND
            % join components of AND-joined component
            % split both the current and previous component by OR
            % and join all components elementwise with AND
            split_current = regexp(match, '\ ?\|\ ?', 'split');
            split_previous = regexp(components{end}, '\ ?\|\ ?', 'split');
            components(end) = [];
            for i=1:numel(split_previous)
                for j=1:numel(split_current)
                    components = vertcat(components, strjoin([split_previous(i), split_current(j)], ' & '));
                end
            end
        else
            components = vertcat(components, match);
        end
        % remove the component from the rule
        rule = regexprep(rule, regexptranslate('escape', match), '', 'once');
        % get the next operator and remove it from the rule as well
        next_operator = regexp(rule, '[\|&]', 'match', 'once');
        rule = regexprep(rule, '[\|&]', '', 'once');
    end
end

% remove brackets from components
components = strtrim(regexprep(components, '[^\w]\)|\([^\w]', ''));

end