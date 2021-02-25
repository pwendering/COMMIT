function new_gpr = mergeGPR(gpr_1, gpr_2)
%% Combine GPR rules for two reactions in COBRA format
% example :
%               GPR = 'x(1) | ( x(2) & x(3) )'
% Input:
%           cellstr/char gpr_1, gpr_2:      GPR rules for the two reactions
% Output:
%           cellstr new_gpr:                new GPR rule
%% Steps:
%   1. find all components that are connected by OR or AND on the first
%      level
%   2. take the union of all equal components
%   3. if complex rules are connected by OR, they are compared recursively
%      by splitting the component by the main operator

if iscellstr(gpr_1)
    gpr_1 = char(gpr_1);
elseif ~ischar(gpr_1)
    error('The first rule is has an incorrect type')
end

if iscellstr(gpr_2)
    gpr_2 = char(gpr_2);
elseif ~ischar(gpr_2)
    error('The second rule is has an incorrect type')
end

% First, remove outer brackets around the rule if present
gpr_1 = regexprep(gpr_1, '(^\()|(\D\)$)', '');
gpr_2 = regexprep(gpr_2, '(^\()|(\D\)$)', '');
if isequal(strtrim(gpr_1), strtrim(gpr_2))
    % equal rules
    new_gpr = gpr_1;
elseif ~contains([gpr_1, gpr_2], '&')
    % only OR rules
    % ==> join all genes by OR
    genes = regexp([gpr_1, gpr_2], 'x\(\d+\)', 'match');
    new_gpr = strjoin(unique(genes), ' | ');
else
    % add the first rule to the new rule
    new_gpr = gpr_1;
    % Split the rule into it components
    components_1 = parseGPRrule(gpr_1);
    components_2 = parseGPRrule(gpr_2);
    
    if isempty(components_1) || isempty(components_2)
        error('No components found for at least one of the rules!')
    end
    
    % compare first rule components to new_gpr element-wise
    for i=1:numel(components_2)
        genes_comp_2 = regexp(components_2{i}, 'x\(\d+\)', 'match');
        
        if any(contains(components_1, components_2(i)))
            % element is exact match ot part of another element
            continue
            %         elseif numel(genes_comp_2) == 1
            %             % single genes and not contained as single in the other
            %             % rule
            %             % ==> add to new rule with OR
            %             new_gpr = char(strcat({new_gpr}, {' | '}, components_2(i)));
        else
            hit = false;
            for j=1:numel(components_1)
                genes_comp_1 = regexp(components_1{j}, 'x\(\d+\)', 'match');
                if all(ismember(genes_comp_1, genes_comp_2)) && numel(genes_comp_1) == numel(genes_comp_2)
                    % same set of genes
                    if numel(unique(regexp(components_1{j}, '[\|&]', 'match'))) == 1  && ...
                            numel(unique(regexp(components_2{i}, '[\|&]', 'match'))) == 1
                        % only AND or only OR
                        % ==> consider them as equal
                        hit = true;
                        continue
                    elseif contains(components_1(j), '&') && ~contains(components_2(i), '&') % numel(regexp(components_1, '\&')) == numel(genes_comp_1) - 1
                        % only the first rule contains at least one AND operator
                        % ==> prefer the AND-joined rule
                        hit = true;
                        continue
                    elseif contains(components_2(i), '&') && ~contains(components_1(j), '&')
                        % only the second rule contains at least one AND operator
                        % ==> prefer the AND-joined rule
                        new_gpr = char(regexprep(new_gpr,...
                            regexptranslate('escape', components_1(j)), components_2{i}));
                        hit = true;
                        break
                    else
                        % both contain AND and contain the same genes
                        % ==> split recursively until identical or not
                        comps_split_1 = parseGPRrule(components_1(j), 0);
                        comps_split_2 = parseGPRrule(components_2(i), 0);
%                         disp(components_1(j))
%                         disp(components_2(i))
                        for k=1:numel(comps_split_2)
                            for l=1:numel(comps_split_1)
                                fprintf('Recursive split\n')
                                new_comp = mergeGPR(comps_split_2(k), comps_split_1(l));
                                new_gpr = char(strcat({new_gpr}, {' | '}, char(new_comp)));
                            end
                        end
                        hit = true;
                    end
                end
            end
            if ~hit
                % no hit was found
                % ==> add this component with an OR
                new_gpr = char(strcat({new_gpr}, {' | '}, components_2(i)));
            end
        end
    end
    
    % now check for nested components
    components = parseGPRrule(new_gpr);
    remove = [];
    for i=1:numel(components)
        if (sum(contains(components, components(i))) > 1)
            % if the component is found in another component, remove it
            % but not if is AND-joined
            remove = vertcat(remove, i);
        end
    end
    components(remove) = [];
    new_gpr = strjoin(components, ' | ');
    
end
new_gpr = {new_gpr};
end

