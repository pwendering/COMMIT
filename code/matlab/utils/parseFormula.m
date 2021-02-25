function [substrates, products, direction] = parseFormula(formula)
% Parse a formula of the form:
% 
% 'Reaction: M1[.] + M2[.] <=> M3[.] + M4[.]'
% 
%  Input:       char formula:           formula
%  Output:      struct substrates:      contains substrate information in
%                                       fields 'id' (IDs), 'tr'
%                                       (istransported), 'coeff'
%                                       (stoichiometric coefficients)
%               struct products:        contains product information
%               double direction:       0 if reversible, 1 of forward, -1
%                                       if backward
                                        
if (contains(formula, 'rxn') || ~any(regexp(formula, '[mnq]', 'once'))) && contains(formula, {'<=>','<=', '=>'})
    
    if sum(contains(formula, '<=>'))
        direction = 0;
    elseif sum(contains(formula, '=>'))
        direction = 1;
    else
        direction = -1;
    end
    
    if contains(formula, '[1]')
        pattern_tr = '[1]';
    else
        pattern_tr = '[e]';
    end
    
    if any(regexp(formula, 'C[0-9]{5}'))
        pattern = 'C[0-9]{5}';
    elseif any(regexp(formula, 'cpd[0-9]{5}'))
        pattern = 'cpd[0-9]{5}';
    else
        pattern = 'MNXM[0-9]*';
    end
    
    
    formula = regexp(formula, '[^:]+$', 'match');
    split = regexp([formula{:}], '<=>|<=|=>', 'split');

    % substrates
    sub = split{1};
    sub = regexp(sub, '+', 'split');
    substrates.id = regexp([sub{:}], pattern, 'match');
    substrates.tr = contains(sub, pattern_tr)*2;
    substrates.coeff = regexp(sub, '[^\w\d][0-9]+\ ', 'match');
    substrates.coeff(~cellfun('isempty', substrates.coeff)) = [substrates.coeff{~cellfun('isempty', substrates.coeff)}];
    substrates.coeff(cellfun('isempty', substrates.coeff)) = {'1'};
    substrates.coeff = cellfun(@(x)str2double(x),substrates.coeff);
    
    % case of multiple occurrences of the same metabolite
    sub_uniq = unique(substrates.id, 'stable');
    if numel(sub_uniq) < numel(substrates.id)
        % get the inidices of the duplicate metabolites
        sub_mult = find(cell2mat(cellfun(@(x)sum(ismember(substrates.id, x))>1, sub_uniq, 'un',0)));
        for i=1:numel(sub_mult)
            j = sub_mult(i);
            idx = find(ismember(substrates.id, substrates.id{j}));
            if numel(unique(substrates.tr(idx))) == 1
                substrates.coeff(j) = sum(substrates.coeff(idx));
                substrates.coeff(idx(2:end)) = [];
                substrates.tr(idx(2:end)) = [];
                substrates.id(idx(2:end)) = [];
            end
        end
    end
    clear sub sub_uniq sub_mult idx j
    
    %products
    prod = split{2};
    prod = regexp(prod, '+', 'split');
    products.id = regexp([prod{:}], pattern, 'match');
    products.tr = contains(prod, pattern_tr);
    products.coeff = regexp(prod, '\ [0-9]+\ ', 'match');
    products.tr = contains(prod, pattern_tr)*2;
    products.coeff(~cellfun('isempty', products.coeff)) = [products.coeff{~cellfun('isempty', products.coeff)}];
    products.coeff(cellfun('isempty', products.coeff)) = {'1'};
    products.coeff = cellfun(@(x)str2double(x),products.coeff);
    
    % case of multiple occurrences of the same metabolite
    prod_uniq = unique(products.id, 'stable');
    if numel(prod_uniq) < numel(products.id)
        % get the inidices of the duplicate metabolites
        prod_mult = find(cell2mat(cellfun(@(x)sum(ismember(products.id, x))>1, prod_uniq, 'un',0)));
        for i=1:numel(prod_mult)
            j = prod_mult(i);
            idx = find(ismember(products.id, products.id{j}));
            if numel(unique(products.tr(idx))) == 1
                products.coeff(j) = sum(products.coeff(idx));
                products.coeff(idx(2:end)) = [];
                products.tr(idx(2:end)) = [];
                products.id(idx(2:end)) = [];
            end
        end
    end
    clear prod prod_uniq prod_mult idx j
    
    %{
    if sum(substrates.tr)~=0
        tr_p = products.tr;
        tr_p(ismember(products.id, substrates.id(logical(substrates.tr)))) = 1;
    end
    
    if sum(products.tr)~=0
        substrates.tr(ismember(substrates.id, products.id(logical(products.tr)))) = 1;
    end
    
    if exist('tr_p', 'var')
        products.tr = tr_p;
    end
    %}
    
    
else
    substrates.id = {};
    substrates.coeff = [];
    substrates.tr = [];
    products.id = {};
    products.coeff = [];
    products.tr = [];
    direction = 0;
end