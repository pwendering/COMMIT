function new_ec = completeEC(ec)
% If an EC number has less than four levels, the missing ones are replaced
% with '-'
% Input:
%           cellstr ec:         array of E.C. numbers

if ~char(ec)
    error('Input must be of type char')
end

level = numel(regexp(ec, '\.')) + 1;
if  level < 4
    new_ec = strcat(ec, repmat('.-',1,4-level));
else
    new_ec = ec;
end

end