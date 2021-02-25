function smiles = getSMILESFromCID(cid)
% PUG REST URL request

if iscell(cid)
    cid = char(cid);
elseif isdouble(cid)
    cid = num2str(cid);
end

[status, response] = ...
    unix(['wget -q -O - https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',...
    cid, '/property/CanonicalSMILES/TXT']);

if status
    warning('Non-zero exit status')
    smiles = {''};
else
    smiles = strtrim(response);
end

end