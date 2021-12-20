%% Run COMMIT on small community using manually-curated models

% set seed for random numers
rng('default');
warning off

dbFile = 'data/gap-filling/database/Universal-model-MNXref-balanced.mat';
load(dbFile)

% lower limit for biomass reaction
epsilon = 1e-3;

% weights
% reactions with sequence evidence (E <= 1E-6)
weights.sequence = 25;
% transport reactions
weights.transport = 50;
% transport reactions operating on highly-permeable metabolites
weights.permeable = 50;
% general weight for metabolic reactions
weights.metabolic = 50;
% weight for making a reaction reversible
weights.reverse = 50;
% weight for exchange reactions from the database
weights.exchange = 50;
% weight for allowed exchange / uptake reactions (from previous solutions)
weights.uptake = 1;
% weight for sink reactions if they should be determined in the gap-filling program
weights.sink = 0;
% weights for reactions already contained in the model
weights.model = 0;

% sequence similarity workspace
seq_sim_workspace = 'data\gap-filling\sequence-similarity\sequence_similarity_dv_mm.mat';
load(seq_sim_workspace);

% sequence similarity weights
dbModel_MNXref_balanced.scores = ones(size(dbModel_MNXref_balanced.rxns));
dbModel_MNXref_balanced.genes = repmat({''}, size(dbModel_MNXref_balanced.rxns));

D_vulgaris = readCbModel('data\reference-models\Desulfovibrio_vulgaris_Flowers_2018\iJF744.xml');
M_maripaludis = readCbModel('data\reference-models\Methanococcus_maripaludis_Richards_2016\iMR539.mat');

% add id field to match with column headers of sequence similarity
D_vulgaris.id = 'D. vulgaris';
M_maripaludis.id = 'M. maripaludis';

% correct format of reaction compartment info in M. maripaludis model and
% remove numbers from compartment ids in reactions and metabolite identifiers
M_maripaludis.rxns = regexprep(M_maripaludis.rxns,'\[(?<comp>\w)\d\]','_$<comp>');
M_maripaludis.mets = regexprep(M_maripaludis.mets,'\[(?<comp>\w)\d\]','[$<comp>]');

% remove universally-blocked reactions from both models
D_vulgaris = reduceModel(D_vulgaris,1e-6,false,false,false,true,false);
M_maripaludis = reduceModel(M_maripaludis,1e-6,false,false,false,true,false);

% add rxnNotes field
D_vulgaris.rxnNotes = repmat({''},size(D_vulgaris.rxns));
M_maripaludis.rxnNotes = repmat({''},size(D_vulgaris.rxns));

% translate models from ModelSEED to MNXref
D_vulgaris_transl = translateModel(D_vulgaris, 'ModelSEED', 'MNXref', [], false, true);   
M_maripaludis_transl = translateModel(M_maripaludis, 'ModelSEED', 'MNXref', [], false, true);

nIter = 50;
removePerc = [1 2 5 10];
precInd_dv = nan(nIter,numel(removePerc));
recallInd_dv = nan(nIter,numel(removePerc));
precInd_mm = nan(nIter,numel(removePerc));
recallInd_mm = nan(nIter,numel(removePerc));
precCond_dv = nan(nIter,numel(removePerc));
recallCond_dv = nan(nIter,numel(removePerc));
precCond_mm = nan(nIter,numel(removePerc));
recallCond_mm = nan(nIter,numel(removePerc));

% get shared medium (intersection of default media)
upt_dv = D_vulgaris_transl.rxns((sum(D_vulgaris_transl.S==-1)==1)&(sum(D_vulgaris_transl.S~=0)==1));
upt_dv = cellfun(@(x)regexp(x,'cpd\d{5}','match'),upt_dv);
upt_mm = M_maripaludis_transl.rxns((sum(M_maripaludis_transl.S==-1)==1)&(sum(M_maripaludis_transl.S~=0)==1));
upt_mm = cellfun(@(x)regexp(x,'cpd\d{5}','match'),upt_mm,'un',0);
upt_mm = [upt_mm{:}]';

medium = translateIDs(intersect(upt_dv,upt_mm),'met',[],'ModelSEED','MNXref');
medium = setdiff(medium,'BIOMASS');
medium = strcat(medium,'[e]');
medium = [medium; {'MNXM78334[c]';'MNXM51417[c]';'MNXM81014[c]'}];

auxo_media = {...
    strcat(translateIDs(upt_dv,'met',[],'ModelSEED','MNXref',false),'[e]'),...
    [strcat(translateIDs(upt_mm,'met',[],'ModelSEED','MNXref',false),'[e]');...
    {'MNXM78334[c]';'MNXM51417[c]';'MNXM81014[c]'}]...
    };

% only allow those reactions to be removed that contain at least one
% metabolite that is contained in the gap-filling database
allowed_rem_dv = findRxnsFromMets(D_vulgaris_transl,...
    D_vulgaris_transl.mets(ismember(D_vulgaris_transl.mets,dbModel_MNXref_balanced.mets)));
allowed_rem_dv = setdiff(allowed_rem_dv,[D_vulgaris_transl.rxns(D_vulgaris_transl.c==1);...
    D_vulgaris_transl.rxns(contains(D_vulgaris_transl.rxns,'EX_'));...
    D_vulgaris_transl.rxns(~contains(D_vulgaris_transl.rxns,'MNXR'))]);

allowed_rem_mm = findRxnsFromMets(M_maripaludis_transl,...
    M_maripaludis_transl.mets(ismember(M_maripaludis_transl.mets,dbModel_MNXref_balanced.mets)));
allowed_rem_mm = setdiff(allowed_rem_mm,[M_maripaludis_transl.rxns(M_maripaludis_transl.c==1);...
    M_maripaludis_transl.rxns(contains(M_maripaludis_transl.rxns,'EX_'));...
    M_maripaludis_transl.rxns(~contains(M_maripaludis_transl.rxns,'MNXR'));...
    {'rxn13784_c'; 'rxn13783_c'; 'rxn13782_c'}]);

%% Run for N iterations
for i=1:nIter
    fprintf('------------------------------------------------------\n')
    fprintf('Iteration #%d\n',i)
    fprintf('------------------------------------------------------\n')
    for j=1:numel(removePerc)
        % remove a random set of X% of reactions from both models
        removeIdx_dv = randsample(numel(allowed_rem_dv),...
            floor(numel(allowed_rem_dv)*removePerc(j)/100));
        D_vulgaris_tmp = removeRxns(D_vulgaris_transl,...
            allowed_rem_dv(removeIdx_dv),'metFlag', false);
        
        removeIdx_mm = randsample(numel(allowed_rem_mm),...
            floor(numel(allowed_rem_mm)*removePerc(j)/100));
        M_maripaludis_tmp = removeRxns(M_maripaludis_transl,...
            allowed_rem_mm(removeIdx_mm),'metFlag', false);
        
        % run individual gap-filling
        dbModel_MNXref_balanced.scores = seq_sim_mat(:,1);
        dbModel_MNXref_balanced.genes = genes{1};
        [~, ~, ~, addedRxns_dv] = ...
            condFastGapFilling(D_vulgaris_tmp, dbModel_MNXref_balanced, [],...
            weights, epsilon, false, false);
        
        dbModel_MNXref_balanced.scores = seq_sim_mat(:,2);
        dbModel_MNXref_balanced.genes = genes{2};
        [~, ~, ~, addedRxns_mm] = ...
            condFastGapFilling(M_maripaludis_tmp, dbModel_MNXref_balanced, [],...
            weights, epsilon, false, false);
        
        % extract plain reaction IDs and compound IDs from exchange
        % reactions for both removed and added reactions
        rem_dv = [...
            regexp(allowed_rem_dv(removeIdx_dv),'MNXR\d+','match');...
            regexp(allowed_rem_dv(removeIdx_dv),'rxn\d+','match');...
            regexp(allowed_rem_dv(removeIdx_dv),'cpd\d+','match')...
            ];
        rem_dv = [rem_dv{:}];
        
        add_dv = [...
            regexp(addedRxns_dv,'MNXR\d+','match');...
            regexp(addedRxns_dv,'MNXM\d+','match');...
            regexp(addedRxns_dv,'rxn\d+','match');...
            regexp(addedRxns_dv,'cpd\d+','match')...
            ];
        add_dv = [add_dv{:}];
        
        % caluclate precision and recall
        if ~isempty(add_dv)
            precInd_dv(i,j) = numel(intersect(rem_dv,add_dv)) / numel(add_dv);
            recallInd_dv(i,j) = numel(intersect(rem_dv,add_dv)) / numel(rem_dv);
        end
        
        % extract plain reaction IDs and compound IDs from exchange
        % reactions for both removed and added reactions
        rem_mm = [...
            regexp(allowed_rem_mm(removeIdx_mm),'MNXR\d+','match');...
            regexp(allowed_rem_mm(removeIdx_mm),'rxn\d+','match');...
            regexp(allowed_rem_mm(removeIdx_mm),'cpd\d+','match')...
            ];
        rem_mm = [rem_mm{:}];
        
        add_mm = [...
            regexp(addedRxns_mm,'MNXR\d+','match');...
            regexp(addedRxns_mm,'MNXM\d+','match');...
            regexp(addedRxns_mm,'rxn\d+','match');...
            regexp(addedRxns_mm,'cpd\d+','match')...
            ];
        add_mm = [add_mm{:}];
        
        % caluclate precision and recall
        if ~isempty(add_mm)
            precInd_mm(i,j) = numel(intersect(rem_mm,add_mm)) / numel(add_mm);
            recallInd_mm(i,j) = numel(intersect(rem_mm,add_mm)) / numel(rem_mm);
        end
        
        % conditional gap-filling
        models = {D_vulgaris_tmp,M_maripaludis_tmp};
        try
            modelsGF = iterativeGapFilling(models, medium, auxo_media,...
                dbModel_MNXref_balanced, weights, epsilon,...
                false, [], 2, seq_sim_workspace);
            
            addedRxns_dv = modelsGF{1}.rxns(cellfun(@(x)ismember({'gf'},x),modelsGF{1}.rxnNotes));
            
            % extract plain reaction IDs and compound IDs from exchange
            % reactions for both removed and added reactions
            add_dv = [...
                regexp(addedRxns_dv,'MNXR\d+','match');...
                regexp(addedRxns_dv,'MNXM\d+','match');...
                regexp(addedRxns_dv,'rxn\d+','match');...
                regexp(addedRxns_dv,'cpd\d+','match')...
                ];
            add_dv = [add_dv{:}]';
            
            % caluclate precision and recall
            if ~isempty(add_dv)
                precCond_dv(i,j) = numel(intersect(rem_dv,add_dv)) / numel(add_dv);
                recallCond_dv(i,j) = numel(intersect(rem_dv,add_dv)) / numel(rem_dv);
            end
            
            addedRxns_mm = modelsGF{2}.rxns(cellfun(@(x)ismember({'gf'},x),modelsGF{2}.rxnNotes));
            
            % extract plain reaction IDs and compound IDs from exchange
            % reactions for both removed and added reactions
            add_mm = [...
                regexp(addedRxns_mm,'MNXR\d+','match');...
                regexp(addedRxns_mm,'MNXM\d+','match');...
                regexp(addedRxns_mm,'rxn\d+','match');...
                regexp(addedRxns_mm,'cpd\d+','match')...
                ];
            add_mm = [add_mm{:}];
            
            % caluclate precision and recall
            if ~isempty(add_mm)
                precCond_mm(i,j) = numel(intersect(rem_mm,add_mm)) / numel(add_mm);
                recallCond_mm(i,j) = numel(intersect(rem_mm,add_mm)) / numel(rem_mm);
            end
        catch
            disp('Iterative gap-filling not successful!')
        end
    end
end

%% Create Barplots
colorsHex = strrep({'#999999','#E69F00'},'#','');
colorsDec = cell2mat(cellfun(@hex2dec,regexp(colorsHex,'\w{2}','match'),'un',0))';
colorsDec = colorsDec/255;

subplot(1,2,1)
label_text = [strcat(strsplit(strtrim(sprintf('-%d%% ',removePerc))),' (ind.)')...
    strcat(strsplit(strtrim(sprintf('-%d%% ',removePerc))),' (cond.)')];
labels = categorical(label_text);
labels = reordercats(labels,label_text);
b = bar(labels,...
    [mean(precInd_dv,'omitnan')', mean(recallInd_dv,'omitnan')';...
    mean(precCond_dv,'omitnan')', mean(recallCond_dv,'omitnan')']);
for i=1:2 b(i).FaceColor = colorsDec(i,:); end

ylim([0 0.6])
legend({'Precision','Recall'},'Box','off')
ylabel('Precision / Recall')
text(.05,.94,'A','FontSize',16,'FontWeight','bold','Units','normalized')
% title('\it D. vulgaris','FontSize',24)
box on
set(gca,'FontSize',14,'LineWidth',1.3)

subplot(1,2,2)
b = bar(labels,...
    [mean(precInd_mm,'omitnan')', mean(recallInd_mm,'omitnan')';...
    mean(precCond_mm,'omitnan')', mean(recallCond_mm,'omitnan')']);
for i=1:2 b(i).FaceColor = colorsDec(i,:); end

ylim([0 0.6])
legend({'Precision','Recall'},'Box','off')
ylabel('Precision / Recall')
text(.05,.94,'B','FontSize',16,'FontWeight','bold','Units','normalized')
% title('\it M. Maripaludis','FontSize',24)
box on
set(gca,'FontSize',14,'LineWidth',1.3)

print('commit_ref_community.png','-dpng','-painters')
save('ref_community_ws')

%% Only predict secreted metabolites separately given original models

D_vulgaris = readCbModel('data\reference-models\Desulfovibrio_vulgaris_Flowers_2018\iJF744.xml');
M_maripaludis = readCbModel('data\reference-models\Methanococcus_maripaludis_Richards_2016\iMR539.mat');

% enable L-lactate uptake
D_vulgaris = changeRxnBounds(D_vulgaris,'EX_cpd00159(e)',-1000,'l');
D_vulgaris = changeRxnBounds(D_vulgaris,'EX_cpd00159(e)',1000,'u');
% find permeable metabolites
comps = cellfun(@(x)regexp(x,'\[.\]$','match'),D_vulgaris.mets);
D_vulgaris.mets = translateIDs(strtok(D_vulgaris.mets,'['),'met',[],'ModelSEED','MNXref');
D_vulgaris.mets = strcat(D_vulgaris.mets,comps);
ex_d_vulgaris = findPotentialExcMets(D_vulgaris,1);
% translate IDs to names
ex_d_vulgaris_names = translateIDs(strtok(ex_d_vulgaris,'['),'met',[],'MNXref','NAMES');

% enable formate update
M_maripaludis = changeRxnBounds(M_maripaludis,'EX_cpd00047[e0]',-1000,'l');
M_maripaludis = changeRxnBounds(M_maripaludis,'EX_cpd00047[e0]',1000,'u');

comps = cellfun(@(x)regexp(x,'\[.0\]$','match'),M_maripaludis.mets);
M_maripaludis.mets = translateIDs(strtok(M_maripaludis.mets,'['),'met',[],'ModelSEED','MNXref');
M_maripaludis.mets = strcat(M_maripaludis.mets,comps);

M_maripaludis.mets = strrep(M_maripaludis.mets,'[c0]','[c]');
M_maripaludis.mets = strrep(M_maripaludis.mets,'[e0]','[e]');

ex_m_maripaludis = findPotentialExcMets(M_maripaludis,1);
ex_m_maripaludis_names = translateIDs(strtok(ex_m_maripaludis,'['),'met',[],'MNXref','NAMES');

fprintf('Unique to D. vulgaris:\n')
disp(setdiff(ex_d_vulgaris_names,ex_m_maripaludis_names));

fprintf('Unique to M. maripaludis:\n')
disp(setdiff(ex_m_maripaludis_names,ex_d_vulgaris_names));



