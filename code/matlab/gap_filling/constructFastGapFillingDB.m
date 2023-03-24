% parse the MNXref universal reaction database and add metabolite
% permeability information to construct a gap-filling database

%% Preparation
% load options
options

compartments = {'c', 'e'};

blackList = {'BIOMASS',...
    'MNXM517',...       % Light/hnu/photon
    'MNXM162292',...    % Light/hnu/photon
    'MNXM22409',...     % UV light
    'MNXM1093995',...   % photons
    'MNXM145920',...    % Photon (281 to 306 nm, UVB)
    'MNXM145816',...    % Photon (378 to 482 nm, violet/blue)
    'MNXM145712',...    % Photon (380 to 750 nm, visible spectrum)
    'MNXM1093982',...   % Photon (400nm-420nm)
    'MNXM145817',...    % Photon (406 to 454 nm, indigo/blue)
    'MNXM145796',...    % Photon (417 to 472 nm, indigo/blue)
    'MNXM1093983',...   % Photon (420nm-440nm)
    'MNXM145818',...    % Photon (451 to 526 nm, blue/green)
    'MNXM1093984',...   % Photon (460nm-480nm)
    'MNXM1093985',...   % Photon (500nm-520nm)
    'MNXM1093986',...   % Photon (520nm-540nm)
    'MNXM1093987',...   % Photon (540nm-560nm)
    'MNXM1093988',...   % Photon (560nm-580nm)
    'MNXM1093989',...   % Photon (580nm-600nm)
    'MNXM1093990',...   % Photon (600nm-620nm)
    'MNXM145673',...    % Photon (608 to 666 nm, orange/red)
    'MNXM1093991',...   % Photon (620nm-640nm)
    'MNXM1093992',...   % Photon (640nm-660nm)
    'MNXM145690',...    % Photon (659 to 684 nm, red)
    'MNXM1093993',...   % Photon (660nm-680nm)
    'MNXM145691',...    % Photon (662 to 691 nm, red)
    'MNXM1093994',...   % Photon (680nm-700nm)
    'MNXM680566',...    % photon_b
    'MNXM1093977',...   % Photon energy loss (heat/fluorescence)
    'MNXM1092938',...   % Photons with 680nm wavelength
    'MNXM1092939',...   % Photons with 700nm wavelength
    'MNXM639',...       % peptide
    'MNXM55103'...      % glucan
    };

% MNXref universal database filtered for balanced reactions
databaseFile = 'data/tables/MNXref/reaction_MNXref_balanced.lst';

% parse database file
dbModel_MNXref_balanced = prepareFastGapFilling(databaseFile, compartments, blackList);

% add metabolite names and permeability info
dbModel_MNXref_balanced = addNamesPermeabilityToGfDb(...
    dbModel_MNXref_balanced,...
    'data/tables/MNXref/');
% save the output as a matlab workspace
save('data/gap-filling/database/Universal-model-MNXref-balanced.mat',...
    'dbModel_MNXref_balanced')



