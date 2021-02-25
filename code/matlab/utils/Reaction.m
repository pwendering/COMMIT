classdef Reaction
    % Reaction object in a metabolic model
    properties
        id              % identifier
        index           % index of the reaction
        mets            % substrates and products
        rev             % reversibility
        dir             % direction (->: 1, <=>: 0, <-: -1)
        massBalance     % mass balance of substrates and products
        lb              % lower bound
        ub              % upper bound
    end

    methods
        function rxn = Reaction(id, index, mets, lb, ub, massBalance)
            % Constructor method
            if nargin == 0
                rxn.id = '';
                rxn.index = 0;
                rxn.mets = {};
                rxn.lb = -1000;
                rxn.ub = 1000;
                rxn.massBalance = 1;
            else
                rxn.id = id;
                rxn.index = index;
                rxn.mets = mets;
                rxn.lb = lb;
                rxn.ub = ub;
                rxn.massBalance = massBalance;
            end
            
            % reversibility
            rxn.rev = rxn.getReversibility;
            % direction
            rxn.dir = rxn.getDirection;
            
        end
        function rev = getReversibility(rxn)
            rev = rxn.getLower && rxn.getUpper;
        end
        function rxn = createRevID(rxn)
            rxn.id = strcat(rxn.id, '_r');
        end
        function dir = getDirection(rxn)
            dir = sign(rxn.getUpper) + sign(rxn.getLower);
        end
        function lb = getLower(rxn)
            lb = rxn.lb;
        end
        function ub = getUpper(rxn)
            ub = rxn.ub;
        end
        function index = getIndex(rxn)
            index = rxn.index;
        end
        function mb = getMassBalance(rxn)
            mb = rxn.massBalance;
        end
        function rxn = setMassBalance(rxn, mb)
            rxn.massBalance = mb;
        end
        function p = containsProton(rxn)
            p = ~isempty(regexp(rxn.mets, 'MNXM1\[.+?\]', 'once'));
        end
    end
end

