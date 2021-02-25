classdef ReactionPair
    properties
        R1              % reaction 1
        R2              % reaction 2
        similarity      % cosine similarity
    end
    properties (Constant)
        epsilon = mergeModelsClass.epsilon;
        t_identity = mergeModelsClass.t_identity;
        t_identity_lower = mergeModelsClass.t_identity_lower;
    end
    methods
        function rp = ReactionPair(R1, R2, s)
            % constructor for ReactionPair
            if nargin == 0
                rp.R1 = Reaction;
                rp.R2 = Reaction;
                rp.similarity = 0;
            else
                if R1.getIndex ~= R2.getIndex
                    rp.R1 = R1;
                    rp.R2 = R2;
                else
                    error('Both reactions have the same index.')
                end
                if nargin == 3
                    rp.similarity = s;
                end
            end
            
        end
        function remove = decideByReversibility(rp, Reaction_1, Reaction_2)
            % Compares the reversibility of both reactions and chooses to remove the one
            % that reversible. If both are irreversible, no decision can be made.
            % Input:
            %       class instances Reaction_1, Reaction_2
            % Output:
            %       double remove:              index of the reaction that
            %                                   should be removed; 0 if no
            %                                   decision could be made
            if nargin == 1
                Reaction_1 = rp.R1;
                Reaction_2 = rp.R2;
            end
            if Reaction_2.getReversibility
                % the second one or both are reversible, remove the second one
                remove = Reaction_2.getIndex;
            elseif Reaction_1.getReversibility && ~Reaction_2.getReversibility
                % only the first one is reversible, remove the first one
                remove = Reaction_1.getIndex;
            else
                % both are irreversible
                remove = 0;
            end
            
        end
        function remove = decideByMassBalance(rp)
            % Compares the mass balance of both reactions and chooses to remove the one
            % that is not mass-balanced. If both are mass-balanced, no decision can be made.
            % Input:
            %       class instances Reaction_1, Reaction_2
            % Output:
            %       double remove:              index of the reaction that
            %                                   should be removed; 0 if no
            %                                   decision could be made
            Reaction_1 = rp.R1;
            Reaction_2 = rp.R2;
            if Reaction_1.getMassBalance && ~Reaction_2.getMassBalance
                % only the first one is mass-balanced, remove the second one
                remove = Reaction_2.getIndex;
            elseif ~Reaction_1.getMassBalance && Reaction_2.getMassBalance
                % only the second one is mass-balanced, remove the first one
                remove = Reaction_1.getIndex;
            else
                % either both are mass-balanced or both are not mass-balanced
                remove = 0;
            end
        end
        function remove = decideByProtons(rp)
            metDiff = rp.getMetDiff;
            diff_is_proton = rp.containsProton(metDiff) & numel(metDiff)==1;
            protons_1 = rp.R1.containsProton;
            protons_2 = rp.R1.containsProton;
            if diff_is_proton && protons_1
                remove = rp.R2.index;
            elseif diff_is_proton && protons_2
                remove = rp.R2.index;
            else 
                remove = 0;
                
            end
        end
        function i = identity(rp)
            i = (rp.similarity >= (rp.t_identity - rp.epsilon));
        end
        function i = identity_lower(rp)
            s = rp.similarity;
            t_i = rp.t_identity - rp.epsilon;
            t_l = rp.t_identity_lower;
            i = (s >= t_l) & (s < t_i) | (s <= -t_l) & (s > -t_i);
        end
        function sameDirection = compareDirections(rp)
            % Compares the direction of the reactions reports whether they have the same
            % or opposite directions.
            % Input:
            %       class instances Reaction_1, Reaction_2
            % Output:
            %       logical sameDirection:           1 if they have the same direction
            sameDirection = rp.getSimilarity > 0;
        end
        function sameIDs = compareIDs(rp)
            id_1 = rp.R1.id;
            id_2 = rp.R2.id;
            sameIDs = isequal(id_1, id_2);
        end
        function metDiff = getMetDiff(rp)
            mets_1 = rp.R1.mets;
            mets_2 = rp.R2.mets;
            metDiff = setxor(mets_1, mets_2);
        end
        function sameMets = compareMets(rp)
            sameMets = isempty(rp.getMetDiff);
        end
        function sameRev = compareReversibility(rp)
            sameRev = rp.R1.rev & rp.R2.rev;
        end
        function rp = setSimilarity(rp, s)
            % set the reaction cosine similarity
            rp.similarity = s;
        end
        function similarity = getSimilarity(rp)
            % get the cosine similarity value
            similarity = rp.similarity;
        end
        function ids = findMNXRxn(rp)
            Reaction_1 = rp.R1;
            Reaction_2 = rp.R2;
            ids = strcat(Reaction_1.id,  Reaction_2.id);
            ids = regexp(char(ids), 'MNXR\d+_\w', 'match');
        end
        function id = selectMNXref(rp)
            ids = rp.findMNXRxn;
            if ~isempty(ids)
                id = ids(1);
            else
                id = {};
            end
        end
    end
    methods (Static)
        function p = containsProton(mets)
            p = ~isempty(regexp(mets, 'MNXM1\[.+?\]', 'once'));
        end
    end
end

