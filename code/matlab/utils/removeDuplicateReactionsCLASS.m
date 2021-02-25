classdef removeDuplicateReactionsCLASS < handle
    properties
        model
        dbCheck
        database
        verbose
        
        epsilon = 10E-6;
        t_identity = 1-epsilon;
        t_identity_lower = 0.95;
        count_identity_lower = 0;
        count_identity = 0;
        count_diff_rev = 0;
        count_opposite = 0;
        count_added = 0;
        
    end
    
    methods (model, dbCheck, database, verbose)
        
        
        
        
        function similarity = cosineSimilarity(A, B)
            % Calculate the cosine similarity
            if isvector(A) && isvector(B)
                similarity = dot( A, B ) / ( norm(A) * norm(B) );
            else
                error('At least one of the inputs is not a vector')
            end
        end
        
        function similarity_matrix = getSimilarityMatrix(model)
            % Cosine similarity matrix
            n = numel(model.rxns);
            stoich_mat = full(model.S);
            similarity_matrix = zeros(n);
            
            parfor i=1:n-1
                A = stoich_mat(:,i);
                row_tmp = zeros(1,n);
                for j=i+1:n
                    B = stoich_mat(:,j);
                    row_tmp(j) = cosineSimilarity(A,B);
                end
                similarity_matrix(i,:) = row_tmp;
            end
        end
        
        function model = removeByStoichiometry(model, dbCheck, database, verbose)
            % get the cosine similarity between each reaction pair
            if verbose
                fprintf('\nFind duplicate reactions by stoichiometry\n')
                fprintf('\n\t> Calculating the Cosine similarity matrix... ')
            end
            tic
            similarity_matrix = getSimilarityMatrix(model);
            if verbose; fprintf('done. (%.2fs)\n', toc); end
            
            
            
        end
    end
end