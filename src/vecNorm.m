function varargout = vecNorm(X, dim, norm_type)
% VECNORM return the norm of the vectors which are along the
% dimension dim. 
%   norms = vecNorm(X, dim, norm_type) returns the norms of the vectors in matrix
%   X. The norm is taken along the dimension dim. Particulary, if dim = 1, 
%   vectors are column-vectors and if dim = 2, vectors are row-vectors. 
%   norm_typecan be anything that matlab's NORM takes (e.g. 1, 2, Inf) or 
%   missing, which defaults to Euclidean (2-Norm).
%
%   [norms, normalizedX] = vecNorm(X, dim, norm_type) also returns a normalized X
%

    if ~exist('norm_type', 'var');
        norm_type = 2;  % euclidean default
    end
    assert(dim == 1 || dim == 2, ...
        'dim can only be 1 (col-vecs) or 2 (row-vecs)');

    % note: this is much faster than looping and calling on norm(.) in
    % cases with low dimensionality and large number of vectors
    if norm_type == Inf
        norms = max(X, [], dim);
    elseif norm_type == -Inf
        norms = min(X, [], dim);
    elseif isnumeric(norm_type) % assuming p-norm
        p = norm_type;
        norms = sum(abs(X) .^ p, dim) .^ (1/p);
    
            
    % for norms like the frobenius norm
    else
        % column-vectors
        if dim == 1
            nr_vecs = size(X, 2);
            norms = zeros(1, nr_vecs);
            for i = 1:nr_vecs
                norms(i) = norm(X(:, i), norm_type);
            end

        % row-vectors
        else
            nr_vecs = size(X, 1);
            norms = zeros(nr_vecs, 1);
            for i = 1:nr_vecs
                norms(i) = norm(X(i, :), norm_type);
            end

        end
    end
    varargout{1} = norms;
    
    % returned normalized vec
    if nargout == 2
        repvec = ones(1, ndims(X));
        repvec(dim) = size(X, dim);
        varargout{2} = X ./ repmat(norms, repvec);
    end
    