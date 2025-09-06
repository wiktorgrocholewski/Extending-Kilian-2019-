function [D, lambda] = decompose_hetero(sigma1, sigma2)
    % Ensure sigma1 is symmetric positive definite
    if ~isequal(sigma1, sigma1') || any(eig(sigma1) <= 0)
        error('Matrix sigma1 must be symmetric positive definite.');
    end

    % Step 1: Cholesky factorization of sigma1 = L*L'
    L = chol(sigma1, 'lower');

    % Step 2: Transform sigma2 into the sigma1 basis
    S = L \ (sigma2 / L');  % Equivalent to inv(L) * sigma2 * inv(L)'

    % Step 3: Eigen-decomposition of the transformed matrix
    [V, D_lambda] = eig(S);  % V is orthogonal, D_lambda is diagonal

    % Step 4: Recover D and lambda
    D = L * V;
    lambda = D_lambda;

    % Optional: sort eigenvalues and reorder columns of D accordingly
    [vals, idx] = sort(diag(lambda));
    lambda = diag(vals);
    D = D(:, idx);
end