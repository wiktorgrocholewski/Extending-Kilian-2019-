function Theta = irf(A, D0, H)
% Computes IRFs from a VAR(p) model
% A   : K × (K·p) coefficient matrix
% D0  : K × K impact matrix (structural shocks)
% H   : scalar horizon
% Theta : K × K × (H+1) IRFs

    [K, Kp] = size(A);
    p = Kp / K;

    % Build companion matrix F: (Kp x Kp)
    F = zeros(K*p, K*p);
    F(1:K, :) = A;
    F(K+1:end, 1:K*(p-1)) = eye(K*(p-1));

    % J selects the first K rows (current values) from the companion state
    J = [eye(K); zeros(K*(p-1), K)];

    % Allocate IRF array
    Theta = zeros(K, K, H+1);

    for h = 0:H
        Fh = F^h;  % (Kp x Kp)
        Theta(:, :, h+1) = Fh(1:K, :) * J * D0;  % (K x K)
    end
end