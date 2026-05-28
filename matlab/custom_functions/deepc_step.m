function [u0, pred, diagnostics] = deepc_step( ...
    U_p, Y_p, U_f, Y_f, ...
    u_ini, y_ini, ...
    r_y, r_u, ...
    Qy, R, ...
    lambda_g, lambda_y, ...
    u_min, u_max)

    % number of columns of hankel matrix
    n_col = size(U_p, 2);

    m = size(R, 1);
    p = size(Qy, 1);

    L_ref = size(U_f, 1) / m;

    g = sdpvar(n_col, 1);
    u_pred = sdpvar(m*L_ref, 1);
    y_pred = sdpvar(p*L_ref, 1);
    sigma_y = sdpvar(size(Y_p, 1), 1);

    % block cost matrices
    Qbar = kron(eye(L_ref), Qy);
    Rbar = kron(eye(L_ref), R);

    % constraints
    Constraints = [];

    Constraints = [Constraints, U_p*g == u_ini];

    % Output past matching with slack
    Constraints = [Constraints, Y_p*g == y_ini + sigma_y];

    % Future prediction
    Constraints = [Constraints, U_f*g == u_pred];
    Constraints = [Constraints, Y_f*g == y_pred];

    % Input constraints
    if ~isempty(u_min) && ~isempty(u_max)
        if isscalar(u_min)
            u_min = u_min * ones(m, 1);
        end
        if isscalar(u_max)
            u_max = u_max * ones(m, 1);
        end

        Constraints = [Constraints, repmat(u_min, L_ref, 1) <= u_pred <= repmat(u_max, L_ref, 1)];
    end

    Objective = ...
        (y_pred - r_y)' * Qbar * (y_pred - r_y) + ...
        (u_pred - r_u)' * Rbar * (u_pred - r_u) + ...
        lambda_g * (g' * g) + ...
        lambda_y * (sigma_y' * sigma_y);

    options = sdpsettings( ...
        'solver', 'quadprog', ...
        'verbose', 0);

    diagnostics = optimize(Constraints, Objective, options);

    if diagnostics.problem ~= 0
        warning("DeePC optimization failed: %s", diagnostics.info);
        u0 = zeros(m, 1);
        pred.u = [];
        pred.y = [];
        pred.g = [];
        pred.sigma_y = [];
        return;
    end

    u_seq = value(u_pred);
    y_seq = value(y_pred);

    % first input
    u0 = u_seq(1:m);

    % store predictions
    pred.u = reshape(u_seq, m, []).';
    pred.y = reshape(y_seq, p, []).';
    pred.g = value(g);
    pred.sigma_y = value(sigma_y);
end