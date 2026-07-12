% Code is taken from https://github.com/iamKaiZhang/DeePC/ and modified to
% our needs

classdef DeePC < handle
    properties
        Up, Uf, Yp, Yf,
        y_f,
        Q, R
        ny, nu, ng,
        L_ini, L_ref
        u_max, y_max,
        yalmip_optimizer
    end

    methods
        function obj = DeePCc(Up, Uf, Yp, Yf, y_f, Q, R, u_max, y_max)
            obj.Up = Up;
            obj.Uf = Uf;
            obj.Yp = Yp;
            obj.Yf = Yf;

            obj.y_f = y_f;
            obj.Q = Q;
            obj.R = R;

            obj.ny = size(Q, 1);
            obj.nu = size(R, 1);
            obj.ng = size(Up, 2);

            obj.L_ini = size(Up, 1) / obj.nu;
            obj.L_ref = size(Uf, 1) / obj.nu;
            
            obj.u_max = u_max;
            obj.y_max = y_max;

            obj.init_optimizer();

        end

        function init_optimizer(obj)
            u_var = sdpvar(obj.nu, obj.L_ref, 'full');
            y_var = sdpvar(obj.ny, obj.L_ref, 'full');
            u_ini = sdpvar(obj.nu, obj.L_ini, 'full');
            y_ini = sdpvar(obj.ny, obj.L_ini, 'full');
            g = sdpvar(obj.ng, 1, 'full');

            objective = 0;
            constraints = [
                obj.Up * g == reshape(u_ini, [], 1), ...
                obj.Yp * g == reshape(y_ini, [], 1), ...
                obj.Uf * g == reshape(u_var, [], 1), ...
                obj.Yf * g == reshape(y_var, [], 1), ...
            ]; 

            for i = 1:obj.L_ref
                objective = objective ...
                    + (y_var(:, i)-obj.y_f)'*obj.Q*(y_var(:, i)-obj.y_f) ...
                    + (u_var(:, i))'*obj.R*(u_var(:, i));
                constraints = [constraints,...
                    u_var(:, i) <= obj.u_max * ones(obj.nu, 1), ...
                    y_var(:, i) <= obj.y_max * ones(obj.ny, 1)
                ];
            end

            opts = sdpsettings('verbose', 1, 'solver', 'mosek');
            opts.mosek.MSK_DPAR_INTPNT_QO_TOL_REL_GAP = 1e-8;
            obj.yalmip_optimizer = optimizer( ...
                constraints, objective, opts, ...
                {u_ini, y_ini}, ...  % params: initial condition
                {u_var, y_var, objective} ...  % outputs
            );            
        end

        function [u, info] = solve(obj, u_ini, y_ini)
            [out, errorcode] = obj.yalmip_optimizer({u_ini, y_ini});
            [U, Y, objective] = out{:};
            u = U(:, 1);
            info = struct('errorcode', errorcode, 'objective', objective, 'U', U, 'Y', Y);
        end        
    end
end