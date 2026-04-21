function [u, t, u_prbs] = generate_filtered_prbs(Ts, duration, amplitude, hold_samples, tau_f)
    N = round(duration / Ts);
    t = (0:N-1)' * Ts;

    n_blocks = ceil(N / hold_samples);
    prbs_blocks = 2 * randi([0, 1], n_blocks, 1) - 1;
    u_prbs = repelem(prbs_blocks, hold_samples);
    u_prbs = amplitude * u_prbs(1:N);

    %Filtering the digital signal -> more realistic     
    if tau_f <= 0
        u = u_prbs;
    else
        alpha = Ts / (tau_f + Ts);
        u = zeros(N, 1);
        u(1) = u_prbs(1);

        for k = 2:N
            u(k) = (1 - alpha) * u(k-1) + alpha * u_prbs(k);
        end
    end
end