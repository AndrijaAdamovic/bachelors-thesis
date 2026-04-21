% Simpler than built in hankel function, better for MIMO systems
function H = create_hankel(data, depth)
    if isvector(data)
        data = data(:);
    end

    [N, m] = size(data);

    n_cols = N - depth + 1;
    H = zeros(m * depth, n_cols);

    for k = 1:depth
        rows = (k - 1)*m + (1:m);
        H(rows, :) = data(k : k + n_cols - 1, :).';
    end
end