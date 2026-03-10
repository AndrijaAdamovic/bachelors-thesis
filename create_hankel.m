function H = create_hankel(data, depth)
    H = zeros(depth, length(data) - depth + 1);
    for k = 1:depth
        H(k, :) = data(k : (length(data) - depth + k));
    end
end