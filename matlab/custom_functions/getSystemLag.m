function l = getSystemLag(A, C)
    if size(A, 1) ~= size(A, 2)
        error("Matrix A not square")
    elseif size(A, 1) ~= size(C, 2)
        error("Matrix A and C dont match")
    end
    
    n = size(A, 1);
    prev = [];
    for i = 0:n-1
        if (rank([prev; C*A^i]) == rank(obsv(A, C)))
            l = i;
            break;
        end
        prev = [prev; C*A^i];
    end
end