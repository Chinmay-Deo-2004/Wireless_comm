%Finding minimum distance of Linear Block Code using Parity Check Matrix

clear;
H = createParityCheckMatrix();
disp(H);
min_add_cols(H);

function min_cols = min_add_cols(A)
    [m,n] = size(A);
    min_cols = n;
    for k = 1:n
        comb = nchoosek(1:n,k);
    
        for i = 1:size(comb,1)
            xor_result = zeros(m,1);
                for j = 1:k
                    xor_result = mod(xor_result + A(:, comb(i,j)), 2);
                end
                if all(xor_result == 0)
                    min_cols = k;
                    disp(k);
                    return;
                end
        end 
    end
end 

function H = createParityCheckMatrix()
    n = input('Enter value for n:');
    k = input('Input value for k:');
        if n<k
            disp('Invalid. n Should be > k');
            H = createParityCheckMatrix();
        end
    H = [randi([0,1], n-k, k), eye(n-k)];
end
            disp('Invalid. n Should be > k');
            H = createParityCheckMatrix();
        end
    H = [randi([0,1], n-k, k), eye(n-k)];
end
