clc; clear;

% Part A

states = States([1, 2]);

% Part B

prob = [0.5, 0.8];
k = size(states, 1);
n = size(states, 2);


function tbl = States(v)

    n = length(v);

    C = {1, n};

    for i = 1 : n
        C{1, i} = 0 : v(i);
    end

    rev_args = flip(C);

    [A{1:n}] = ndgrid(rev_args{:});

    B = cellfun(@(x) x(:), A, 'UniformOutput', false);
    C = flip(B);

    tbl = table(C{:});

    tbl = tbl{:,:}.';
end

function transition = Transition(states, prob)
    
    k = size(states, 1);
    n = size(states, 2);
    coins = states(:, n);
    transition = zeros(n, n);

    for j = 1 : n

        for i = 1 : n

            if j == 1 || j == n
                
            end


        end
    end

end
