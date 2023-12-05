clc; clear;

% Part A

states = States([1, 2]);

% Part B

prob = [0.5, 0.8];
transition = Transition(states, prob);

% Part C

stationary = Stationary(transition);
disp(stationary)
% Prat D

probFullUrn = ProbFullUrn(transition, 6);
probFullUrn2 = ProbFullUrn(transition, 20);
disp(probFullUrn2)

% Part E

expTimeToFill = ExpTimeToFill(transition, 0.001);

function states = States(v)

    n = length(v);

    C = {1, n};

    for i = 1 : n
        C{1, i} = 0 : v(i);
    end

    rev_args = flip(C);

    [A{1:n}] = ndgrid(rev_args{:});

    B = cellfun(@(x) x(:), A, 'UniformOutput', false);
    C = flip(B);

    states = table(C{:});

    states = states{:,:}.';
end

function transition = Transition(states, prob)

    n = size(states, 2); 
    full = states(:,6);
    fullCoins = sum(full);
    transition = zeros(n, n);

    for i = 1 : n
        urn = states(:,i);
        urnCoins = sum(urn);
        out = full - urn;
        outCoins = sum(out);

        for j = 1 : n

            diff = states(:,j) - urn;
            posSum = 0;
            negSum = 0;

            for g = 1 : length(diff)
                if diff(g) > 0
                    posSum = posSum + 1;
                elseif diff(g) < 0
                    negSum = negSum + 1;
                end
            end

            diffSum = sum(diff);

            if i == 1
                
                if diffSum == 0
                    transition(j,i) = dot(full,1 - prob) / fullCoins;
                elseif diffSum == 1
                    idx = find(diff == 1);
                    transition(j,i) = full(idx) / fullCoins * prob(idx);
                else
                    transition(j,i) = 0;
                end

                continue
            end

            if i == n
                
                if diffSum == 0
                    transition(j,i) = dot(full, prob) / fullCoins;
                elseif diffSum == -1
                    idx = find(diff == -1);
                    transition(j,i) = full(idx) / fullCoins * (1 - prob(idx));
                else
                    transition(j,i) = 0;
                end

                continue
            end

            if diffSum == 0 & posSum == 0 & negSum == 0
                transition(j,i) = 0.5 * dot(urn, prob) / urnCoins ...
                                + 0.5 * dot(out, 1 - prob) / outCoins;
            elseif diffSum == 1 & posSum == 1 & negSum == 0
                idx = find(diff == 1);
                transition(j,i) = 0.5 * out(idx) / outCoins * prob(idx);
            elseif diffSum == -1 & posSum == 0 & negSum == 1
                idx = find(diff == -1);
                transition(j,i) = 0.5 * urn(idx) / urnCoins * (1 - prob(idx));
            else
                transition(j,i) = 0;

            end

        end
    end

end

function stationary = Stationary(M)
    n = length(M);
    Q = M - eye(n);
    Q(n+1,:) = ones(1,n);
    ret = zeros(n+1,1);
    ret(n+1) = 1;

    stationary = Q \ ret;
end

function probFullUrn = ProbFullUrn(M, N)

    probFullUrn = zeros(1,N);

    U = M;
    n = length(U);
    U(:,n) = zeros(n,1);
    U(n,n) = 1;

    p = zeros(n,1);
    p(1) = 1;

    for i = 1 : N
        p = U * p;
        probFullUrn(i) = p(n);
    end

end

function expTimeToFill = ExpTimeToFill(M, eps)

    U = M;
    n = length(U);
    U(:,n) = zeros(n,1);
    U(n,n) = 1;

    p = zeros(n,1);
    p(1) = 1;

    expTimeToFill = 0;

    p0 = 0;
    idx = 0;
    bound = 1;

    while bound > eps
        idx = idx + 1;
        p = U * p;
        p1 = p(n);

        diff = p1 - p0;
        expTimeToFill = expTimeToFill + idx * diff;
        bound = idx * (1 - p1);
        p0 = p1;
    end

end
