clc; clear;
format short g;

% Part A

cps = ProbSuccess(0.2, 0.7, 0.4, 1, 1);
disp("Part A");
disp("Given ProbSuccess(0.2, 0.7, 0.4, 1, 1), we have");
disp(cps);

% Part B

ev = PubValue1(50000, 10000, 500, 0.2, 0.7, 0.4, 1, 1);
disp("Part B");
disp("Given PubValue1(50000, 10000, 500, 0.2, 0.7, 0.4, 1, 1), we have");
disp(ev);

% Part C

aggEv = PubValue2(50000, 10000, 200, 0.2, 0.7, 0.4, 10);
disp("Part C");
disp("Given PubValue2(50000, 10000, 500, 0.2, 0.7, 0.4, 2), we have");
disp(aggEv);

% Part D
opt = OptimalN(50000, 10000, 200, 0.2, 0.7, 0.4);
disp("Part D");
disp("Given OptimalN(50000, 10000, 500, 0.2, 0.7, 0.4), we have");
disp(opt);

function P = ProbSuccess(probSuc, probPosSuc, probPosFail, Npos, Nneg)

    % k -> Npos; q -> Nneg
    pskq = probSuc * (probPosSuc) ^ Npos * (1 - probPosSuc) ^ Nneg;
    pfkq = (1 - probSuc) * (probPosFail) ^ Npos * ((1 - probPosFail)) ^ Nneg;
    P = pskq / (pskq + pfkq);
end

function V = PubValue1(ValueSuc, ValueFail, Fee, probSuc, probPosSuc, probPosFail, Npos, Nneg)

    ps = ProbSuccess(probSuc, probPosSuc, probPosFail, Npos, Nneg);
    n = Npos + Nneg;

    exp = ValueSuc * ps - ValueFail * (1 - ps);

    % Check if decide to publish (*reviewers fees regarded as sunk cost)
    threshold = -n * Fee;

    if exp > threshold
        V = ["True", exp + threshold];
    else
        V = ["False", threshold];
    end
end

function V = PubValue2(ValueSuc, ValueFail, Fee, probSuc, probPosSuc, probPosFail, N)

    aggExp = 0;
    probA = probPosSuc * probSuc + probPosFail * (1 - probSuc);

    for i = 0 : N
        % calculate the Prob of getting i approves
        comb = nchoosek(N, i);
        pips = comb * probPosSuc ^ i * (1 - probPosSuc) ^ (N - i);
        pipf = comb * probPosFail ^ i * (1 - probPosFail) ^ (N - i);
        pipt = probSuc * pips + (1 - probSuc) * pipf;

        ps = ProbSuccess(probSuc, probPosSuc, probPosFail, i, N - i);
        exp = ValueSuc * ps - ValueFail * (1 - ps);

        if exp < 0
            exp = 0;
        end

        aggExp = aggExp + pipt * exp;
    end

    V = aggExp - N * Fee;
end

function O = OptimalN(ValueSuc, ValueFail, Fee, probSuc, probPosSuc, probPosFail)
    val = 0;
    reviewers = 0;

    for i = 0 : probSuc * ValueSuc / Fee
        exp = PubValue2(ValueSuc, ValueFail, Fee, probSuc, probPosSuc, probPosFail, i);
        if (exp > val)
            val = exp;
            reviewers = i;
        end
    end

    O = [val, reviewers];
end



