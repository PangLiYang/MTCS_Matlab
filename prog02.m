clc; clear;

disp("Problem 1");
% For a 2*k matrix P and its 2*k midpoints matrix Q
% With hint 1, we will have C * Pt = Qt
% Multiple C-1 to each side, we get Pt = C-1 * Qt

Q = [1 3 1 -1;
     4 4 3 6];

[found, P] = InvertMidpoints(Q);

disp(found);
disp(P);

disp("-------------------------------------")
disp("Problem 2.A");

Pa = [0;
      1];
Pb = [4 1 -3;
      4 0 1];
F = ForceMatrix(Pa, Pb);

disp("F = ")
disp(round(F, 4));

Qb = [50;
      -6;
      36];

disp("ans = ")
disp(round(F * Qb, 4));

disp("-------------------------------------")
disp("Problem 2.B");

TF = TotalForce(Pa,Pb,Qb);

disp("TF = ")
disp(round(TF, 4));

C = PossibleCharge(Pa, Pb, TF);
disp("C = ")
disp(round(C, 4))

disp("ans = ")
disp(round(F * C, 4))

disp("-------------------------------------")
disp("Problem 2.C");

Paa = [0 2;
       1 0];

Tf1 = TotalForce(Paa(:, 1), Pb, Qb);
Tf2 = TotalForce(Paa(:, 2), Pb, Qb);

Tff = [Tf1, Tf2];

C = FindCharge(Paa, Pb, Tff);

disp("C = ")
disp(C);




function [found, P] = InvertMidpoints(Q)
    found = "true";
    k = length(Q(1, :));
    C = zeros(k, k);
    e = 1e-8;
    
    for i = 1 : k - 1
        C(i, i) = 1 / 2;
        C(i, i + 1) = 1 / 2;
    end

    C(k, k) = 1 / 2;
    C(k, 1) = 1 / 2;

    if mod(k, 2) == 1
        C = inv(C);
        P = C * Q.';
        P = P.';
    else
        cc = C;
        cc(k, :) = [];
        cc = pinv(cc);
        qq = Q;
        qq(:, k) = [];
        P = cc * qq.';
        P = P.';

        if abs(P(1, 1) + P(1, k) - 2 * Q(1,k)) < e && ...
           abs(P(2, 1) + P(2, k) - 2 * Q(2,k)) < e
            P = round(P, 8);
        else
            found = "false";
            P = zeros(2, k);
        end
    end

end

function F = ForceMatrix(Pa, Pb)
    n = length(Pb(1, :));
    container = ones(2, n);
    for i = 1 : n
        container(:, i) = Pa;
    end

    container = container - Pb;
    distant = container;

    for i = 1 : n
        temp = sqrt(distant(1, i) ^ 2 + distant(2, i) ^ 2) ^ 3;
        distant(:, i) = [temp, temp];
    end

    container = container ./ distant;
    F = container;
end

function TF = TotalForce(Pa, Pb, Qb)
    F = ForceMatrix(Pa, Pb);
    TF = F * Qb;
end

function C = PossibleCharge(Pa, Pb, TF)
    F = ForceMatrix(Pa, Pb);
    C = pinv(F) * TF;
end

function C = FindCharge(Pa, Pb, TF)
    F1 = ForceMatrix(Pa(:, 1), Pb);
    F2 = ForceMatrix(Pa(:, 2), Pb);
    F = vertcat(F1, F2);
    Fv = pinv(F);

    TF = reshape(TF, [4,1]);
    C = Fv * TF;
end


