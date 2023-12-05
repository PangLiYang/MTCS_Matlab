clc;
clear;
warning('off');

disp("1.");
disp(EulerSieve1(27));

disp("2.");
disp(EulerSieve2(27));

disp("3.");
P.px = [0 2 4 -2 6];
P.py = [0 8 0 6 6];
[qx, qy] = ConnectMidPoints(P);
disp(qx);
disp(qy);

disp("4.");
P.px = 7 + [1 3 5 7 2.5 4 6].*cos(10 * pi * (1 : 7) / 7);
P.py = 7 + [1 3 5 7 2.5 4 6].*sin(10 * pi * (1 : 7) / 7);
ConvergingPolygons(P, 6);


function l = EulerSieve1(N)
    l = 2 : N;
    p = 2;
    while p^2 < N
        l1 = p : (N / p);
        l2 = p * l1;
        l = setdiff(l, l2);
        idx = find(l == p);
        p = l(idx + 1);
    end
end

function l = EulerSieve2(N)
    l = ones(1, N);
    l(1) = 0;
    p = 2;
    while p^2 < N

        % We actually donnot need l1 here (will be redundant)
        % l1 = zeros(1, N);
        % for idx = p : (N / p)
        %     l1(idx) = 1;
        % end

        l2 = zeros(1, N);
        for idx = p : (N / p)
            l2(p * idx) = 1;
        end

        l = l - l2;
        idx = p;

        for i = p + 1 : N
            if l(i) == 1
                idx = i;
                break;
            end
        end

        if (idx == p)
            break;
        end

        p = idx;
    end

    l = max(l,0);
end

function [qx, qy] = ConnectMidPoints(P)
    n = length(P.px);
    qx = zeros(1, n);
    qy = zeros(1, n);

    for i = 1 : n
        if i == n
            qx(i) = (P.px(i) + P.px(1)) / 2;
            qy(i) = (P.py(i) + P.py(1)) / 2;
        else
            qx(i) = (P.px(i + 1) + P.px(i)) / 2;
            qy(i) = (P.py(i + 1) + P.py(i)) / 2;
        end
    end
end

function ConvergingPolygons(P, N)
    n = length(P.px);

    for c = 1 : N

        subplot(N, 2, 2 * c - 1);
        pgon = polyshape(P.px, P.py);
        pg = plot(pgon);
        pg.FaceColor = 'white';
        axis('equal');

        subplot(N, 2, 2 * c);
        pgon = polyshape(P.px, P.py);
        [qx, qy] = ConnectMidPoints(P);
        P.px = qx;
        P.py = qy;

        hold on;
        pg = plot(pgon);

        for j = 1 : n
            plot(P.px(j), P.py(j), 'r.');
        end
        hold off;

        pg.FaceColor = 'white';
        axis('equal');
    end
end
