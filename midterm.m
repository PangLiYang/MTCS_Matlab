clc; clear;

disp("Problem 10-A")

heat = [0 .25 0 .25 0 0 0 0 0 50/4;
        .25 0 .25 0 .25 0 0 0 0 25/4;
        0 .25 0 0 0 .25 0 0 0 50/4;
        .25 0 0 0 .25 0 .25 0 0 25/4;
        0 .25 0 .25 0 .25 0 .25 0 0;
        0 0 .25 0 .25 0 0 0 .25 25/4;
        0 0 0 .25 0 0 0 .25 0 325/4;
        0 0 0 0 .25 0 .25 0 .25 300/4;
        0 0 0 0 0 .25 0 .25 0 325/4;
        0 0 0 0 0 0 0 0 0 1];

base = [0 0 0 0 0 0 0 0 0 1];

final = heat^1000 * base.';

final = final(1:9);

disp("From t1 to t9:")
disp(final);

disp("Problem 10-B");

disp("The function TempDistribution belows");

ansB = TempDistribution(5, 300, 25, 25, 25);

disp(ansB);

disp("Problem 10-C");

ansC = TempDistribution(3, 300, 25, 25, 25);

disp(ansC);

function a = TempDistribution(n, ta, tb, tc, td)
    a = zeros(n, n);
    
    a(1, :) = tc;
    a(:, 1) = tb;
    a(:, n) = td;
    a(n, :) = ta;

    % Iteration 1000 times (asymptotic to converge)
    m = 1000;

    while m > 0
        for i = 2 : n - 1
            for j = 2 : n -1
                a(i,j) = (a(i-1, j)+a(i+1, j)+a(i, j-1)+a(i, j+1))/4;
            end
        end
        m = m -1;
    end

end

