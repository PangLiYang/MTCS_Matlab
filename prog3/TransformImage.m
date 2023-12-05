% TransformImage(M,A) applies  transformation M to image A.
% M is a 3 x 3 transformation matrix in homogeneous coordinates
% A is a gray scale image. See p. 174.

function B = TransformImage(M,A);
  white = 255; % gray level for white.
  [n,m] = size(A);
  for i = 1:n
    for j = 1:m
      B(i,j) = white; 
     end 
  end
  for i = 1:n
    for j = 1:m
      v = floor(M*[i;j;1])
      if (1 <= v(1) & v(1) <= n & 1 <= v(2) & v(2) <= m & A(i,j) ~= white) 
        B(v(1),v(2)) = A(i,j);
      end
    end
  end 
end 

