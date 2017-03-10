function  [H] = homography2d(x1, x2)
    [x1norm, T1]=normalization(x1);
    [x2norm, T2]=normalization(x2);
    x2 = x2norm(1,:);
    y2 = x2norm(2,:);
    w2 = x2norm(3,:);
    A = [];
    for i=1:size(x1norm,2)
        A = [A; zeros(3,1)'     -w2(i)*x1norm(:,i)'   y2(i)*x1norm(:,i)'; ...
                w2(i)*x1norm(:,i)'   zeros(3,1)'     -x2(i)*x1norm(:,i)'];
    end
    [U,D,V] = svd(A);
    H = reshape(V(:,9),3,3)';
    H = inv(T2) * H * T1;
end