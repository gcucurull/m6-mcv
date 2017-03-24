function [F_es] = fundamental_matrix(x1_test, x2_test)
    [x1Norm, t1] = normalise2dpts(x1_test); 
    [x2Norm, t2] = normalise2dpts(x2_test); 
    x1 = x1Norm(1,:)';
    y1 = x1Norm(2,:)';
   
    x2 = x2Norm(1,:)';
    y2 = x2Norm(2,:)';

    a = [x2.*x1 x2.*y1 x2 y2.*x1 y2.*y1 y2 x1 y1 ones(8,1)];

    [u, d, v] = svd(a);

    F_es = v(:,9);
    F_es = [F_es(1) F_es(2) F_es(3); F_es(4) F_es(5) F_es(6); F_es(7) F_es(8) F_es(9)];

    [u, d, v] = svd(F_es);
    d(3,3) = 0;
    F_es = u * d * v';

    F_es = t2' * F_es * t1;
end

