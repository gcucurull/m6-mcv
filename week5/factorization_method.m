function [Pproj, Xproj] = factorization_method(x1,x2)

    [norm_x1, T1] = normalise2dpts(x1);
    [norm_x2, T2] = normalise2dpts(x2);
    
    lambda = ones (2, size(x1,2));
    
    d= Inf;
    flag = true;
    while flag
        
        for i=1:2
            for j = 1:2
                lambda(j,:) = lambda(j,:) / norm(lambda(j,:));
            end
            for j = 1:size(x1,2)
                lambda(:,j) = lambda(:,j) / norm(lambda(:,j));
            end
        end
        
        
        
        M = zeros(3*2, size(x1,2));
        M(1,:) = lambda(1,:) .* norm_x1(1,:);
        M(2,:) = lambda(1,:) .* norm_x1(2,:);
        M(3,:) = lambda(1,:) .* norm_x1(3,:);
        M(4,:) = lambda(2,:) .* norm_x2(1,:);
        M(5,:) = lambda(2,:) .* norm_x2(2,:);
        M(6,:) = lambda(2,:) .* norm_x2(3,:);
        
        [U,D,V] = svd(M);
        
        Pm = U * D(:,1:4);
        Xproj = V(:,1:4)';
        
        d_old= d;
       
        for i=1:2
            if i==1
                 Px = Pm(1:3,:) * Xproj;
                 x = norm_x1;
            elseif i==2
                 Px = Pm(4:6,:) * Xproj;
                 x = norm_x2;
            end
            for j=1:size(x1,2)
                 d = d + norm(x(:,j) - Px(:,j))^2; 
            end
        end
        d = d / (2 * size(x1,2));
        if ~((abs(d - d_old)/d) < 0.1)
            flag = false;
        end
    end
    
    Pproj(1:3,:) = inv(T1) * Pm(1:3,:);
    Pproj(4:6,:) = inv(T2) * Pm(4:6,:);
    
end