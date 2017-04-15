function [Pproj, Xproj] = factorization_method(x1,x2, init) %Pproj, Xproj
%       - 'init' can be either "ones" or "sturm"

    [norm_x1, T1] = normalise2dpts(x1);
    [norm_x2, T2] = normalise2dpts(x2);
    
    lambda = ones (2, size(x1,2));
    
    if isequal(init, 'sturm')
        % camera 1
        F1 = fundamental_matrix(x1, x1);
        [U, D, V] = svd(F1);
        e1 = V(:,3) / V(3,3);
            
        % camera 2
        F2 = fundamental_matrix(x2, x1);
        [U, D, V] = svd(F2);
        e2 = V(:,3) / V(3,3);
            
        for j = 1:size(x1,2)
            num = x1(:, j)'*F1*cross(e1, x1(:,j));
            denom = norm(cross(e1, x1(:,j))).^2*lambda(1, j);
            lambda(1,j) = num/denom;
        end
        for j = 1:size(x2,2)
            num = x1(:, j)'*F2*cross(e2, x2(:,j));
            denom = norm(cross(e2, x2(:,j))).^2*lambda(1, j);
            lambda(2,j) = num/denom;
        end
    end
    
    d= Inf;
    flag = true;
    while flag
        
        rescale = true;
        counter = 0;
        lambda_diff = Inf;
        while rescale
            old_lambda_diff = lambda_diff;
            old_lambda = lambda;
            if mod(counter, 2)
                % normalize rows
                lambda(1,:) = lambda(1,:) ./ norm(lambda(1,:));
                lambda(2,:) = lambda(2,:) ./ norm(lambda(2,:));
            else
                % normalize columns
                for col = 1:size(x1,2)
                    lambda(:,col) = lambda(:,col) / norm(lambda(:,col));
                end
            end
            
            % compute euclidan difference with old lambda
            lambda_diff = (old_lambda - lambda).^2;
            lambda_diff = sum(lambda_diff(:));
            lambda_diff = sqrt(lambda_diff);
                        
            if ((abs(lambda_diff - old_lambda_diff)/lambda_diff) < 0.1)
                rescale = false;
            end
            counter = counter +1;
            
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
        d = 0;
       
        for i=1:2
            if i==1
                 Px = Pm(1:3,:) * Xproj;
                 x = norm_x1;
            elseif i==2
                 Px = Pm(4:6,:) * Xproj;
                 x = norm_x2;
            end
            for j=1:size(x1,2)
                 %d = d + sqrt(sum((x(:,j) - Px(:,j)).^2));
                 d = d + sum((x(:,j) - Px(:,j)).^2);
            end
        end
        
        d
        
        if ((abs(d - d_old)/d) < 0.1)
            flag = false;
        else
            % If it has not converged update lambdas
            temp = Pm*Xproj;
            lambda(1,:) = temp(3,:);
            lambda(2,:) = temp(6,:);
        end
    end
    
    Pproj(1:3,:) = inv(T1) * Pm(1:3,:);
    Pproj(4:6,:) = inv(T2) * Pm(4:6,:);







%     Normalize the points
%     [x1, H1] = normalise2dpts(x1);
%     [x2, H2] = normalise2dpts(x2);
%     x{1} = x1;
%     x{2} = x2;
% 
%     Get F and the epipole
%     [F, ~] = ransac_fundamental_matrix(x1, x2, 2.0); 
%     [~, ~, V] = svd(F');
%     e = V(:,end);
% 
%     lambda1 = ones(1,length(x1));
%     lambda2 = [];
% 
%     for i=1:length(x1)
%         lambda2(i) = (x1(:,i)' * F' * cross(e,x2(:,i)))/ norm(cross(e,x2(:,i)))^2;
%     end
% 
%     d_old = 0;
%     d = 15;
%     while (abs(d - d_old)/d) >= 0.1
%         4_Alternate rescaling the rows of the depth matrix lambda(2x24)(formed by lambdas) to
%         have unit norm and the columns of lambda(2x24) to have unit norm until lambda(2x24)
%         stops changing significantly (usually two loops).
%         lambda = [lambda1;lambda2]; % 2x24
%         for j=1:2
%             normalize per row
%             for i=1:size(lambda,2)
%                 a = sum(lambda(:,i));
%                 lambda(:,i) = lambda(:,i) ./ a;
%             end
%             normalize per column
%             for i=1:size(lambda,1)
%                 a = sum(lambda(i,:));
%                 lambda(i,:) = lambda(i,:) ./ a;
%             end
%         end
%         lambda1 = lambda(1,:);
%         lambda2 = lambda(2,:);
% 
%         build matrix M and decompose
%         M = [lambda1(1,:).*x1(1,:);
%             lambda1(1,:).*x1(2,:);
%             lambda1(1,:).*x1(3,:);
%             lambda2(1,:).*x2(1,:);
%             lambda2(1,:).*x2(2,:);
%             lambda2(1,:).*x2(3,:)];
% 
%         [U,D,V] = svd(M);
% 
%         P and X
%         Pproj = U*D(:,1:4); % 3m x 4 = 6x4
%         P{1} = Pproj(1:3,:);
%         P{2} = Pproj(4:6,:);
%         V_T = V';
%         X = V_T(1:4,:); % 4 x n = 4x24
% 
%         Convergence
%         d_old = d;
%         d = 0;
%         for j=1:2
%             m = euclid(P{j}*X);
%             m1 = euclid(x{j});
%             d = d + sum(sum(norm(m1(:,:) - m(:,:)).^2));
%         end
% 
%         Normalize 3D point
%         X = X ./ repmat(X(end,:), size(X,1), 1);
% 
%         lambda1 = P{1}*X;
%         lambda1 = lambda1(end,:);
%         lambda2 = P{2}*X;
%         lambda2 = lambda2(end,:);
%     end
% 
%     De-normalise
%     Pproj(1:3,:) = H1\P{1};
%     Pproj(4:6, :) = H2\P{2};
end