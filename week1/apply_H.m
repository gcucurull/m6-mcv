function [ out_im ] = apply_H( im, T )
% Receives an image 'im' and a matrix transofrmation 'T'
    
    out_im=zeros(size(im));
    
    from = im;
    to = out_im;
    
    % based on http://dsp.stackexchange.com/questions/21703/how-to-apply-a-2d-2d-homography-matrix-to-an-image
    % It is quite bad, should be improved
    for i=1:size(out_im,1)
        for j=1:size(out_im,2)
            v1 = [j;i;1];
            v2 = T*v1;
            v3 =v2/v2(3,1);
            if (v3(1,1) <0 || v3(2,1)<0 || v3(3,1)<0)
               continue;
            else     
                to(round(1+v3(2,1)),round(1+v3(1,1)),:)=from(i,j,:);
            end
        end
    end
    
    out_im = to;

end

