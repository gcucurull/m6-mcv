function [ dangle ] = angle( l1, l2 )
%ANGLE computes angle between 2 lines

    dangle=acosd( dot(l1,l2)/norm(l1)/norm(l2) ); 

end

