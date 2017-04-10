function [v1] = vanishing_point(xo1, xf1, xo2, xf2)
    v1 = cross(cross(xo1, xf1),cross(xo2, xf2)); 
end