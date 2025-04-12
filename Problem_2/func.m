function [res_r, res_v] = func(r, v, mu, a)

    J2 = 1.5 * 1082.8*10^(-6) * mu * a.^2 * [(5 * r(1) * r(3).^2 / norm(r).^2) - r(1);
                                         (5 * r(2) * r(3).^2 / norm(r).^2) - r(2);
                                         (5 * r(3).^3 / norm(r).^2) - 3*r(3)] / norm(r).^5;
    res_v = (- mu * r / norm(r).^3);
    res_r = v;
    
end