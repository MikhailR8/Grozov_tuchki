function [r, v] = runge(f, mu, r0, v0, step, a)

    arguments
        f function_handle
        mu double
        r0 (3, 1) double
        v0 (3, 1) double
        step double
        a double
    end

    r = r0;
    v = v0;

    [k1r, k1v] = f(r, v, mu, a);

    [k2r, k2v] = f(r + k1r * step / 2, v + k1v * step / 2, mu, a);

    [k3r, k3v] = f(r + k2r * step / 2, v + k2v * step / 2, mu, a);

    [k4r, k4v] = f(r + k3r * step, v + k3v * step, mu, a);

    r = r + (k1r + 2*k2r + 2*k3r + k4r) * step / 6;
    v = v + (k1v + 2*k2v + 2*k3v + k4v) * step / 6;

end