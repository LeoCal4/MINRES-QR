function [x] = diag_inverter(D, b)
    x = 1./diag(D) .* b;
end
