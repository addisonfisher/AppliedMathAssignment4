function [EMBT_C, EMBT_A, EMBT_B] = explicit_midpoint
    EMBT_C = [0,   1/2];
    EMBT_A = [0,   0
              1/2, 0];
    EMBT_B = [0,   1];
end