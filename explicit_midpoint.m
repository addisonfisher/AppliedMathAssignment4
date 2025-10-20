function EMBT = explicit_midpoint
    EMBT = struct();

    EMBT.C = [0,   1/2];
    EMBT.A = [0,   0
              1/2, 0];
    EMBT.B = [0,   1];
end