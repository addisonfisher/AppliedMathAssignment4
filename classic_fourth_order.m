function CFOBT = classic_fourth_order()
    CFOBT = struct();

    CFOBT.C = [0,   1/2, 1/2, 0];
    CFOBT.A = [0,   0,   0,   0;
               1/2, 0,   0,   0;
               0,   1/2, 0,   0;
               0,   0,   1,   0];
    CFOBT.B = [1/6, 1/3, 1/3, 1/6];
end