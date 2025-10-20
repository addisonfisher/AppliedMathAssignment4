function [CFOBT_C, CFOBT_A, CFOBT_B]  = classic_fourth_order()
    CFOBT_C = [0,   1/2, 1/2, 0];
    CFOBT_A = [0,   0,   0,   0;
               1/2, 0,   0,   0;
               0,   1/2, 0,   0;
               0,   0,   1,   0];
    CFOBT_B = [1/6, 1/3, 1/3, 1/6];
end