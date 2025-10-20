function [SSPRK3BT_C, SSPRK3BT_A, SSPRK3BT_B] = strong_stability_preserving_RK()
    SSPRK3BT_C = [0,   1,   1/2];
    SSPRK3BT_A = [0,   0,   0
                  1,   0,   0
                  1/4, 1/4, 0];
    SSPRK3BT_B = [1/6, 1/6, 2/3];
end