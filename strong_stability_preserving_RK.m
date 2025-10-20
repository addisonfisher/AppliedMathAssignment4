function SSPRK3BT = strong_stability_preserving_RK();
    SSPRK3BT = struct();

    SSPRK3BT.C = [0,   1,   1/2];
    SSPRK3BT.A = [0,   0,   0
                1,   0,   0
                1/4, 1/4, 0];
    SSPRK3BT.B = [1/6, 1/6, 2/3];
end