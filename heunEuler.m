function Heun = HeunEuler()
    Heun.C = [0,1];
    Heun.B = [1/2,1/2;1,0];
    Heun.A = [0,0;1,0];
end