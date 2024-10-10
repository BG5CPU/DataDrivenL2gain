function [isSolved, vK, vgam2] = findKgama2(dimN, dimM, Adz, Bdz, Cdz, lmda, Gama)

tune = 1e-10;

% dimAll = dimN + dimN + dimN + dimN+dimM + dimN+dimM;


zerosN_NM = zeros(dimN,dimN+dimM);
zerosNM_N = zeros(dimN+dimM,dimN);

cvx_begin sdp quiet
    variable K(dimM, dimN);
    variable gam2;

    K1 = [eye(dimN); K];
    K2 = [zeros(dimN), zeros(dimN,dimM); K, eye(dimM)];

    blockS = [ -Gama-Cdz,      zeros(dimN),     zeros(dimN),  zerosN_NM,                  Bdz;
                zeros(dimN),  -lmda*eye(dimN),  Gama,         zerosN_NM,                  zerosN_NM;
                zeros(dimN),   Gama,           -Gama,         zerosN_NM,                 -Gama*K1'; 
                zerosNM_N,     zerosNM_N,       zerosNM_N,   -lmda*gam2*eye(dimN+dimM),  -lmda*K2';
                Bdz',          zerosNM_N,      -K1*Gama,     -lmda*K2,                   -Adz ];
    
    minimize(gam2);

    subject to   
        blockS <= 0;
        gam2 >= tune;
cvx_end


if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    isSolved = true;
    vK = K; 
    vgam2  = gam2;
else
    isSolved = false;
    vK = 0; 
    vgam2  = 0;
end

end