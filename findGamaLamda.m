function [isSolved, vLamda, vGama] = findGamaLamda(dimN, dimM, Adz, Bdz, Cdz, K, gam2)

tune = 1e-5;

% dimAll = dimN + dimN + dimN + dimN+dimM + dimN+dimM;

K1 = [eye(dimN); K];
K2 = [zeros(dimN), zeros(dimN,dimM); K, eye(dimM)];
zerosN_NM = zeros(dimN,dimN+dimM);
zerosNM_N = zeros(dimN+dimM,dimN);

cvx_begin sdp quiet
    variable Gama(dimN, dimN) semidefinite;
    variable lmda;

    blockS = [ -Gama-Cdz,      zeros(dimN),     zeros(dimN),  zerosN_NM,                  Bdz;
                zeros(dimN),  -lmda*eye(dimN),  Gama,         zerosN_NM,                  zerosN_NM;
                zeros(dimN),   Gama,           -Gama,         zerosN_NM,                 -Gama*K1'; 
                zerosNM_N,     zerosNM_N,       zerosNM_N,   -lmda*gam2*eye(dimN+dimM),  -lmda*K2';
                Bdz',          zerosNM_N,      -K1*Gama,     -lmda*K2,                   -Adz ];

%     minimize(-lmda);

    subject to   
        blockS <= 0; %-tune*eye(dimAll);
        Gama >= tune*eye(dimN);
        lmda >= tune;
cvx_end

if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    isSolved = true;
    vLamda = lmda; 
    vGama  = Gama;
else
    isSolved = false;
    vLamda = 0; 
    vGama  = 0;
end

end