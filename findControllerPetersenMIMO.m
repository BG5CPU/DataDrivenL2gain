function [isSolved, vP, vK] = findControllerPetersenMIMO(dimN, dimM, Adz, Bdz, Cdz)

tune = 1e-8;

dimAll = dimN + dimN + dimN + dimM;

cvx_begin sdp quiet
    variable P(dimN, dimN) semidefinite;
    variable Y(dimM, dimN);

    blockS = [ -P-Cdz,        zeros(dimN),  Bdz;
                zeros(dimN), -P,           [P, Y'];
                Bdz',        [P; Y],       -Adz ];
        
    % minimize( -trace(P) );
    
    subject to   
        blockS <= -tune*eye(dimAll);
cvx_end

if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    isSolved = true;
    vP = P; 
    vK = Y/P;
else
    isSolved = false;
    vP = 0;
    vK = 0;
end

end