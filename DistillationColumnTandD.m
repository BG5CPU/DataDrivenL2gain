clear; clc;

td    =  1.5; % sampling time
Ts    =  td;
dimN  =  7;   % dimension of state
dimM  =  3;   % dimension of input




%% system model ===========================================================
Ad = [0,      0,      0,      0,      0,      0,      0;
      1,      0,      0,      0,      0,      0,      0;
      0,      0.0548, 0.9217, 0,      0,      0,      0;
      0.4882, 0,      0,      0.9024, 0,      0,      0;
      0,      0,      0,      0,      0.8607, 0,      0;
      0,      0,      0,      0,      0,      0.8607, 0;
      0,      0,      0,      0,      0,      0,      0.9024 ];

Bd = [1,      0,      0;
      0,      0,      0;
      0,      0,      0;
      0,      0,      0;
      0.0348, 0,      0.2298;
      0,      0.0808, 0;
      0,      0.3564, 0 ];

Cd = eye(dimN);






%% loop for noise energy ==================================================
% settings for algorithm
tolNorm = 1e-6;
tol1  = 0.01;
tol2  = 0.01;

% number of noise bound
Nnoise = 100; 
% energyNoise = logspace(log10(0.000054),log10(0.00054),Nnoise);
energyNoise = linspace(0.0, 0.00034, Nnoise);
% energyNoise = [0, energyNoise];

% record the performance with respect to noise
gamma_noise0 = cell(1, Nnoise); hinfy_noise0 = cell(1, Nnoise); % performance for K0
gamma_noise  = cell(1, Nnoise); hinfy_noise  = cell(1, Nnoise); % performance for final K

% number of data sets
NdataSet = 50; 

% settings for datasets
eu = 2;
Nsample =  100;
trajt   =  1;
dataMatrixLength = Nsample*trajt;
NoiseBoundPrior = calcuNoiseBound(dimN, dimM, dataMatrixLength, energyNoise);
tm1  = dimN;
tm2  = dimN+dimN+dimM;

% start parallel
parfor iq = 1:Nnoise
    energyNoiseNow = energyNoise(iq);
    noisBondPriNow = NoiseBoundPrior(iq);
    % record the performance for each dataset
    gammaRecord0 = zeros(1,NdataSet); HinfyRecord0 = zeros(1,NdataSet);
    gammaRecord  = zeros(1,NdataSet); HinfyRecord  = zeros(1,NdataSet);
    for ik = 1:NdataSet
        ed = energyNoiseNow; % state noise bound
        ew = energyNoiseNow; % input noise bound
        U0 = []; X0 = []; X1 = []; D0 = []; 
        % generate data
        for it = 1:trajt
            uu = (rand(dimM, Nsample)-0.5)*2* eu;
            xx = zeros(dimN, Nsample+1); xx(:,1) = randn(dimN,1) *0.001;
            dd = (rand(dimN, Nsample+1)-0.5)*2* ed;
            dw = (rand(dimM, Nsample)-0.5)*2* ew;    
            for ii = 1:Nsample
                xx(:, ii+1) =  Ad*xx(:, ii) + Bd*uu(:, ii);
            end
            xd = xx+dd; uw = uu+dw; % data with noise
            % creat data matrices
            U0_ = uw; X0_ = xd(:, 1:Nsample); X1_ = xd(:, 2:Nsample+1);
            D0_ = [ dd(:, 2:Nsample+1);
                    dd(:, 1:Nsample);
                    dw(:, 1:Nsample) ];
            U0  = [U0, U0_]; X0  = [X0, X0_]; 
            X1  = [X1, X1_]; D0  = [D0, D0_];
        end
        W0 = [X0; U0];
        % noise bound
        De   = noisBondPriNow * eye(dimN+dimN+dimM);
        De11 = De(1:tm1,     1:tm1); De12 = De(1:tm1,     tm1+1:tm2);  
        De21 = De(tm1+1:tm2, 1:tm1); De22 = De(tm1+1:tm2, tm1+1:tm2); 
        % data consistent set
        Adz  =  W0*W0' - De22;
        Bdz  = -X1*W0' + De12;
        Cdz  =  X1*X1' - De11;
        
        % find the initial controller 
        [isSolved0, vP, vK] = findControllerPetersenMIMO(dimN, dimM, Adz, Bdz, Cdz);
        if isSolved0
            HinfyRecord0(ik) = calcuHinftyNorm(Ad, Bd, Cd, vK, Ts, tolNorm);
        else
            gammaRecord0(ik) = 0; HinfyRecord0(ik) = 0;
            gammaRecord(ik)  = 0; HinfyRecord(ik)  = 0;
            vK = 0;
            disp("warning0")
        end

        if isSolved0
            gam0  = 100000;
            gam2Upper = gam0; gam2Lower = 0;
            gamNow = gam0; gamPre = gam0 + tol1*2;
            
            % run the algorithm, find the performance
            isFirstIteration = true;
            while abs(gamPre-gamNow) > tol1
                vKpre  = vK;
                gamPre = gamNow;        
                % bisection
                while abs(gam2Upper-gam2Lower) > tol2
                    gamNow = (gam2Upper+gam2Lower)/2;
                    [isSolved, ~, ~] = findGamaLamda(dimN, dimM, Adz, Bdz, Cdz, vK, gamNow);
                    if isSolved
                        gam2Upper = gamNow;
                    else
                        gam2Lower = gamNow;
                    end
                end
                disp("gam2Upper = " + gam2Upper);  
                if isFirstIteration
                    gammaRecord0(ik) = gam2Upper;
                    isFirstIteration = false;
                end
                % ADMM
                [isSolvedA, vLmd, vGam] = findGamaLamda(dimN, dimM, Adz, Bdz, Cdz, vK, gam2Upper);
                if ~isSolvedA
                    disp("warningA") 
                    gamNow = gamPre;
                    vK = vKpre;
                end
                if isSolvedA
                    [isSolvedB, vK, gamNow] = findKgama2(dimN, dimM, Adz, Bdz, Cdz, vLmd, vGam);
                    if ~isSolvedB
                        disp("warningB")
                        gamNow = gamPre;
                        vK = vKpre;
                    end
                    disp("gamNow = " + gamNow);
                    gam2Upper = gamNow;
                    gam2Lower = 0;
                end
            end % end while abs(gamPre-gamNow) > tol1
            gammaRecord(ik) = gamNow;
            HinfyRecord(ik) = calcuHinftyNorm(Ad, Bd, Cd, vK, Ts, tolNorm);     
        end % end if isSolved0             
    end % end for ik = 1:NdataSet
    gamma_noise0{iq} = gammaRecord0; hinfy_noise0{iq} = HinfyRecord0;
    gamma_noise{iq}  = gammaRecord;  hinfy_noise{iq}  = HinfyRecord;
end









%% data pre_process =======================================================
% load my test data to save time
% clear; clc;
% load("data\DistillationTandD1.mat");

gamma_ini = cell(1, Nnoise); hinfy_ini = cell(1, Nnoise);
gamma_fnl = cell(1, Nnoise); hinfy_fnl = cell(1, Nnoise);
badValue_hinfy = 0;
badValue_gama1 = 0; 
badValue_gama2 = 100000;
noisBound_axis = [];
gamma_ini_max  = []; gamma_ini_min  = [];
gamma_ini_mean = []; gamma_ini_medi = [];
hinfy_ini_max  = []; hinfy_ini_min  = [];
hinfy_ini_mean = []; hinfy_ini_medi = [];
gamma_fnl_max  = []; gamma_fnl_min  = [];
gamma_fnl_mean = []; gamma_fnl_medi = [];
hinfy_fnl_max  = []; hinfy_fnl_min  = [];
hinfy_fnl_mean = []; hinfy_fnl_medi = [];
for ip = 1:Nnoise
    index1 = gamma_noise0{ip} ~= badValue_gama1;
    index2 = gamma_noise0{ip} ~= badValue_gama2;
    index3 = gamma_noise{ip} ~= badValue_gama1;
    index4 = gamma_noise{ip} ~= badValue_gama2;
    index5 = hinfy_noise0{ip} ~= badValue_hinfy;
    index6 = hinfy_noise{ip} ~= badValue_hinfy;
    indexA = index1 & index2 & index3 & index4 & index5 & index6;

    gamma_ini{ip} = gamma_noise0{ip}(indexA);
    hinfy_ini{ip} = hinfy_noise0{ip}(indexA);
    gamma_fnl{ip} = gamma_noise{ip}(indexA);
    hinfy_fnl{ip} = hinfy_noise{ip}(indexA);

    if any(indexA)
        noisBound_axis = [noisBound_axis, NoiseBoundPrior(ip)];
        gamma_ini_max  = [gamma_ini_max, max(gamma_ini{ip})];
        gamma_ini_min  = [gamma_ini_min, min(gamma_ini{ip})];
        gamma_ini_mean = [gamma_ini_mean, mean(gamma_ini{ip})];
        gamma_ini_medi = [gamma_ini_medi, median(gamma_ini{ip})];

        hinfy_ini_max  = [hinfy_ini_max, max(hinfy_ini{ip})];
        hinfy_ini_min  = [hinfy_ini_min, min(hinfy_ini{ip})];
        hinfy_ini_mean = [hinfy_ini_mean, mean(hinfy_ini{ip})];
        hinfy_ini_medi = [hinfy_ini_medi, median(hinfy_ini{ip})];

        gamma_fnl_max  = [gamma_fnl_max, max(gamma_fnl{ip})];
        gamma_fnl_min  = [gamma_fnl_min, min(gamma_fnl{ip})];
        gamma_fnl_mean = [gamma_fnl_mean, mean(gamma_fnl{ip})];
        gamma_fnl_medi = [gamma_fnl_medi, median(gamma_fnl{ip})];

        hinfy_fnl_max  = [hinfy_fnl_max, max(hinfy_fnl{ip})];
        hinfy_fnl_min  = [hinfy_fnl_min, min(hinfy_fnl{ip})];
        hinfy_fnl_mean = [hinfy_fnl_mean, mean(hinfy_fnl{ip})];
        hinfy_fnl_medi = [hinfy_fnl_medi, median(hinfy_fnl{ip})];
    end
end





%% plot the figures =======================================================
close all;

plotDataL = length(noisBound_axis)-25;
figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 500, 150]);
set(gca,'position',[0.05 0.33 0.93 0.65]);
plot(noisBound_axis(1:plotDataL), sqrt(gamma_ini_mean(1:plotDataL)),...
     'b--', 'LineWidth', 1); 
hold on;
plot(noisBound_axis(1:plotDataL), sqrt(gamma_fnl_mean(1:plotDataL)),...
     'r-', 'LineWidth', 1); 
grid on;
set(gca,'fontsize',12,'fontname','Times');
axis([0 1e-4 0 75]);
xlabel('noise energy bound $\theta$', 'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',12);
legend('$\gamma_{13}^{(1)}$',... 
       '$\gamma$, output of Algorithm', ...
       'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',11);




figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 500, 150]);
set(gca,'position',[0.05 0.33 0.93 0.65]);
plot(noisBound_axis(1:plotDataL), hinfy_ini_mean(1:plotDataL),...
     'b--', 'LineWidth', 1); 
hold on;
plot(noisBound_axis(1:plotDataL), hinfy_fnl_mean(1:plotDataL),...
     'r-', 'LineWidth', 1); 
grid on;
set(gca,'fontsize',12,'fontname','Times');
axis([0 1e-4 0 50]);
xlabel('noise energy bound $\theta$', 'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',12);
legend('$\mathcal{H}_{\infty}$ norm of $K_0$',... 
       '$\mathcal{H}_{\infty}$ norm of $K$', ...
       'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',11);















%% useful function ========================================================
function Ninf = calcuHinftyNorm(A, B, C, K, st, tol)
    Aclo = A+B*K;
    Bclo = B*[K, eye(size(B,2))];
    Cclo = C;
    Dclo = zeros(size(C,1),size(Bclo,2));
    Gclo = ss(Aclo,Bclo,Cclo,Dclo,st);
    [ninf,~] = norm(Gclo,Inf,tol);
    Ninf = ninf;
end

function noiseb = calcuNoiseBound(n, m, L, bound)
    noiseb = zeros(size(bound));
    for ii = 1:length(bound)
        Dd   =  ones(n+n, L)*bound(ii);
        Dw   =  ones(m, L)*bound(ii);
        DAll =  [Dd;Dw];
        noiseb(ii) = max(svd(DAll*DAll'));
    end
end














