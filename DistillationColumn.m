clear; clc;

td    =  1.5; % sampling time
Ts    =  td;
dimN  =  7;    % dimension of state
dimM  =  3;    % dimension of input


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




%% generate data ==========================================================
Nsample =  100;
trajt   =  1;

eu = 2;
% this is for noise free case
ed = 0.00143 *0;
ew = 0.00143 *0;

% this is for noisy case
% ed = 0.00143 *1;
% ew = 0.00143 *1;


U0 = []; X0 = [];
X1 = []; D0 = [];

for it = 1:trajt
    uu = (rand(dimM, Nsample)-0.5)*2* eu;
    xx = zeros(dimN, Nsample+1); xx(:,1) = randn(dimN,1) *0.001;
    dd = (rand(dimN, Nsample+1)-0.5)*2* ed;
    dw = (rand(dimM, Nsample)-0.5)*2* ew;    
    for ii = 1:Nsample
        xx(:, ii+1) =  Ad*xx(:, ii) + Bd*uu(:, ii);
    end
    xd = xx+dd;
    uw = uu+dw;
    % creat data matrices
    U0_ = uw;
    X0_ = xd(:, 1:Nsample);
    X1_ = xd(:, 2:Nsample+1);
    D0_ = [ dd(:, 2:Nsample+1);
            dd(:, 1:Nsample);
            dw(:, 1:Nsample) ];
    U0  = [U0, U0_];
    X0  = [X0, X0_];
    X1  = [X1, X1_];
    D0  = [D0, D0_];
end


W0 = [X0; U0];

Tesd = [Ad Bd]*W0 + [eye(dimN), -Ad, -Bd]*D0 - X1;

Ndat =  size(U0,2);
Dd   =  ones(dimN+dimN, Ndat)*ed;
Dw   =  ones(dimM, Ndat)*ew;
DAll =  [Dd;Dw];
% De   =  max(svd(DAll*DAll')) * eye(dimN+dimN+dimM);
de   =  ceil( max(svd(D0*D0'))* 1e5 )* 1e-5;
% de   =  10e-5;
% de   = 0;
De   =  de * eye(dimN+dimN+dimM);




tm1  = dimN;
tm2  = dimN+dimN+dimM;
De11 = De(1:tm1,     1:tm1); De12 = De(1:tm1,     tm1+1:tm2);  
De21 = De(tm1+1:tm2, 1:tm1); De22 = De(tm1+1:tm2, tm1+1:tm2); 



Adz  =  W0*W0' - De22;
Bdz  = -X1*W0' + De12;
Cdz  =  X1*X1' - De11;

eig_X0  = eig(W0*W0');
eig_De  = eig(De22);
eig_Adz = eig(Adz);






%% find the controller ====================================================
[vP, vK] = findControllerPetersenMIMO(dimN, dimM, Adz, Bdz, Cdz);
disp(vK);
eig_vK = eig(Ad+Bd*vK);
disp(eig_vK);
vK0 = vK;

disp( "min_eig_Adz = " + min(eig_Adz) );
disp( "max_eig_Adz = " + max(eig_Adz) );
disp( "noise energy = " + max(svd(D0*D0')) );
disp( "energy bound = " + de );



%% ADMM test ==============================================================
tolNorm = 1e-6;
tol1  = 0.01;
tol2  = 0.01;
gam0  = 100000;

gam2Upper = gam0;
gam2Lower = 0;

gamNow = gam0;
gamPre = gam0 + tol1*2;

gamRecord = [];
vKRecord  = [];
InfRecord = [];

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
    gamRecord = [gamRecord; [0, gam2Upper]];
    vKRecord  = [vKRecord; [zeros(size(vK,1),1), vK]];
    InfRecord = [InfRecord; [0, calcuHinftyNorm(Ad, Bd, Cd, vK, Ts, tolNorm)]];
    disp("gam2Upper = " + gam2Upper);
    
    % ADMM
    [isSolvedA, vLmd, vGam] = findGamaLamda(dimN, dimM, Adz, Bdz, Cdz, vK, gam2Upper);
    [isSolvedB, vK, gamNow] = findKgama2(dimN, dimM, Adz, Bdz, Cdz, vLmd, vGam);
    if ~isSolvedB
        disp("warning")
        gamNow = gamPre;
        vK = vKpre;
    end
    disp("gamNow = " + gamNow);
    gam2Upper = gamNow;
    gam2Lower = 0;
    gamRecord = [gamRecord; [1, gam2Upper]];
    vKRecord  = [vKRecord; [ones(size(vK,1),1), vK]];
    InfRecord = [InfRecord; [1, calcuHinftyNorm(Ad, Bd, Cd, vK, Ts, tolNorm)]];
end


%% load my test data to save time =========================================
% clear; clc;
% load("data\DistillationNoise1.mat");
% load("data\DistillationNoisefree1.mat");



%% plot figures ===========================================================
gamPlot = [];
infPlot = [];
for ii = 1 : size(gamRecord,1)
    if( gamRecord(ii,1) == 1 )
        gamPlot = [gamPlot; gamRecord(ii,2)];
        infPlot = [infPlot; InfRecord(ii,2)];
    end
end

gamPlot = sqrt(gamPlot);

ite = ( 1:size(gamPlot,1) )';

figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 500, 200]);
subplot(2,1,1);
plot(ite, gamPlot, '*', 'MarkerSize', 8, 'Color', 'k'); grid on;
legend('$\gamma$', 'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',11);
set(gca,'fontsize',12,'fontname','Times');
set(gca,'position',[0.072 0.68 0.90 0.30]);
axis([ite(1) length(gamPlot)-1 min(gamPlot) max(gamPlot)]);
subplot(2,1,2);
plot(ite, infPlot, 'o', 'MarkerSize', 8, 'Color', 'r'); grid on;
set(gca,'fontsize',12,'fontname','Times');
xlabel('iterations', 'Fontname', 'Times New Roman','FontSize',12);
legend('$\mathcal{H}_{\infty}$ norm', 'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',11);
set(gca,'position',[0.072 0.22 0.90 0.30]);
axis([ite(1) length(infPlot)-1 min(infPlot) max(infPlot)]);


figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 500, 120]);
plot(ite, gamPlot, '*', 'MarkerSize', 8, 'Color', 'k'); grid on;
set(gca,'fontsize',12,'fontname','Times');
set(gca,'position',[0.05 0.35 0.93 0.60]);
legend('$\gamma^{(i)}_{15}$',... 
       'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',11);
axis([ite(1) length(gamPlot)-1 min(gamPlot) max(gamPlot)]);
axis([ite(1) length(gamPlot)-1 min(gamPlot) max(gamPlot)]);
xlabel('$i$th-iteration', ...
       'Fontname', 'Times New Roman','Interpreter','latex','FontSize',12);





%% test the controller Implementation =====================================
euIm = 0.001; % input noise
exIm = 0.001; % output noise

N = 1000;
t = 0:td:(N-1)*td;
tx = repmat(t, dimN, 1);
tu = repmat(t, dimM, 1);

% fre = 2*pi*0.5; % rad/s
fre = InfRecord(end,3)/td;

Unoise = (rand(dimM, N)-0.5)*2* euIm*0 +... 
         ones(dimM, N)* euIm*0 +... 
         sin(fre*tu)* euIm *1 +... 
         (randi([0, 1], dimM, N)-0.5)*2* euIm*0;

Xnoise = (rand(dimN, N)-0.5)*2* exIm*0 +... 
         ones(dimN, N)* exIm*0 +...
         sin(fre*tx)* exIm *1 +...
         (randi([0, 1], dimN, N)-0.5)*2* exIm*0;

deltXU = [ Xnoise; -Unoise ];
NoiSum = zeros(1,size(deltXU,2));
NoiTem = NoiSum;
for ik = 1:size(deltXU,2)
    NoiTem(ik) = deltXU(:,ik)'*deltXU(:,ik);
    NoiSum(ik) = sum(NoiTem(1:ik));
end


ic = randn(dimN,1)*0; % initial condition


UactualOD = zeros(dimM, N); % actual input
UmeasurOD = zeros(dimM, N); % measured input

XactualOD = zeros(dimN, N); % actual output
XmeasurOD = zeros(dimN, N); % measured output
XsumOD    = zeros(1, N);

xo = ic; % actual system state with randn initial value

% implementation of the ordinary controller
% closed-loop system
for k = 1:N
    XactualOD(:,k) = xo;
    XmeasurOD(:,k) = xo + Xnoise(:,k);

    XsumOD(k) = trace(XactualOD(:,1:k)'*XactualOD(:,1:k));

    UmeasurOD(:,k) =  vK0 * XmeasurOD(:,k);
    UactualOD(:,k) =  UmeasurOD(:,k) - Unoise(:,k);

    xo = Ad*xo+Bd*UactualOD(:,k);
end

% ========================================

UactualOP = zeros(dimM, N); % actual input
UmeasurOP = zeros(dimM, N); % measured input

XactualOP = zeros(dimN, N); % actual output
XmeasurOP = zeros(dimN, N); % measured output
XsumOP    = zeros(1, N);

xo = ic; % actual system state with randn initial value

% implementation of the optimal controller
% closed-loop system
for k = 1:N
    XactualOP(:,k) = xo;
    XmeasurOP(:,k) = XactualOP(:,k) + Xnoise(:,k);

    XsumOP(k) = trace(XactualOP(:,1:k)'*XactualOP(:,1:k));

    UmeasurOP(:,k) =  vK * XmeasurOP(:,k);
    UactualOP(:,k) =  UmeasurOP(:,k) - Unoise(:,k);

    xo = Ad*xo+Bd*UactualOP(:,k);
end

figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 500, 200]);
% plot(NoiSum, 'k'); hold on;
plot(XsumOD(1:N), 'b--', 'LineWidth', 1.5); hold on;
plot(XsumOP(1:N), 'r', 'LineWidth', 1.5); grid on;
set(gca,'fontsize',12,'fontname','Times');
xlabel('horizon $N$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize',12);
ylabel('$\sum_{k=0}^{N} |x(k)|^2$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize',12);
legend('$K_0$ returned from (19)', ...
       '$K$ returned from Algorithm', ...
       'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',11);
set(gca,'position',[0.12 0.22 0.83 0.74]);












%% useful function ========================================================
function NinfandFpeak = calcuHinftyNorm(A, B, C, K, st, tol)
    Aclo = A+B*K;
    Bclo = B*[K, eye(size(B,2))];
    Cclo = C;
    Dclo = zeros(size(C,1),size(Bclo,2));
    Gclo = ss(Aclo,Bclo,Cclo,Dclo,st);
    [ninf,fpeak] = norm(Gclo,Inf,tol);
    NinfandFpeak = [ninf,fpeak];
end














