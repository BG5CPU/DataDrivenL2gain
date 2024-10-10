clear; clc;

dimN  =  7;   % dimension of state
dimM  =  3;   % dimension of input


Nnoise = 5; % number of noise bound
energyNoise = logspace(-3,-2.2,Nnoise);

NdataSet = 10; % number of data sets

RecordNoiseBound = cell(1,Nnoise);
RecordBoundUpper = cell(1,Nnoise);


for iq = 1:Nnoise
    energyNoiseNow = energyNoise(iq);
    recordNoiseBound = zeros(1,NdataSet);
    recordBoundUpper = zeros(1,NdataSet);
    for ik = 1:NdataSet
        Nsample =  100;
        trajt   =  1;
        eu = 2;
        ed = energyNoiseNow; % state noise bound
        ew = energyNoiseNow; % input noise bound
        D0 = [];
        % generate data
        for it = 1:trajt
            dd = (rand(dimN, Nsample+1)-0.5)*2* ed;
            dw = (rand(dimM, Nsample)-0.5)*2* ew;    
            D0_ = [ dd(:, 2:Nsample+1);
                    dd(:, 1:Nsample);
                    dw(:, 1:Nsample) ];
            D0  = [D0, D0_];
        end
        Ndat =  size(D0,2)
        Dd   =  ones(dimN+dimN, Ndat)*ed;
        Dw   =  ones(dimM, Ndat)*ew;
        DAll =  [Dd;Dw];
        det  =  max(svd(D0*D0'));
        de   =  round_up_to_two_significant_figures(det);
        % disp("true noise bound is " + det);
        % disp("upper bound of boud is " + de);
        recordNoiseBound(ik) = det;
        recordBoundUpper(ik) = de;
    end
    RecordNoiseBound{iq} = recordNoiseBound;
    RecordBoundUpper{iq} = recordBoundUpper;
end






% useful functions ========================================================
function rounded_value = round_up_to_two_significant_figures(value)
    if value == 0
        rounded_value = 0;
    else
        magnitude = floor(log10(abs(value)));
        scaled_value = value * 10^(-magnitude + 1);
        rounded_scaled_value = ceil(scaled_value);
        rounded_value = rounded_scaled_value * 10^(magnitude - 1);
    end
end
























