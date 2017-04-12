function varargout = getCobraComParams(param2get, options, modelCom)
%get the required default parameters
%[param_1, ..., param_N] = getCobraComParams({'param_1',...,'param_N'},options,modelCom)
%Input:
%   'param_1',...,'param_N': parameter names
%   options: option structure
%            If the required parameter is a field in options, take from
%            options. Otherwise, return the default value.
%   modelCom: the community model for which parameters are constructed.
   
if nargin < 3
    modelCom = struct('rxns',[]);
    modelCom.infoCom.spAbbr = {};
    modelCom.infoCom.rxnSps = {};
end
if nargin < 2 || isempty(options)
    options = struct();
end
if ischar(param2get)
    param2get = {param2get};
end
paramNeedTransform = {'GRfx', 'BMlb', 'BMub', 'BMfx'};
    
varargout = cell(numel(param2get), 1);
for j = 1:numel(param2get)
    if any(strcmp(param2get{j}, paramNeedTransform))
        %if need transformation
        varargout{j} = transformOptionInput(options, param2get{j}, numel(modelCom.infoCom.spAbbr));
    elseif isfield(options, param2get{j})
        %if provided in the call
        varargout{j} = options.(param2get{j});
    else
        %use default if default exist and not provided
        %return empty if no default
        varargout{j} = paramDefault(param2get{j}, modelCom);
    end
    %if calling for a directory, make sure to return a new directory
    if strcmp(param2get{j}, 'directory')
        k = 0;
        while exist(varargout{j}, 'file')
            k = k + 1;
            varargout{j} = [paramDefault.directory num2str(k)];
        end
    end
end


end

function param = paramDefault(paramName,modelCom)

switch paramName
    %general parameters
    case 'threads',     param = 1; %threads for general computation, 0 or -1 to turn on maximum no. of threads
    case 'verbFlag',	param = 3;%verbal dispaly
    case 'loadModel',   param = '';
    case 'CplexParam',  %default Cplex parameter structure
        [param.simplex.display, param.tune.display, param.barrier.display,...
            param.sifting.display, param.conflict.display] = deal(0);
        [param.simplex.tolerances.optimality, param.simplex.tolerances.feasibility] = deal(1e-9,1e-8);
        param.read.scale = -1;
    %parameters for createCommModel
    case 'metExId',     param = '[e]';
    %parameters for optimizeCbModelCom
    case 'GRguess',     param = 0.2;%initial guess for growth rate
    case 'BMtol',       param = 0.8;%tolerance for relative biomass amount (used only for feasCrit=3)
    case 'BMtolAbs',    param = 1e-5;%tolerance for absolute biomass amount
    case 'GR0',         param = 0.001;%small growth rate to test growth
    case 'GRtol',       param = 1e-5;%gap for growth rate convergence
    case 'GRdev',       param = 1e-5;%percentage deviation from the community steady state allowed
    case 'maxIter',     param = 1e3;%maximum no. of iterations
    case 'feasCrit',    param = 1; %feasibility critierion
    case 'algorithm',   param = 1;%1:invoke Fzero after getting bounds; 2:simple guessing algorithm
    case 'BMgdw',       param = ones(numel(modelCom.infoCom.spAbbr), 1);%relative molecular weight of biomass. For scaling the relative abundance
    case 'saveModel',   param = '';
    case 'BMobj',       param = ones(numel(modelCom.infoCom.spBm),1); %objective coefficient for each species
    case 'BMweight',    param = 1; %sum of biomass for feasibility criterion 1
    case 'LPonly',      param = false;%true to return LP only but not calculate anything
    case 'solveGR0',    param = false;%true to solve the model at very low growth rate (GR0)
    case 'resultTmp',   param = struct('GRmax',[],'vBM',[],'BM',[],'Ut',[],...
                                'Ex',[],'flux',[],'iter0',[],'iter',[],'stat','');%result template
    %parameters for optimizeCbModelComCplexFixGr
    case 'GR',          param = 0.2; %given growth rate
    %parameters for fluxVarCom
    case 'optBMpercent',param = 99.99;
    case 'rxnNameList', if isfield(modelCom.infoCom,'spBm'), param = modelCom.rxns(findRxnIDs(modelCom,modelCom.infoCom.spBm));else param = modelCom.rxns;end
    case 'rxnFluxList', if isfield(modelCom.infoCom,'spBm'), param = modelCom.rxns(findRxnIDs(modelCom,modelCom.infoCom.spBm));else param = modelCom.rxns;end
    case 'BMmaxLB',     param = 1; %maximum biomass when it is unknown
    case 'BMmaxUB',     param = 1; %maximum biomass when it is unknown
    case 'optGRpercent',param = 99.99;
    case 'NstepGR',     param = 1;
    case 'saveFVA',     param = '';
    case 'removeCycles',param = true;
    %parameters for fluxCoupleCom
    case 'Nstep',       param = 10;
    case 'NstepScale',  param = 'lin';
    case 'tol',         param = 1e-8;
    case 'symmetric',   param = true; %treat it as symmetric, optimize for only j > k
    case 'saveTmp',     param = false;
    case 'savePOA',     param = 'POAtmp/POA';
    %parameters for fluxSampleCom
    case 'maxTime',     param = 20 * 60;
    case 'maxTime0',    param = 60 * 60 * 2;
    case 'maxStep',     param = 250;
    case 'saveSampling',param = 'sampling_temp';
    case 'saveFre',     param = 100;
    case 'warmupPtsLoad', param = false;
    case 'sampleThreads', param = 1; %threads for sampling, -1 for using matlab parallel toolbox
    case 'cont',        param = false;
    case 'Nmix',        param = 3;
    case 'mixFrac',     param = 0.7;
    case 'bias',        param = struct('method','uniform',...
                                'index',((numel(modelCom.lb) + 1) : ...
                                (numel(modelCom.lb) + numel(modelCom.infoCom.spAbbr)))');
    case 'FVAbound',    param = true;
    % for MMCbCom
    case 'xLB',         param = 0.001;
    case 'atpm',        param = 0.1;
    % for dSteadyCom
    case 'X0', %if not given initial biomass, get a random small one
        param = rand(numel(modelCom.infoCom.spAbbr), 1);
        param = param * 0.001 / sum(param);
    case 'C0',          param = zeros(size(modelCom.infoCom.EXcom,1),1); %concentration vector for community metabolites
    case 'dt',          param = 1/60; %time step (hour)
    %case 'T',           param = 10; %maximum simulation time (hour)
    case 'concTol',     param = 0; %tolerance for concentration below which assume no uptake is possible
    case 'saveDsteady', param = 'dSteadyTmp'; %save name for dSteady
    case 'saveFlux',    param = false;
    % for SteadyComOxy
    case 'X_muc',       param = 0.1;
    case 'X_lum',       param = 0.1;
    case 'O_muc',       param = 0.5;
    case 'O_lum',       param = 0.1;
    case 'C_ex',        param = 100;
    case 'mu_muc',      param = 0.02;
    case 'O2',          param = 'o2[u]';
    case 'H',           param = 'h[u]';
    case 'L',           param = 10; %number of discretization step for calculating oxygen diffusion
    % ESScc
    case 'u',           param = modelCom.ub(findRxnIDs(modelCom,modelCom.infoCom.EXcom(:,1)));%minimum % difference of r and v_neg
    case 'eps1',        param = 1e-3;%minimum % difference of r and v_neg
    case 'M',           param = 1000;%|ub| and |lb| of modelCom >= M, treat the reaction as unbounded
    case 'indicator',   param = true;
    case 'spBmId',      [~,param] = ismember(findRxnIDs(modelCom,modelCom.infoCom.spBm),find(~strcmp(modelCom.infoCom.rxnSps,'com')));
    case 'r',           param = ones(size(modelCom.infoCom.EXcom,1),numel(modelCom.infoCom.spAbbr));
    case 'kcat',        param = 10*ones(size(modelCom.infoCom.EXcom,1),numel(modelCom.infoCom.spAbbr));
    case 'Km',          param = ones(size(modelCom.infoCom.EXcom,1),numel(modelCom.infoCom.spAbbr));
    case 'Rc',          param = ones(numel(modelCom.infoCom.spAbbr),1);
    case 'Rm',          param = zeros(numel(modelCom.infoCom.spAbbr),1);
    %SteadyGut
    case 'Conc0',       param = 1000*ones(size(modelCom.infoCom.EXcom,1),1); %initial concentration (diet)
    case 'Xm',          param = [6.75e-6,2.7e-6,1.5e-5,1.95e-4,5.25e-4,7.5e-6,9e-6]'; %mucosal biomass density
    case 'volMuc',      param = ones(7,1);
    case 'volLum',      param = 10*ones(7,1);
    case 'T',           param = [2 2 2 3 3 3 3]';
    case 'O2sum',       param = (2666+2/3) / 6.02214086e23 * 1000 / 1e-13 * 60 * 60 * 1e-5;
    case 'O2muc',       param = [0.8 0.85 0.85 0.99 0.95 0.95 1]';
    case 'O2spRateMuc', param = 2*ones(7, numel(modelCom.infoCom.spAbbr));
    case 'O2spRateLum', param = 2*ones(7, numel(modelCom.infoCom.spAbbr));
    case 'deathRate',   param = 0; %deathRate of biomass not able to satisfy ATPM
    case 'fluxFlag',    param = false;
    case 'steadyTol',   param = 1e-4;
    case 'iterLimit',   param = 1e3;
    case 'epsC0',       param = 1e-7;
    case 'M_unlim',     param = 1000;
    case 'sf',          param = 1e4; %scaling factor
    case 'multiSS',     param = 5; %number of local minimum to find
    %optimizeCbModel's original parameters
    case 'osenseStr',   param = 'max';
    case 'minNorm',     param = 0; % 0 usual LP
                                    % 1 minimal sum of flux
                                    % >1 or n x 1 vector, 
                                    % minimize Euclidean norm with uniform weight or weights in the vector
                                    % (to be implemented)
    case 'allowLoops',  param = true;
    otherwise,          param = [];
end
end

function x = transformOptionInput(options, field, nSp)
if isfield(options, field)
    if size(options.(field), 2) == 2
        x = NaN(nSp, 1);
        x(options.(field)(:,1)) = options.(field)(:,2);
    else
        x = options.(field);
    end
else
    x = NaN(nSp, 1);
end

end
