function [POAtable, fluxRange, Stat, GRvector] = SteadyComPOACplex(modelCom,options,solverParam)
%Pairwise POA for community model at community steady-state for a range of growth rates
%[POAtable, fluxRange, Stat, GRvector] = SteadyComPOACplex(modelCom,options,solverParam)
%
%INPUT
% modelCom       A community COBRA model structure with the following extra fields:
% (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
% (at least one of the below two is needed)
%   infoCom      structure containing community reaction info 
%                (returned along with the community model created with createCommModel)
%   indCom       the index structure corresponding to infoCom
%
% options (optional) option structure with the following fields:
%   GRmax           Maximum growth rate of the model (default to be found
%                     SteadyComCplex.m)
%   optGRpercent    A vector of percentages. Perform POA at these percents
%                     of max. growth rate respectively (Default = 99.99)
%   optBMpercent    Only consider solutions that yield at least a certain
%                     percentage of the optimal biomass (Default = 95)
%   rxnNameList     List of reactions (IDs or .rxns) to be analyzed.
%                     Use a (N_rxns + N_organism) x K matrix for POA of K
%                     linear combinations of fluxes and/or abundances
%                     (Default = biomass reaction of each species)
%   pairList        Pairs in rxnNameList to be analyzed. N_pair by 2 array of:
%                     - indices referring to the rxns in rxnNameList, e.g.,
%                       [1 2] to analyze rxnNameList{1} vs rxnNameList{2}
%                     - rxn names which are members of rxnNameList, e.g.,
%                       {'EX_glc-D(e)','EX_ac(e)'}
%                     If not supplied, analyze all K(K-1) pairs from the K
%                     targets in rxnNameList.
%   symmetric       Used only when pairList is not supplied. Avoid running 
%                     symmetric pairs (e.g. j vs k and k vs j)
%   Nstep           Number of steps for fixing one flux at a value between 
%                     the min. and the max. possible fluxes. Default 10.
%                   Can also be a vector indicating the fraction of intermediate value to be analyzed
%                   e.g. [0 0.5 1] means computing at minFlux, 0.5(minFlux + maxFlux) and maxFlux
%   NstepScale      Used only when Nstep is a single number. 
%                     -'lin' for a linear (uniform) scale of step size 
%                     -'log' for a log scaling of the step sizes
%   fluxRange       Flux range for each entry in rxnNameList. K x 2 matrix.
%                     Defaulted to be found by SteadyComFVACplex.m
%  (other parameters)
%   savePOA         Must be non-empty. The filename to save the POA results
%                     (default 'POAtmp/POA')
%   threads         > 1 for explicitly stating the no. of threads used,
%                     0 or -1 for using all available threads. Default 1.
%   verbFlag        Verbose output. 0 or 1.
%   loadModel       String of filename to be loaded. If non-empty, load the 
%                     cplex model ('loadModel.mps'), basis ('loadModel.bas') 
%                     and parameters ('loadModel.prm').
%  May add also other parameters in SteadyComCplex for calculating the maximum growth rate.
%
%OUTPUT
% (N_gr = numel(optGRpercent))
% POAtable          K x K cells.  
%                     Each (i,i)-cell contains a Nstep x 1 x N_gr matrix of
%                     the fluxes at which rxnNameList{i} is fixed
%                     (i,j)-cell contains a Nstep x 2 x N_gr matrix,  
%                     (p,:,q)-entry being the range of rxnNameList{j} 
%                     when rxnNameList{i} is fixed at POAtable{i,i}(p,1,q)
%                     at growth rate = GRvector(q)
% fluxRange         K x 2 x N_gr matrix of flux range for each entry in rxnNameList 
% Stat              K x K x N_gr structure array with fields:
%                     -'cor': the slope from linear regression between the
%                             fluxes of a pair
%                     -'r2':  the corresponding coefficient of determination (R-square)
% GRvector          Vector of growth rates being analyzed

%% Initialization
%check required fields for community model
if ~isfield(modelCom,'indCom')
    if ~isfield(modelCom,'infoCom') || ~isstruct(modelCom.infoCom) || ...
            ~all(isfield(modelCom.infoCom,{'spBm','EXcom','EXsp','spAbbr','rxnSps','metSps'}))
        error('infoCom must be provided for calculating the max. community growth rate.\n');
    end
    %get useful reaction indices
    modelCom.indCom = infoCom2indCom(modelCom);
end

%get paramters
if ~exist('options', 'var')
    options = struct();
end
if ~exist('solverParam', 'var') || isempty(solverParam)
    %default Cplex parameters
    solverParam = getCobraComParams('CplexParam');
end
param2get = {'GRmax', 'optGRpercent', 'Nstep', 'GRfx','BMmaxLB','BMmaxUB',...
             'rxnNameList', 'verbFlag', 'savePOA','loadModel'};
eval(sprintf('[%s] = getCobraComParams(param2get, options, modelCom);', ...
            strjoin(param2get, ',')...
            )...
    );
if isempty(savePOA)
    error('A non-empty file name must be provided to save the POA results.');
end
[feasTol, ~] = getCobraSolverParams('LP',{'feasTol'; 'optTol'}, solverParam);
if isfield(solverParam,'simplex') && isfield(solverParam.simplex, 'tolerances')...
        && isfield(solverParam.simplex.tolerances,'feasibility')
    %override the feasTol in CobraSolverParam if given in solverParam
    feasTol = solverParam.simplex.tolerances.feasibility;
else
    %otherwise use the feasTol in COBRA toolbox
    solverParam.simplex.tolerances.feasibility = feasTol;
end

init = true;
if exist([savePOA '_MasterModel.mat'], 'file')
    load([savePOA '_MasterModel.mat'],'LPstart','LPmodel','GRvector','kDisp','idRow');
    LP = Cplex('POA');
    LP.Model = LPmodel;
    LP.Start = LPstart;
    LP = setCplexParam(LP, solverParam);
    if ~isfield(options, 'BMmaxLB')
        options.BMmaxLB = LP.Model.lhs(idRow);
    end
    if ~isfield(options, 'BMmaxUB')
        options.BMmaxUB = LP.Model.rhs(idRow);
    end
    init = false;
end
if numel(Nstep) > 1
    Nstep = numel(Nstep);
end
if ischar(rxnNameList)
    rxnNameList = {rxnNameList};
end
if iscell(rxnNameList)
    Ncheck = numel(rxnNameList);
else
    Ncheck = size(rxnNameList,2);
end
if init
    addRow = false;
    %get maximum growth rate
    if isempty(GRmax)
        if exist('Cplex.p','file') == 6
            options.minNorm = false;
            [~, result,LP] = SteadyComCplex(modelCom, options,solverParam);
        else
            %need further achitecture for using COBRA solver
            warning('Currently support Cplex only.');
            return
        end
        if strcmp(result.stat,'infeasible')
            %infeasible model
            warning('Model is infeasible.');
            POAtable = cell(Ncheck);
            fluxRange = NaN(Ncheck,2); 
            Stat = repmat(struct('cor',[],'r2',[]),Ncheck,Ncheck);
            GRvector = NaN(numel(optGRpercent), 1);
            return
        end
        GRmax = result.GRmax;
        idRow = size(LP.Model.A,1);
    else
        %If GRmax is given, BMmaxLB and BMmaxUB should be included in options in this case to ensure feasibility
        if ~isempty(loadModel)
            [m, n] = size(modelCom.S);
            nRxnSp = sum(modelCom.indCom.rxnSps > 0); %number of species-specific rxns
            nSp = numel(modelCom.indCom.spBm); %number of species
            % load solution if given and growth rate is known
            LP = Cplex('fluxSampling');
            LP.readModel([loadSol '.mps']);
            LP.readBasis([loadSol '.bas']);
            LP.readParam([loadSol '.prm']);
            fprintf('Load model ''%s'' successfully.\n', loadModel);
            addRow = true;
            if size(LP.Model.A,1) > m + 2*nRxnSp + nSp
                [ynRow,idRow] = ismember(sparse(ones(nSp,1),n+1:n+nSp,ones(nSp,1),1,n+nSp),...
                    LP.Model.A(m+2*nRxnSp+nSp+1:end,1:n+nSp),'rows');
                if ynRow
                    idRow = m + 2*nRxnSp + nSp + idRow;
                end
                addRow = ~ynRow;
            end
        else
            %get LP using SteadyComCplex if only growth rate is given 
            options2 = options;
            options2.LPonly = true;
            [~, ~, LP] = SteadyComCplex(modelCom, options2, solverParam);
            %no constraint on total biomass using LPonly option
            addRow = true;
        end
    end
    if addRow
        %add a row for constraining the sum of biomass if not exist
        LP.addRows(BMmaxLB, ...
            sparse(ones(1, nSp), n + 1: n + nSp, ones(1, nSp), 1, size(LP.Model.A,2)),...
            BMmaxUB, 'UnityBiomass');
        idRow = size(LP.Model.A,1);
    else
        %using BMmaxLB and BMmaxUB stored in the LP if not given in options
        if ~isfield(options,'BMmaxLB') %take from LP if not supplied
            BMmaxLB = LP.Model.lhs(idRow);
        end
        if ~isfield(options,'BMmaxUB') %take from LP if not supplied
            BMmaxUB = LP.Model.rhs(idRow);
        end
        LP.Model.lhs(idRow) = BMmaxLB;
        %not allow the max. biomass to exceed the one at max growth rate,
        %can happen if optBMpercent < 100. May dismiss this constraint or 
        %manually supply BMmaxUB in the options if sum of biomass should be variable
        LP.Model.rhs(idRow) = BMmaxUB;
    end
    %set Cplex parameters
    LP = setCplexParam(LP, solverParam);
    LP.Model.A = updateLPcom(modelCom, GRmax, GRfx, [], LP.Model.A, []);
    LP.Model.sense = 'minimize';
    LP.Model.obj(:) = 0;
    LP.solve();
    %check and adjust for feasibility
    dev = checkSolFeas(LP);
    kBMadjust = 0;
    BMmaxLB = LP.Model.lhs(idRow);
    while (~isfield(LP.Solution, 'x') || dev > feasTol) && kBMadjust < 10
        kBMadjust = kBMadjust + 1;
        %row of biomass constraint should at the end
        LP.Model.lhs(idRow) = BMmaxLB * (1 - feasTol/(11 - kBMadjust));
        LP.solve();
        dev = checkSolFeas(LP);
        if verbFlag
            fprintf('BMmax adjusment: %d\n',kBMadjust);
        end
    end
    if (~isfield(LP.Solution, 'x') || dev > feasTol)
        error('Model not feasible.')
    end
    if ~isfield(options, 'BMmaxLB')
        options.BMmaxLB = LP.Model.lhs(idRow);
    end
    if ~isfield(options, 'BMmaxUB')
        options.BMmaxUB = LP.Model.rhs(idRow);
    end

    GRvector = GRmax * optGRpercent/100;
    if numel(optGRpercent) == 1
        kDisp = 2;
    else
        d = max(GRvector(2:end) - GRvector(1:end-1));
        if d < 1
            kDisp = abs(floor(log10(abs(d))));
        else
            kDisp = 0;
        end
    end
end

nVar = size(LP.Model.A,2);
%handle objective matrix
if isnumeric(rxnNameList)
    if size(rxnNameList,1) >= size(modelCom.S,2) && size(rxnNameList,1) <= nVar
        %it is a matrix of objective vectors
        objList = [sparse(rxnNameList); sparse(nVar - size(rxnNameList,1), size(rxnNameList,2))];
    elseif size(rxnNameList,1) == 1 || size(rxnNameList,2) == 1 
        %reaction index
        objList = sparse(rxnNameList, 1:numel(rxnNameList), ones(numel(rxnNameList),1),...
            nVar, max(size(rxnNameList)));
    else
        error('Invalid numerical input of rxnNameList.');
    end
elseif iscell(rxnNameList)
    objList = sparse(nVar, numel(rxnNameList));
    for jRxnName = 1:numel(rxnNameList)
        rJ = findRxnIDs(modelCom,rxnNameList{jRxnName});
        if ~all(rJ)
            error('Invalid names in rxnNameList');
        end
        objList(rJ,jRxnName) = 1;
    end
else
    error('Invalid input of rxnNameList');
end
options.rxnNameList = objList;
clear objList

if isempty(savePOA)
    %always use save option to reduce memory need
    savePOA = ['POAtmp' filesep 'POA'];
end
directory = strsplit(savePOA,filesep);
if numel(directory) > 1
    %not saving in the current directory. Check existence
    directory = strjoin(directory(1:end-1),filesep);
    if ~exist(directory,'dir')
        mkdir(directory);
    end
end
if init
    LPstart = LP.Start;
    LPmodel = LP.Model;
    save([savePOA '_MasterModel.mat'],'LPstart','LPmodel','GRvector', 'kDisp','idRow','GRmax');
end
for j = 1:numel(GRvector)
    optionsJ = options;
    optionsJ.GR = GRvector(j);
    optionsJ.savePOA = sprintf(['%s_GR%.' num2str(kDisp) 'f'], savePOA, GRvector(j));  
    %better reset the model to ensure feasibility
    if j > 1
        load([savePOA '_MasterModel.mat'],'LPstart','LPmodel');
        LP.Model = LPmodel;
        LP.Start = LPstart;
        LP = setCplexParam(LP, solverParam);
    end
    clear LPmodel LPstart
    SteadyComPOAgrCplex(modelCom,optionsJ,solverParam,LP);
end

%collect output from save file
POAtable = cell(Ncheck, Ncheck);
Stat = repmat(struct('cor',0,'r2',0), [Ncheck, Ncheck, numel(GRvector)]);
fluxRange = zeros(Ncheck, 2, numel(GRvector));

for j = 1:numel(GRvector)
    data = load(sprintf(['%s_GR%.' num2str(kDisp) 'f.mat'],savePOA,GRvector(j)), 'POAtable', 'fluxRange', 'Stat');
    for p = 1:Ncheck
        for q = 1:Ncheck
            if isempty(data.POAtable{p,q})
                POAtable{p,q} = [];
            elseif size(data.POAtable{p,q}, 1) == 1
                %single point, no range
                POAtable{p,q}(:, :, j) = repmat(data.POAtable{p,q}, Nstep, 1);
            else
                POAtable{p,q}(:, :, j) = data.POAtable{p,q};
            end
            Stat(p,q,j).cor = data.Stat(p,q).cor;
            Stat(p,q,j).r2 = data.Stat(p,q).r2;
            fluxRange(:,:, j) = data.fluxRange;
        end
    end
end

end
