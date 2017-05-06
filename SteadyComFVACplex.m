function [minFlux,maxFlux,minFD,maxFD, GRvector, result,LP] = SteadyComFVACplex(modelCom,options,solverParam)
%Flux variability analysis for community model at community steady-state for a range of growth rates. 
%The function is capable of saving intermediate results and continuing from previous results 
%if the file path is given in options.saveFVA. It also allows switch from single thread to parallel 
%computation from intermediate results (but not the reverse).
%
%[minFlux,maxFlux,Vmin,Vmax] = SteadyComFVACplex(modelCom,options,solverParam)
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
% options (optional) structure with the following fields:
%   optGRpercent    A vector of percentages. Perform FVA at these percents
%                     of max. growth rate respectively (Default = 99.99)
%   optBMpercent    Only consider solutions that yield at least a certain
%                     percentage of the optimal biomass (Default = 99.99)
%   rxnNameList     List of reactions (IDs or .rxns) for which FVA is performed.
%                     Use a (N_rxns + N_organism) x K matrix for FVA of K
%                     linear combinations of fluxes and/or abundances
%                     (Default = biomass reaction of each species)
%   rxnFluxList     List of reactions (IDs or .rxns) whose fluxes are also returned
%                     (Default = biomass reaction of each species)
%   GRmax           maximum growth rate of the model (default to be found
%                     SteadyComCplex.m)
%  (the two parameters below are usually determined by solving the problem
%      during the program. Provide them only if you want to constrain the
%      total biomass to a particular value)
%   BMmaxLB         lower bound for the total biomass (default 1)
%   BMmaxUB         upper bound for the total biomass
%   (other parameters)
%   saveFVA         If non-empty, become the filename to save the FVA results
%                   (default empty, not saving)
%   threads         > 1 for explicitly stating the no. of threads used,
%                   0 or -1 for using all available threads. Default 1.
%   verbFlag        Verbose output. 1 to have waitbar, >1 to have stepwise output
%                   (default 3)
%   loadModel       String of filename to be loaded. If non-empty, load the 
%                   cplex model ('loadModel.mps'), basis ('loadModel.bas') 
%                   and parameters ('loadModel.prm').
%  May add also other parameters in SteadyComCplex for calculating the maximum growth rate.
%
% solverParam       Cplex parameter structure. E.g., struct('simplex',struct('tolerances',struct('feasibility',1e-8)))
%
%OUTPUT
% minFlux       Minimum flux for each reaction
% maxFlux       Maximum flux for each reaction
%
% OPTIONAL OUTPUT
% minFD         #rxnFluxList x #rxnNameList matrix containing the fluxes in
%               options.rxnFluxList corresponding to minimizing each reaction in
%               options.rxnNameList
% maxFD         #rxnFluxList x #rxnNameList matrix containing the fluxes in
%               options.rxnFluxList corresponding to maximizing each reaction in
%               options.rxnNameList
% GRvector      A vector of growth rates at which FVA has been performed
% result        result structure from SteadyComCplex
% LP            Cplex LP object
%

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
param2get = {'GRmax', 'optGRpercent', 'rxnNameList', 'rxnFluxList',...
             'GRfx','BMmaxLB','BMmaxUB', ...
             'verbFlag', 'loadModel','saveFVA','threads'};
eval(sprintf('[%s] = getCobraComParams(param2get, options, modelCom);', ...
            strjoin(param2get, ',')...
            )...
    );

[feasTol, ~] = getCobraSolverParams('LP',{'feasTol'; 'optTol'}, solverParam);
if isfield(solverParam,'simplex') && isfield(solverParam.simplex, 'tolerances')...
        && isfield(solverParam.simplex.tolerances,'feasibility')
    %override the feasTol in CobraSolverParam if given in solverParam
    feasTol = solverParam.simplex.tolerances.feasibility;
else
    %otherwise use the feasTol in COBRA toolbox
    solverParam.simplex.tolerances.feasibility = feasTol;
end

[m, n] = size(modelCom.S);
nSp = numel(modelCom.indCom.spBm); %number of species
nRxnSp = sum(modelCom.indCom.rxnSps > 0); %number of species-specific rxns

if ischar(rxnNameList)
    rxnNameList = {rxnNameList};
end
if iscell(rxnNameList)
    nRxnFVA = numel(rxnNameList);
else
    nRxnFVA = size(rxnNameList,2);
end
if ischar(rxnFluxList)
    rxnFluxList = {rxnFluxList};
end
%get maximum growth rate
addRow = false;
GRgiven = false;
if isempty(GRmax)
    if exist('Cplex.p','file') == 6
        [~, result,LP] = SteadyComCplex(modelCom, options, solverParam);
        
    else
        warning('Support Cplex only right now.');
        return
        %need further achitecture for using COBRA solver
    end
    if strcmp(result.stat,'infeasible')
        %infeasible model
        warning('Model is infeasible.');
        [minFlux,maxFlux] = deal(NaN(nRxnFVA,1));
        [minFD,maxFD] = deal(NaN(numel(rxnFluxList), nRxnFVA));
        GRvector = NaN(numel(optGRpercent), 1);
        return
    end
    GRmax = result.GRmax;
    idRow = size(LP.Model.A,1); %row that constrains total biomass
else
    %If GRmax is given, BMmaxLB and BMmaxUB should be included in options in this case to ensure feasibility
    if ~isempty(loadModel)
        % load solution if given and growth rate is known
        LP = Cplex('fluxSampling');
        LP.readModel([loadSol '.mps']);
        LP.readBasis([loadSol '.bas']);
        LP.readParam([loadSol '.prm']);
        fprintf('Load model ''%s'' successfully.\n', loadModel);
        addRow = true;
        if size(LP.Model.A,1) > m + 2*nRxnSp + nSp
            %try to find the row that constrains total biomass
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
    result = struct('GRmax',GRmax,'vBM',[],'BM',[],'Ut',[],'Ex',[],'flux',[],'iter0',[],'iter',[],'stat','optimal');
    GRgiven = true;
end
if addRow
    %add a row for constraining the sum of biomass if not exist
    %using default BMmaxLB and BMmaxUB if not given in options
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
%update the LP to ensure the current growth rate is constrained
LP.Model.A = updateLPcom(modelCom, GRmax, GRfx, [], LP.Model.A, []);
LP.Model.sense = 'minimize';
LP.Model.obj(:) = 0;
LP.solve();
%check and adjust for feasibility
%(LP from SteadyComCplex should pass this automatically as the row has
% been added in SteadyComCplex)
dev = checkSolFeas(LP);
kBMadjust = 0;
while (~isfield(LP.Solution, 'x') || dev > feasTol) && kBMadjust < 10
    kBMadjust = kBMadjust + 1;
    %the row of biomass constraint should at the end
    LP.Model.lhs(idRow) = BMmaxLB * (1 - feasTol/(11 - kBMadjust));
    LP.solve();
    dev = checkSolFeas(LP);
    if verbFlag
        fprintf('BMmax adjusment: %d\n',kBMadjust);
    end
end
if (~isfield(LP.Solution, 'x') || dev > feasTol)
    warning('Model not feasible.')
    [minFlux,maxFlux] = deal(NaN(nRxnFVA,1));
    [minFD,maxFD] = deal(NaN(numel(rxnFluxList), nRxnFVA));
    GRvector = NaN(numel(optGRpercent), 1);
    result.stat = 'infeasible';
    return
end
if GRgiven
    %assign result structure if in the rare case of given maximum growth
    %rate
    result.vBM = LP.Solution.x(modelCom.indCom.spBm);
    result.BM = LP.Solution.x(n+1:n+nSp);
    result.Ut = LP.Solution.x(modelCom.indCom.EXcom(:,1));
    result.Ex = LP.Solution.x(modelCom.indCom.EXcom(:,2));
    result.flux = LP.Solution.x(1:n);
end
if ~isfield(options, 'BMmaxLB')
    options.BMmaxLB = LP.Model.lhs(idRow);
end
if ~isfield(options, 'BMmaxUB')
    options.BMmaxUB = LP.Model.rhs(idRow);
end

GRvector = GRmax * optGRpercent/100;
if ~isempty(saveFVA)
    %decide number of digits in the save name
    if numel(optGRpercent) == 1
        kDisp = 2;
    else
        d = min(GRvector(2:end) - GRvector(1:end-1));
        if d < 1
            kDisp = abs(floor(log10(d)));
        else
            kDisp = 0;
        end
    end
end

[minFlux, maxFlux] = deal(zeros(nRxnFVA, numel(GRvector)));
[minFD, maxFD] = deal(zeros(numel(rxnFluxList), nRxnFVA, numel(GRvector)));

%parallel computation
p = gcp('nocreate');
if isempty(p)
    if threads > 1
        %given explicit no. of threads
        parpool(ceil(threads));
    elseif threads ~= 1
        %default max no. of threads (input 0 or -1 etc)
        parpool;
    end
end
%perform FVA at each growth rate
for j = 1:numel(GRvector)
    optionsJ = options;
    optionsJ.GR = GRvector(j);
    if ~isempty(saveFVA)
        optionsJ.saveFVA = sprintf(['%s_GR%.' num2str(kDisp) 'f'], saveFVA, GRvector(j));  
    end
    [minFluxJ,maxFluxJ,minFDj,maxFDj,LP] = SteadyComFVAgrCplex(modelCom,optionsJ, solverParam,LP);
    minFlux(:, j) = minFluxJ;
    maxFlux(:, j) = maxFluxJ;
    minFD(:,:, j) = minFDj;
    maxFD(:,:, j) = maxFDj;
end

end
