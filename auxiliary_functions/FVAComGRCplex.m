function [minFlux,maxFlux,minFD,maxFD,LP,GR] = FVAComGRCplex(modelCom,options,solverParam,LP)
%Flux variability analysis for community model at community steady-state at a given growth rate. 
%The function is capable of saving intermediate results and continuing from previous results 
%if the file path is given in options.saveFVA. It also allows switch from single thread to parallel 
%computation from intermediate results (but not the reverse).
%
%[minFlux,maxFlux,minFD,maxFD,LP,GR] = FVAComGRCplex(modelCom,options,solverParam,LP)
%
%INPUT
% modelCom        a community COBRA model structure with the following extra fields:
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
%   GR              The growth rate at which FVA is performed. If not
%                     given, find the maximum growth rate by SteadyComCplex.m
%   optBMpercent    Only consider solutions that yield at least a certain
%                     percentage of the optimal biomass (Default = 99.99)
%   rxnNameList     List of reactions (IDs or .rxns) for which FVA is performed.
%                     Use a (N_rxns + N_organism) x K matrix for FVA of K
%                     linear combinations of fluxes and/or abundances
%                     (Default = biomass reaction of each species)
%   rxnFluxList     List of reactions (IDs or .rxns) whose fluxes are also returned
%                     (Default = biomass reaction of each species)
%  (the two parameters below are usually determined by solving the problem
%      during the program. Provide them only if you want to constrain the
%      total biomass to a particular value)
%   BMmaxLB         lower bound for the total biomass (default 1)
%   BMmaxUB         upper bound for the total biomass
%  (other parameters)
%   saveFVA         If non-empty, become the filename to save the FVA results
%                   (default empty, not saving)
%   threads         > 1 for explicitly stating the no. of threads used in FVA,
%                   0 or -1 for using all available threads. Default 1.
%   verbFlag        Verbose output. 1 to have waitbar, >1 to have stepwise output
%                   (default 3)
%   loadModel       String of filename to be loaded. If non-empty, load the 
%                   cplex model ('loadModel.mps'), basis ('loadModel.bas') 
%                   and parameters ('loadModel.prm').
%  May add also other parameters in SteadyComCplex for calculating the maximum growth rate.
%
% solverParam       Cplex parameter structure. E.g., struct('simplex',struct('tolerances',struct('feasibility',1e-8)))
% LP                Cplex LP object from, e.g., calling FVAComCplex
%
%OUTPUT
% minFlux       Minimum flux for each reaction
% maxFlux       Maximum flux for each reaction
%
%OPTIONAL OUTPUT
% minFD         #rxnFluxList x #rxnNameList matrix containing the fluxes in
%               options.rxnFluxList corresponding to minimizing each reaction in
%               options.rxnNameList
% maxFD         #rxnFluxList x #rxnNameList matrix containing the fluxes in
%               options.rxnFluxList corresponding to maximizing each reaction in
%               options.rxnNameList
% LP            Cplex LP object
% GR            the growth rate at which FVA is performed

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
param2get = {'GR', 'optBMpercent', 'rxnNameList', 'rxnFluxList', ...
             'GRfx', 'BMmaxLB','BMmaxUB',...
             'threads','verbFlag','loadModel','saveFVA'};
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
    %otherwise use that in the solver
    solverParam.simplex.tolerances.feasibility = feasTol;
end

[m, n] = size(modelCom.S);
nRxnSp = sum(modelCom.indCom.rxnSps > 0); %number of species-specific rxns
nSp = numel(modelCom.indCom.spBm); %number of species

if (~isfield(modelCom,'b'))
    modelCom.b = zeros(size(modelCom.S,1),1);
end

%% setup LP structure
checkBMrow = false;
if isempty(GR)
    %if max growth rate not given, find it and get the LP problem
    options2 = options;
    options2.minNorm = false;
    [~, result,LP] = SteadyComCplex(modelCom, options2,solverParam);
    GR = result.GRmax;
    idRow = size(LP.Model.A,1);
    addRow = false;
elseif nargin < 4
    if ~isempty(loadModel)
        % load solution if given and growth rate is known
        LP = Cplex('fva');
        LP.readModel([loadModel '.mps']);
        LP.readBasis([loadModel '.bas']);
        LP.readParam([loadModel '.prm']);
        fprintf('Load model ''%s'' successfully.\n', loadModel);
        checkBMrow = true;
    else
        %get LP using SteadyComCplex if only growth rate is given 
        options2 = options;
        options2.LPonly = true;
        [~, ~, LP] = SteadyComCplex(modelCom, options2, solverParam);
        %no constraint on total biomass using LPonly option
        addRow = true;
    end
else
    %GR given as input and LP is supplied, expected when called by fluxVarComCplex
    checkBMrow = true;
end
%Check if a row constraining the sum of biomass exists
if checkBMrow
    if size(LP.Model.A,1) > m + 2*nRxnSp + nSp
        [ynRow,idRow] = ismember(sparse(ones(nSp,1),n+1:n+nSp,ones(nSp,1),1,n+nSp),...
            LP.Model.A(m+2*nRxnSp+nSp+1:end,1:n+nSp),'rows');
        if ynRow
            idRow = m + 2*nRxnSp + nSp + idRow;
        end
        addRow = ~ynRow;
    else
        addRow = true;
    end
end
%add a row for constraining the sum of biomass if not exist
if addRow
    %using default BMmaxLB and BMmaxUB if not given in options
    LP.addRows(BMmaxLB * optBMpercent / 100, ...
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
    LP.Model.lhs(idRow) = BMmaxLB * optBMpercent / 100;
    %not allow the max. biomass to exceed the one at max growth rate,
    %can happen if optBMpercent < 100. May dismiss this constraint or 
    %manually supply BMmaxUB in the options if sum of biomass should be variable
    LP.Model.rhs(idRow) = BMmaxUB;
end
%set Cplex parameters
LP = setCplexParam(LP, solverParam);
nVar = size(LP.Model.A,2);

BMmax0 = LP.Model.lhs(idRow);
%update the LP to ensure the current growth rate is constrained
LP.Model.A = updateLPcom(modelCom, GR, GRfx, [], LP.Model.A, []);
LP.Model.sense = 'minimize';
LP.Model.obj(:) = 0;
LP.solve();
%check and adjust for feasibility
dev = checkSolFeas(LP);
kBMadjust = 0;
while (~isfield(LP.Solution, 'x') || dev > feasTol) && kBMadjust < 10
    kBMadjust = kBMadjust + 1;
    LP.Model.lhs(end) = BMmax0 * (1 - feasTol/(11 - kBMadjust));
    LP.solve();
    dev = checkSolFeas(LP);
    if verbFlag
        fprintf('BMmax adjusment: %d\n',kBMadjust);
    end
end

if (~isfield(LP.Solution, 'x') || dev > feasTol)
    error('Model not feasible.')
end
BMmax0 = LP.Model.lhs(idRow);

%% handle variables for FVA (objective matrix) and fluxes to return
%fluxes to return
if isnumeric(rxnFluxList)
    rxnFluxId = rxnFluxList;
elseif iscell(rxnFluxList) || ischar(rxnFluxList)
    rxnFluxId = findRxnIDs(modelCom,rxnFluxList);
    if any(rxnFluxId) == 0
        error('Invalid names in rxnFluxList.');
    end
end

%objective matrix
if ischar(rxnNameList)
    %if input is a string, make it a cell
    rxnNameList = {rxnNameList};
end
if isnumeric(rxnNameList)
    %if input is numeric
    if size(rxnNameList,1) >= n && size(rxnNameList,1) <= nVar
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
    %handle cell input
    objList = sparse(nVar, numel(rxnNameList));
    for jRxnName = 1:numel(rxnNameList)
        %each rxnNameList{jRxnName} can be a cell array of reactions
        %treat it as unweighted sum of the reactions
        rJ = findRxnIDs(modelCom,rxnNameList{jRxnName});
        if ~all(rJ)
            error('Invalid names in rxnNameList');
        end
        objList(rJ,jRxnName) = 1;
    end
else
    error('Invalid input of rxnNameList');
end

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
if ~isempty(saveFVA)
    directory = strsplit(saveFVA,filesep);
    if numel(directory) > 1
        %not saving in the current directory. Check existence
        directory = strjoin(directory(1:end-1),filesep);
        if ~exist(directory,'dir')
            mkdir(directory);
        end
    end
end
if verbFlag
    fprintf('\nFVA for %d sets of fluxes/biomass at growth rate %.6f :\n',...
        size(objList,2), GR);
end
%% main loop of FVA
if threads == 1 
    % single-thread FVA
    if (verbFlag == 1)  
        h = waitbar(0,'Flux variability analysis in progress ...');
    end
    if (verbFlag > 1)
        fprintf('%4s\t%4s\t%10s\t%9s\t%9s\n','No','%','Name','Min','Max');
    end

    m = 0;
    maxFlux = zeros(size(objList,2), 1);
    minFlux = zeros(size(objList,2), 1);
    [minFD, maxFD] = deal(sparse(numel(rxnFluxId), size(objList,2)));
    i0 = 0;
    if ~isempty(saveFVA)
        %save the master model
        LPmodel = LP.Model;
        LPstart = LP.Start;
        optionsFVA = options;
        save([saveFVA '_model.mat'], 'LPmodel','LPstart','optionsFVA');
        clear LPmodel LPstart optionsFVA
        %continue from previous saved file
        if ~isempty(saveFVA)
            if exist([saveFVA '.mat'], 'file')
                load([saveFVA '.mat'],'i0','minFlux','maxFlux','minFD','maxFD');
                if i0 == size(objList,2)
                    fprintf('FVA was already finished previously and saved in %s.mat.\n', saveFVA);
                    return
                else
                    fprintf('Continue FVA from i = %d.\n', i0);
                end
            end
        end
    end
    for i = (i0 + 1):size(objList,2)
        if (verbFlag == 1)
            fprintf('iteration %d.  skipped %d\n', i, round(m));
        end
        %maximize
        LP.Model.obj = -full(objList(:,i)); 
        LP.solve();
        while ~isprop(LP, 'Solution')
            try
                LP.solve();
            catch
                %should not happen
                fprintf('Error in solving!')
            end
        end
        %Infeasibility can occur occasionally. Keep trying to relax the maximum biomass
        eps0 = 1e-8;
        while ~(checkSolFeas(LP) <= feasTol) && eps0 * 10 <= 1e-3 %LP.Solution.status ~= 1
            eps0 = eps0 * 10;
            LP.Model.lhs(idRow) = BMmax0 * (1 - eps0);
            LP.solve();
        end
        if ~(checkSolFeas(LP) <= feasTol)
            %infeasible (suggests numerical issues)
            maxFlux(i) = NaN;
            maxFD(:,i) = NaN;
        else
            %LP.Solution.fval can sometimes return NaN even if a solution is found
            maxFlux(i) = -LP.Model.obj' * LP.Solution.x;
            maxFD(:,i) = LP.Solution.x(rxnFluxId);
        end
        %restore the original BMmax0
        LP.Model.lhs(idRow) = BMmax0;
        
        %minimize
        LP.Model.obj = full(objList(:,i));
        LP.solve();
        while ~isprop(LP, 'Solution')
            try
                LP.solve();
            catch
                %should not happen
                fprintf('Error in solving!')
            end
        end
        %Infeasibility can occur occasionally. Keep trying to relax the maximum biomass
        eps0 = 1e-8;
        while ~(checkSolFeas(LP) <= feasTol) && eps0 * 10 <= 1e-3 %LP.Solution.status ~= 1 
            eps0 = eps0 * 10;
            LP.Model.lhs(idRow) = BMmax0 * (1 - eps0);
            LP.solve();
        end
        if ~(checkSolFeas(LP) <= feasTol)
            minFlux(i) = NaN;
            minFD(:,i) = NaN;
        else
            minFlux(i) = LP.Model.obj' * LP.Solution.x;
            minFD(:,i) = LP.Solution.x(rxnFluxId);
        end
        LP.Model.lhs(idRow) = BMmax0;
        
        if (verbFlag == 1)
            waitbar(i/length(rxnNameList),h);
        end
        if (verbFlag > 1)
            rxnNameDisp = strjoin(cellstr(LP.Model.colname(objList(:,i)~=0,:)),' + ');
            fprintf('%4d\t%4.0f\t%10s\t%9.6f\t%9.6f\n',i,100*i/size(objList,2),rxnNameDisp,minFlux(i),maxFlux(i));
        end
        if mod(i, floor(size(objList,2)/50)) == 0 
            if ~isempty(saveFVA)
                %save intermediate results
                i0 = i;
                save([saveFVA '.mat'],...
                    'i0','minFlux','maxFlux','minFD','maxFD')
            end
        end
    end
    if ~isempty(saveFVA)
        %save final results
        i0 = i;
        save([saveFVA '.mat'],...
            'i0','minFlux','maxFlux','minFD','maxFD')
    end
    if (verbFlag == 1)
        if ( regexp( version, 'R20') )
            close(h);
        end
    end
else
    %% parallel FVA
    i0P = 0;
    if ~isempty(saveFVA)
        %check if previous results from single-thread computation exist
        if exist([saveFVA '.mat'], 'file')
            load([saveFVA '.mat'],'i0');
            i0P = i0;
            clear i0
            if i0P == size(objList,2)
                fprintf('FVA was already finished previously and saved in %s.mat.\n', saveFVA);
                load([saveFVA '.mat'],'minFlux','maxFlux','minFD','maxFD');
                return
            else
                fprintf('Continue FVA from i = %d.\n', i0P);
            end
        end
    end
    
    p = gcp;
    numPool = p.NumWorkers;
    %assign reactions to each thread, from reaction i0P to size(objList,2)
    rxnRange = cell(numPool,1);
    remainder = mod(size(objList,2)-i0P,numPool);
    kRxnDist = i0P;
    for jP = 1:numPool
        if jP <= remainder
            rxnRange{jP} = (kRxnDist + 1):(kRxnDist + floor((size(objList,2)-i0P) / numPool) + 1);
        else
            rxnRange{jP} = (kRxnDist + 1):(kRxnDist + floor((size(objList,2)-i0P) / numPool));
        end
        kRxnDist = kRxnDist + numel(rxnRange{jP});
    end
    LPmodel = LP.Model;
    LPstart = LP.Start;
    if ~isempty(saveFVA)
        %save the master model
        optionsFVA = options;
        save([saveFVA '_model.mat'], 'LPmodel','LPstart','optionsFVA','rxnRange','numPool');
    end
    
    [maxFluxCell, minFluxCell, minFDCell, maxFDCell] = deal(cell(numPool,1));
    fprintf('%s\n',saveFVA);
    %for technical reasons, declare the variables for Matlab
    save('FVAparallelTmpVar.mat','saveFVA','verbFlag');
    data = load('FVAparallelTmpVar.mat');
    saveFVA = data.saveFVA;
    verbFlag = data.verbFlag;
    parfor jP = 1:numPool
        LPp = Cplex('subproblem');
        LPp.Model = LPmodel;
        LPp = setCplexParam(LPp, solverParam);
        LPp.Start = LPstart;
        maxFluxP = zeros(numel(rxnRange{jP}), 1);
        minFluxP = zeros(numel(rxnRange{jP}), 1);
        [minFDP, maxFDP] = deal(sparse(numel(rxnFluxId), numel(rxnRange{jP})));
        iCount = 0; %counter of reactions to go
        iSkip = 0; %previously finished reactions
        if ~isempty(saveFVA)
            %check if previous save of parallel computation exists
            if exist([saveFVA '_thread' num2str(jP) '.mat'], 'file')
                [minFluxP,maxFluxP,minFDP,maxFDP,iSkip] = parLoad(jP,saveFVA);
                fprintf('Thread %d: continue FVA from i = %d.\n', jP, iSkip);
            end
        end
        for i = rxnRange{jP}
            iCount = iCount + 1;
            if i > iSkip
                %maximize
                LPp.Model.obj = -full(objList(:,i));
                LPp.solve();
                while ~isprop(LPp, 'Solution')
                    try
                        LPp.solve();
                    catch
                        %should not happen
                        fprintf('Error in solving!')
                    end
                end
                %Infeasibility can occur occasionally. Keep trying to relax the maximum biomass
                eps0 = 1e-8;
                while ~(checkSolFeas(LPp) <= feasTol) && eps0 * 10 <= 1e-3 %LP.Solution.status ~= 1
                    eps0 = eps0 * 10;
                    LPp.Model.lhs(idRow) = BMmax0 * (1 - eps0);
                    LPp.solve();
                end
                if ~(checkSolFeas(LPp) <= feasTol)
                    %infeasible (suggests numerical issues)
                    maxFluxP(iCount) = NaN;
                    maxFDP(:,iCount) = NaN;
                else
                    %LP.Solution.fval can sometimes return NaN even if a solution is found
                    maxFluxP(iCount) = -LPp.Model.obj' * LPp.Solution.x;
                    maxFDP(:,iCount) = LPp.Solution.x(rxnFluxId);
                end
                LPp.Model.lhs(idRow) = BMmax0;
                
                %minimize
                LPp.Model.obj = full(objList(:,i));
                LPp.solve();
                while ~isprop(LPp, 'Solution')
                    try
                        LPp.solve();
                    catch
                        %should not happen
                        fprintf('Error in solving!')
                    end
                end
                %Infeasibility can occur occasionally. Keep trying to relax the maximum biomass
                eps0 = 1e-8;
                while ~(checkSolFeas(LPp) <= feasTol) && eps0 * 10 <= 1e-3 %LP.Solution.status ~= 1
                    eps0 = eps0 * 10;
                    LPp.Model.lhs(idRow) = BMmax0 * (1 - eps0);
                    LPp.solve();
                end
                if ~(checkSolFeas(LPp) <= feasTol)
                    minFluxP(iCount) = NaN;
                    minFDP(:,iCount) = NaN;
                else
                    minFluxP(iCount) = LPp.Model.obj' * LPp.Solution.x;
                    minFDP(:,iCount) = LPp.Solution.x(rxnFluxId);
                end
                LPp.Model.lhs(idRow) = BMmax0;

                if mod(iCount, ceil(numel(rxnRange{jP})/10)) == 0 
                    if (verbFlag)
                        fprintf('Thread %d:\t%.2f%% finished. %04d-%02d-%02d %02d:%02d:%02.0f\n',...
                            jP, iCount / numel(rxnRange{jP}) * 100, clock);
                    end
                    %save intermediate data
                    if ~isempty(saveFVA)
                        parSave(minFluxP,maxFluxP,minFDP,maxFDP,i,jP,saveFVA)
                    end
                end
            end
        end
        if ~isempty(saveFVA)
            %save finished data for each thread
            parSave(minFluxP,maxFluxP,minFDP,maxFDP,i,jP,saveFVA)
        end
        maxFluxCell{jP} = maxFluxP;
        minFluxCell{jP} = minFluxP;
        maxFDCell{jP} = maxFDP;
        minFDCell{jP} = minFDP;
    end
    %collect all results
    [maxFlux, minFlux, minFD, maxFD] = deal([]);
    if exist([saveFVA '.mat'], 'file')
        load([saveFVA '.mat'],'minFlux','maxFlux','minFD','maxFD'); 
    end
    for jP = 1:numPool
        maxFlux = [maxFlux; maxFluxCell{jP}];
        minFlux = [minFlux; minFluxCell{jP}];
        maxFD = [maxFD maxFDCell{jP}];
        minFD = [minFD minFDCell{jP}];
        maxFluxCell{jP} = [];
        minFluxCell{jP} = [];
        maxFDCell{jP} = [];
        minFDCell{jP} = [];
    end
    if ~isempty(saveFVA)
        i0 = size(objList,2);
        save([saveFVA '.mat'],'i0','minFlux','maxFlux','minFD','maxFD');
    end
end

maxFlux = columnVector(maxFlux);
minFlux = columnVector(minFlux);
end

function parSave(minFluxP,maxFluxP,minFDP,maxFDP,i0,jP,saveFVA)
save([saveFVA '_thread' num2str(jP) '.mat'],...
    'i0','minFluxP','maxFluxP','minFDP','maxFDP','jP')
end

function [minFluxP,maxFluxP,minFDP,maxFDP,i0] = parLoad(jP,saveFVA)
load([saveFVA '_thread' num2str(jP) '.mat'],...
    'i0','minFluxP','maxFluxP','minFDP','maxFDP')
end
