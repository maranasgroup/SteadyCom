function [POAtable, fluxRange, Stat, pairList] = SteadyComPOAgrCplex(modelCom,options,solverParam,LP)
%Pairwise POA for community model at community steady-state at a given growth rate
%[POAtable, fluxRange, Stat, pairList] = SteadyComPOAgrCplex(modelCom,options,solverParam,LP)
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
%   GR              The growth rate at which POA is performed. If not
%                     given, find the maximum growth rate.
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
%                   Defaulted to be found by SteadyComFVACplex.m
%  (other parameters)
%   savePOA         Must be non-empty. The filename to save the POA results
%                   (default 'POAtmp/POA')
%   threads         > 1 for explicitly stating the no. of threads used,
%                   0 or -1 for using all available threads. Default 1.
%   verbFlag        Verbose output. 0 or 1.
%   loadModel       String of filename to be loaded. If non-empty, load the 
%                   cplex model ('loadModel.mps'), basis ('loadModel.bas') 
%                   and parameters ('loadModel.prm').
%  May add also other parameters in SteadyComCplex for calculating the maximum growth rate.
%
%OUTPUT
% POAtable          K x K cells. 
%                   (i,i)-cell contains the flux range of rxnNameList{i}
%                   (i,j)-cell contains a Nstep x 2 matrix, with (k,1)-entry 
%                   being the min of rxnNameList{j} when rxnNameList{i} is 
%                   fixed at the k-th value, (k,2)-entry being the max.
% fluxRange         K x 2 matrix of flux range for each entry in rxnNameList 
% Stat              K x K structure array with fields:
%                     -'cor': the slope from linear regression between the
%                             fluxes of a pair
%                     -'r2':  the corresponding coefficient of determination (R-square)
% pairList          pairList after transformation from various input formats

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
param2get = {'GRfx','GR', 'BMmaxLB','BMmaxUB','optBMpercent',... %parameters for finding maximum growth rate
             'symmetric','rxnNameList','pairList', 'fluxRange', 'Nstep', 'NstepScale',... %parameters for POA
             'verbFlag', 'threads', 'savePOA','loadModel',...
             };
eval(sprintf('[%s] = getCobraComParams(param2get, options, modelCom);', ...
    strjoin(param2get, ',')...
    )...
    );
if isempty(savePOA)
    %always use save option to reduce memory need
    savePOA = 'POAtmp/POA';
end
directory = strsplit(savePOA,filesep);
if numel(directory) > 1
    %not saving in the current directory. Check existence
    directory = strjoin(directory(1:end-1),filesep);
    if ~exist(directory,'dir')
        mkdir(directory);
    end
end
[feasTol, ~] = getCobraSolverParams('LP',{'feasTol'; 'optTol'}, solverParam);
if isfield(solverParam,'simplex') && isfield(solverParam.simplex, 'tolerances')...
        && isfield(solverParam.simplex.tolerances,'feasibility')
    %override the feasTol in CobraSolverParam if given in solverParam
    feasTol = solverParam.simplex.tolerances.feasibility;
else
    %otherwise use that in the solver
    solverParam.simplex.tolerances.feasibility = feasTol;
end

%Check if the whole computation been finished before
if ~isempty(savePOA)
    if exist(sprintf('%s.mat',savePOA), 'file')
        data0 = load(sprintf('%s.mat',savePOA));
        if isfield(data0, 'finished')
            if verbFlag
                fprintf('Already finished. Results loaded from %s.mat\n',savePOA);
            end
            load(sprintf('%s.mat',savePOA), 'POAtable', 'fluxRange', 'Stat')
            return
        else
            clear data0
        end
    end
end

%parallel computation
if isempty(gcp('nocreate'))
    if threads > 1
        %given explicit no. of threads
        parpool(ceil(threads));
    elseif threads ~= 1
        %default max no. of threads (input 0 or -1 etc)
        parpool;
    end
end
%sizes
[m, n] = size(modelCom.S);
nRxnSp = sum(modelCom.indCom.rxnSps > 0); %number of species-specific rxns
nSp = numel(modelCom.indCom.spBm); %number of species

%% handle LP structure
checkBMrow = false;
if isempty(GR)
    %if max growth rate not given, find it and get the LP problem
    options2 = options;
    options2.minNorm = false;
    [~, result,LP] = SteadyComCplex(modelCom, options2, solverParam);
    GR = result.GRmax;
    idRow = size(LP.Model.A,1);
    addRow = false;
elseif nargin < 4
    if ~isempty(loadModel)
        % load solution if given and growth rate is known
        LP = Cplex('poa');
        LP.readModel([loadSol '.mps']);
        LP.readBasis([loadSol '.bas']);
        LP.readParam([loadSol '.prm']);
        fprintf('Load model ''%s'' successfully.\n', loadModel);
        checkBMrow = true;
    else
        %get LP using optimizeCbModelComCplex if only growth rate is given 
        options2 = options;
        options2.LPonly = true;
        [~, ~, LP] = SteadyComCplex(modelCom, options2, solverParam);
        addRow = true;
    end
else
    checkBMrow = true;
end

if checkBMrow && size(LP.Model.A,1) > m + 2*nRxnSp + nSp
    %if a row constraining the sum of biomass exists
    [ynRow,idRow] = ismember(sparse(ones(nSp,1),n+1:n+nSp,ones(nSp,1),1,n+nSp),...
            LP.Model.A(m+2*nRxnSp+nSp+1:end,1:n+nSp),'rows');
    if ynRow
        idRow = m + 2*nRxnSp + nSp + idRow;
    end
    addRow = ~ynRow;
end
if addRow
    %add a row for constraining the sum of biomass if not exist
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
BMmax0 = LP.Model.lhs(idRow);
%update the LP to ensure the current growth rate is constrained
LP.Model.A = updateLPcom(modelCom, GR, GRfx, [], LP.Model.A, []);
LP.Model.sense = 'minimize';
LP.Model.obj(:) = 0;
%number of variables
nVar = size(LP.Model.A,2);

%% handle objective matrix
if ischar(rxnNameList)
    %treat single character vector as cell with one element
    rxnNameList = {rxnNameList};
end
if isnumeric(rxnNameList)
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
%number of variables for POA
Ncheck = size(objList,2);

%% Handle pairList
if isempty(pairList)
    %if pairList not given, run for all pairs
    pairList = [reshape(repmat(1:Ncheck,Ncheck,1),Ncheck^2,1), repmat((1:Ncheck)',Ncheck,1)];
    if symmetric 
        %the option 'symmetric' is only used here when pairList is not
        %supplied to avoid running symmetric pairs (e.g. j vs k and k vs j)
        pairList(pairList(:,1) > pairList(:,2),:) = [];
    end
elseif numel(size(pairList)) ~= 2 || size(pairList,2) < 2 
    error('pairList must be an N-by-2 array denoting the pairs (rxn names or indices in rxnNameList) to analyze!')
else
    if iscell(pairList) 
        if ~iscell(rxnNameList)
            error('rxnNameList must be cell array of rxn names to match the rxn names in pairList!');
        else
            if iscellstr(pairList) && iscellstr(rxnNameList)
                %both are character arrays
                [~,pairId] = ismember(pairList, rxnNameList);
            else
                %Need to use isequal
                pairId = zeros(size(pairList));
                for jP = 1:numel(pairList)
                    for jR = 1:numel(rxnNameList)
                        if isequal(pairList{jP},rxnNameList{jR})
                            pairId(jP) = jR;
                            break
                        end
                    end
                end
            end
            if ~all(all(pairId))
                error('Some entries in pairList cannot be mapped to rxnNameList!')
            end
            pairList = pairId;
            clear pairId
        end
    end  
end
Npair = size(pairList,1);

%% find pairs in the pairList not analyzed yet
undone = true(Npair,1);
if ~isempty(savePOA) 
    if exist(sprintf('%s.mat',savePOA),'file')
        fprintf('Already finished. Results were already saved to %s.mat\n',savePOA);
        load(sprintf('%s.mat',savePOA), 'POAtable', 'fluxRange', 'Stat');
        return
    elseif exist(sprintf('%s_POApre.mat',savePOA),'file')
        fluxRange = load(sprintf('%s_POApre.mat',savePOA),'fluxRange');
        fluxRange = fluxRange.fluxRange;
        for jP = 1:Npair
            undone(jP) = ~exist(sprintf('%s_j%d_k%d.mat',savePOA,pairList(jP,1),pairList(jP,2)), 'file');
        end
        fprintf('Unfinished pairs: %d\n', sum(undone));
    end
end

if any(undone)
    %% initial solve to ensure feasibility
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
    %terminate if not feasible
    if (~isfield(LP.Solution, 'x') || dev > feasTol)
        error('Model not feasible.')
    end
    %% Find flux range by FVA
    if isempty(fluxRange)
        optionsFVA = options;
        optionsFVA.rxnNameList = objList;
        optionsFVA.GR = GR;
        fluxRange = zeros(size(objList,2),2);
        [fluxRange(:,1), fluxRange(:,2)] = SteadyComFVAgrCplex(modelCom,optionsFVA,solverParam,LP);
        printFluxRange = false;
        if ~isempty(savePOA)
            save(sprintf('%s_POApre.mat',savePOA),'fluxRange','GR','options',...
                'objList','pairList','options','solverParam');
        end
    else
        printFluxRange = true;
    end
    %print flux range
    if verbFlag && printFluxRange
        fprintf('Flux range:\n');
        fprintf('rxn\tmin\tmax\n');
        for jRxn = 1:size(objList,2)
            strPrint = strjoin(strtrim(cellstr(LP.Model.colname(objList(:,jRxn)~=0,:))),',');
            fprintf('%s\t%.6f\t%.6f\n', strPrint, fluxRange(jRxn,1), fluxRange(jRxn,2));
        end
    end
    if (~isfield(modelCom,'b'))
        modelCom.b = zeros(size(modelCom.S,1),1);
    end
    
    Npair =sum(undone);
    undoneId = find(undone);
    
    %print starting point
    if verbFlag
        fprintf('\nPOA for %d pairs of reactions at growth rate %.6f\n', ...
            Npair - sum(pairList(undoneId,1)==pairList(undoneId,2)), GR);
        strPrintj = strjoin(strtrim(cellstr(LP.Model.colname(objList(:,pairList(undoneId(1),1))~=0,:))),',');
        strPrintk = strjoin(strtrim(cellstr(LP.Model.colname(objList(:,pairList(undoneId(1),2))~=0,:))),',');
        fprintf('Start from #%d %s vs #%d %s.\n', pairList(undoneId(1),1), ...
            strPrintj, pairList(undoneId(1),2), strPrintk);
        fprintf('%15s%15s%10s%10s%10s%10s   %s\n','Rxn1','Rxn2','corMin','r2','corMax','r2','Time')
    end
    lb0 = LP.Model.lb;
    ub0 = LP.Model.ub;
    addConstraint = false;

    if threads == 1
        %% single thread computation
        %Get data points to be computated
        
        for jP = 1:Npair
            [j,k] = deal(pairList(undoneId(jP),1),pairList(undoneId(jP),2));
            if ~exist(sprintf('%s_j%d_k%d.mat',savePOA,j,k),'file')
                if abs(fluxRange(j,2) - fluxRange(j,1)) < 1e-8
                    fluxRangeJ = fluxRange(j, 1);
                    NstepJK = 1;
                else
                    if numel(Nstep) > 1
                        %manually supply Nstep vector (% from min to max)
                        NstepJK = numel(Nstep);
                        fluxRangeJ = fluxRange(j,1) + (fluxRange(j,2) - fluxRange(j,1)) * Nstep;
                    else
                        %uniform step or log-scaled step
                        NstepJK = Nstep;
                        if strcmp(NstepScale, 'log')
                            if sign(fluxRange(j,2)) == sign(fluxRange(j,1))
                                if fluxRange(j,1) > 0
                                    [a, b] = deal(fluxRange(j,1), fluxRange(j,2));
                                    fluxRangeJ = exp(log(a) + ((log(b) - log(a))/(Nstep - 1)) * (0:(Nstep - 1)));
                                else
                                    [b, a] = deal(-fluxRange(j,1), -fluxRange(j,2));
                                    fluxRangeJ = exp(log(a) + ((log(b) - log(a))/(Nstep - 1)) * (0:(Nstep - 1)));
                                    fluxRangeJ = -fluxRangeJ(end:-1:1);
                                end
                            else
                                %Not an ideal situation. Flux ranges containing zero
                                %should not use step size at log scale 
                                a = [-inf, (1:(Nstep-1))/(Nstep-1)];
                                fluxRangeJ = fluxRange(j,1) + (fluxRange(j,2) - fluxRange(j,1)) * 0.01 * (100 .^ a);
                            end
                        else
                            %uniform step size
                            fluxRangeJ = (fluxRange(j,1) : (fluxRange(j,2) - fluxRange(j,1)) / (Nstep - 1) : fluxRange(j,2))';
                        end
                    end
                    fluxRangeJ = fluxRangeJ(:);
                end
                %delete constraint added in the previous round if any
                if addConstraint
                    LP.delRows(size(LP.Model.A,1));
                    addConstraint = false;
                end
                %if not a single flux, but a linear combination, add explicit constraint.
                if nnz(objList(:,j)) > 1
                    LP.addRows(-inf, objList(:,j)',inf, 'POArow');
                    addConstraint = true;
                else
                    rxnNameId = objList(:,j) ~= 0;
                end
                
                if j == k
                    %Nothing to analyze. Just record the flux range
                    POAtableJK = fluxRangeJ;
                    StatJK.cor = [1 1];
                    StatJK.r2 = [1 1];
                else
                    fluxPOAvalue = zeros(NstepJK, 2);
                    %reset LP bounds
                    LP.Model.lb = lb0;
                    LP.Model.ub = ub0;
                    for p = 1:NstepJK
                        %fix flux of the j-th reaction
                        if addConstraint
                            LP.Model.lhs(end) = fluxRangeJ(p) - 1e-12;
                            LP.Model.rhs(end) = fluxRangeJ(p) + 1e-12;
                        else
                            LP.Model.lb(rxnNameId) = fluxRangeJ(p) - 1e-12;
                            LP.Model.ub(rxnNameId) = fluxRangeJ(p) + 1e-12;
                        end
                        LP.Model.obj(:) = 0;
                        %minimize flux of the k-th reaction
                        LP.Model.obj = objList(:, k);
                        LP.Model.sense = 'minimize';
                        LP.solve();
                        %handle possible tolerance issues
                        dev = checkSolFeas(LP);
                        eps0 = 1e-9;
                        while dev > feasTol && eps0 < 1e-3 %largest acceptable tolerance set to be 0.001
                            eps0 = eps0 * 10;
                            if addConstraint
                                LP.Model.lhs(end) = fluxRangeJ(p) - eps0;
                                LP.Model.rhs(end) = fluxRangeJ(p) + eps0;
                            else
                                LP.Model.lb(rxnNameId) = fluxRangeJ(p) - eps0;
                                LP.Model.ub(rxnNameId) = fluxRangeJ(p) + eps0;
                            end
                            LP.solve();
                            dev = checkSolFeas(LP);
                        end
                        if dev <= feasTol
                            fluxPOAvalue(p, 1) = LP.Model.obj'*LP.Solution.x;
                        else
                            fluxPOAvalue(p, 1) = NaN; %shoud not happen
                        end
                        %maximize flux of the k-th reaction
                        LP.Model.sense = 'maximize';
                        LP.solve();
                        dev = checkSolFeas(LP);
                        eps0 = 1e-9;
                        while dev > feasTol && eps0 < 1e-3
                            eps0 = eps0 * 10;
                            if addConstraint
                                LP.Model.lhs(end) = fluxRangeJ(p) - eps0;
                                LP.Model.rhs(end) = fluxRangeJ(p) + eps0;
                            else
                                LP.Model.lb(rxnNameId) = fluxRangeJ(p) - eps0;
                                LP.Model.ub(rxnNameId) = fluxRangeJ(p) + eps0;
                            end
                            LP.solve();
                            dev = checkSolFeas(LP);
                        end
                        if dev <= feasTol
                            fluxPOAvalue(p, 2) = LP.Model.obj'*LP.Solution.x;
                        else
                            fluxPOAvalue(p, 2) = NaN; %shoud not happen
                        end
                        
                    end
                    POAtableJK = fluxPOAvalue;
                    %simple linear regression to check correlations
                    notNan = ~isnan(fluxPOAvalue(:,1));
                    [bMin,~,~,~,statMin] = regress(fluxPOAvalue(notNan,1),...
                        [fluxRangeJ(notNan) ones(sum(notNan),1)]);
                    notNan = ~isnan(fluxPOAvalue(:,2));
                    [bMax,~,~,~,statMax] = regress(fluxPOAvalue(notNan,2),...
                        [fluxRangeJ(notNan) ones(sum(notNan),1)]);
                    StatJK.cor = [bMin(1) bMax(1)];
                    StatJK.r2 = [statMin(1) statMax(1)];
                end
                if verbFlag
                    if j ~= k
                        strPrintj = strjoin(strtrim(cellstr(LP.Model.colname(objList(:,j)~=0,:))),',');
                        strPrintk = strjoin(strtrim(cellstr(LP.Model.colname(objList(:,k)~=0,:))),',');
                        fprintf('%15s%15s%10.4f%10.4f%10.4f%10.4f   %04d-%02d-%02d %02d:%02d:%02.0f\n',...
                            strPrintj, strPrintk, StatJK.cor(1), StatJK.r2(1), ...
                            StatJK.cor(2), StatJK.r2(2), clock);
                    end
                end
                %save
                if ~isempty(savePOA)
                    iSave(savePOA,POAtableJK,StatJK,GR,j,k);
                end
            end
        end
    else
        %% parallel
        fprintf('POA in parallel...\n');
        
        %save and load the variables to ensure the parallel code processor can
        %recognize the variables
        tmpSave = 'POAtmp.mat';
        kTemp = 0;
        while exist(tmpSave,'file')
            kTemp = kTemp + 1;
            tmpSave = ['POAtmp' num2str(kTemp) '.mat'];
        end
        save(tmpSave, 'symmetric','verbFlag','savePOA','fluxRange','Nstep');
        tmpLoad = load(tmpSave);
        verbFlag = tmpLoad.verbFlag;
        savePOA = tmpLoad.savePOA;
        fluxRange = tmpLoad.fluxRange;
        Nstep = tmpLoad.Nstep;
        delete(tmpSave);
        clear tmpLoad
        LPmodel = LP.Model;
        LPstart = LP.Start;
        %parallelization using spmd to allow redistribution of jobs upon 
        %completion of any of the workers to avoid being idle
        numPool = gcp;
        numPool = numPool.NumWorkers;
        while Npair > 0
            %mannually distribute jobs
            remainder = mod(Npair,numPool);
            nJ = floor(Npair/numPool);
            kJ = 0;
            nRange = cell(numPool,1);
            undoneCur = Composite();
            %Composite object is assigned as one cell/worker outside 
            %spmd blocks but called as the cell content inside spmd blocks
            for kP = 1:numPool
                if kP <= remainder
                    nRange{kP} = (kJ + 1) : (kJ + nJ + 1);
                    kJ = kJ + nJ + 1;
                else
                    nRange{kP} = (kJ + 1) : (kJ + nJ);
                    kJ = kJ + nJ;
                end
                nRange{kP} = undoneId(nRange{kP});
                undoneCur{kP} = undone(nRange{kP});
            end
            spmd
                %setup local LP
                LPp = Cplex('subproblem');
                LPp.Model = LPmodel;
                LPp = setCplexParam(LPp, solverParam);
                LPp.Start = LPstart;
                addConstraint = false;
                %denote wether the current flag is the first to finish
                first = true;
                for jP = 1:numel(nRange{labindex}) %labindex = #thread
                    [j,k] = deal(pairList(nRange{labindex}(jP),1),pairList(nRange{labindex}(jP),2));
                    if ~exist(sprintf('%s_j%d_k%d.mat',savePOA,j,k),'file')
                        if abs(fluxRange(j,2) - fluxRange(j,1)) < 1e-8
                            fluxRangeJ = fluxRange(j, 1);
                            NstepJK = 1;
                        else
                            if numel(Nstep) > 1
                                %mannually supplied Nstep vector (% from min to max)
                                NstepJK = numel(Nstep);
                                fluxRangeJ = fluxRange(j,1) + (fluxRange(j,2) - fluxRange(j,1)) * Nstep;
                            else
                                %uniform steps or log-scaled steps
                                NstepJK = Nstep;
                                if strcmp(NstepScale, 'log')
                                    if sign(fluxRange(j,2)) == sign(fluxRange(j,1))
                                        if fluxRange(j,1) > 0
                                            [a, b] = deal(fluxRange(j,1), fluxRange(j,2));
                                            fluxRangeJ = exp(log(a) + ((log(b) - log(a))/(Nstep - 1)) * (0:(Nstep - 1)));
                                        else
                                            [b, a] = deal(-fluxRange(j,1), -fluxRange(j,2));
                                            fluxRangeJ = exp(log(a) + ((log(b) - log(a))/(Nstep - 1)) * (0:(Nstep - 1)));
                                            fluxRangeJ = -fluxRangeJ(end:-1:1);
                                        end
                                    else
                                        %Not an ideal situation. Flux ranges containing zero
                                        %should not use step size at log scale
                                        a = [-inf, (1:(Nstep-1))/(Nstep-1)];
                                        fluxRangeJ = fluxRange(j,1) + (fluxRange(j,2) - fluxRange(j,1)) * 0.01 * (100 .^ a);
                                    end
                                else
                                    %linear step size
                                    fluxRangeJ = (fluxRange(j,1) : (fluxRange(j,2) - fluxRange(j,1)) / (Nstep - 1) : fluxRange(j,2))';
                                end
                            end
                            fluxRangeJ = fluxRangeJ(:);
                        end
                        StatJK = struct();
                        if j == k
                            %Nothing to analyze. Just record the flux range
                            POAtableJK = fluxRangeJ;
                            StatJK.cor = [1 1];
                            StatJK.r2 = [1 1];
                        else
                            %reset LP bounds
                            LPp.Model.lb = lb0;
                            LPp.Model.ub = ub0;
                            %delete constraint added in the previous round if any
                            if addConstraint
                                LPp.delRows(size(LPp.Model.A,1));
                                addConstraint = false;
                            end
                            %if not a single flux, but a linear combination, add explicit constraint.
                            rxnNameId = [];
                            if nnz(objList(:,j)) > 1
                                LPp.addRows(-inf, objList(:,j)',inf, 'POArow');
                                addConstraint = true;
                            else
                                rxnNameId = objList(:,j) ~= 0;
                            end
                            
                            fluxPOAvalue = zeros(NstepJK, 2);
                            for p = 1:NstepJK
                                %fix flux of the j-th reaction
                                if addConstraint
                                    LPp.Model.lhs(end) = fluxRangeJ(p) - 1e-12;
                                    LPp.Model.rhs(end) = fluxRangeJ(p) + 1e-12;
                                else
                                    LPp.Model.lb(rxnNameId) = fluxRangeJ(p) - 1e-12;
                                    LPp.Model.ub(rxnNameId) = fluxRangeJ(p) + 1e-12;
                                end
                                %minimize flux of the k-th reaction
                                LPp.Model.obj = objList(:, k);
                                LPp.Model.sense = 'minimize';
                                LPp.solve();
                                %handle possible tolerance issues
                                dev = checkSolFeas(LPp);
                                eps0 = 1e-9;
                                while dev > feasTol && eps0 < 1e-3
                                    eps0 = eps0 * 10;
                                    if addConstraint
                                        LPp.Model.lhs(end) = fluxRangeJ(p) - eps0;
                                        LPp.Model.rhs(end) = fluxRangeJ(p) + eps0;
                                    else
                                        LPp.Model.lb(rxnNameId) = fluxRangeJ(p) - eps0;
                                        LPp.Model.ub(rxnNameId) = fluxRangeJ(p) + eps0;
                                    end
                                    LPp.solve();
                                    dev = checkSolFeas(LPp);
                                end
                                if dev <= feasTol
                                    minf = LPp.Model.obj'*LPp.Solution.x;
                                else
                                    minf = NaN; %shoud not happen
                                end
                                %maximize flux of the k-th reaction
                                LPp.Model.sense = 'maximize';
                                LPp.solve();
                                %handle possible tolerance issues
                                dev = checkSolFeas(LPp);
                                eps0 = 1e-9;
                                while dev > feasTol && eps0 < 1e-3
                                    eps0 = eps0 * 10;
                                    if addConstraint
                                        LPp.Model.lhs(end) = fluxRangeJ(p) - eps0;
                                        LPp.Model.rhs(end) = fluxRangeJ(p) + eps0;
                                    else
                                        LPp.Model.lb(rxnNameId) = fluxRangeJ(p) - eps0;
                                        LPp.Model.ub(rxnNameId) = fluxRangeJ(p) + eps0;
                                    end
                                    LPp.solve();
                                    dev = checkSolFeas(LPp);
                                end
                                if dev <= feasTol
                                    maxf = LPp.Model.obj'*LPp.Solution.x;
                                else
                                    maxf = NaN; %shoud not happen
                                end
                                fluxPOAvalue(p,:) = [minf maxf];
                            end
                            
                            POAtableJK = fluxPOAvalue;
                            %simple linear regression to check correlations
                            notNan = ~isnan(fluxPOAvalue(:,1));
                            [bMin,~,~,~,statMin] = regress(fluxPOAvalue(notNan,1),...
                                [fluxRangeJ(notNan) ones(sum(notNan),1)]);
                            notNan = ~isnan(fluxPOAvalue(:,2));
                            [bMax,~,~,~,statMax] = regress(fluxPOAvalue(notNan,2),...
                                [fluxRangeJ(notNan) ones(sum(notNan),1)]);
                            StatJK.cor = [bMin(1) bMax(1)];
                            StatJK.r2 = [statMin(1) statMax(1)];
                        end
                        if verbFlag
                            if j ~= k
                                strPrintj = strjoin(strtrim(cellstr(LPp.Model.colname(objList(:,j)~=0,:))),',');
                                strPrintk = strjoin(strtrim(cellstr(LPp.Model.colname(objList(:,k)~=0,:))),',');
                                fprintf('%15s%15s%10.4f%10.4f%10.4f%10.4f   %04d-%02d-%02d %02d:%02d:%02.0f\n',...
                                    strPrintj, strPrintk,...
                                    StatJK.cor(1), StatJK.r2(1), ...
                                    StatJK.cor(2), StatJK.r2(2), clock);
                            end
                        end
                        if ~isempty(savePOA)
                            iSave(savePOA,POAtableJK,StatJK,GR,j,k);
                        end
                    end
                    undoneCur(jP) = false;
                    %check if any of workers has finished its loop, break the
                    %loop and redistribute if finished
                    if labProbe('any',0);
                        first = false;
                        break
                    end
                end
                if first
                    %finish of one worker, call off other workers
                    if verbFlag
                        fprintf('Current loop finished. Stop other workers...\n');
                    end
                    labSend(true,setdiff(1:numlabs,labindex),0);
                    if verbFlag
                        fprintf('All workers have ceased. Redistributing...\n');
                    end
                end
                %avoid warning of missed message
                pause(1e-8);
                while labProbe('any',0);
                    pause(1e-8);
                    labReceive('any',0);
                end
            end
            %update undone
            for kP = 1:numPool
                undone(nRange{kP}) = undoneCur{kP};
            end
            undoneId = find(undone);
            Npair = numel(undoneId);
        end
    end
end
if verbFlag
    fprintf('Finished. Save final results to %s.mat\n',savePOA);
end
POAtable = cell(Ncheck, Ncheck);
Stat = repmat(struct('cor',0,'r2',0), Ncheck, Ncheck);
for jP = 1:size(pairList,1)
    data = load(sprintf('%s_j%d_k%d.mat',savePOA,pairList(jP,1),...
        pairList(jP,2)),'POAtableJK', 'StatJK');
    POAtable{pairList(jP,1),pairList(jP,2)} = data.POAtableJK;
    Stat(pairList(jP,1),pairList(jP,2)) = data.StatJK;
end
save(sprintf('%s.mat',savePOA),'POAtable', 'Stat', 'GR', ...
    'fluxRange');
end

function iSave(savePOA,POAtableJK,StatJK,GR,j0,k0)
save(sprintf('%s_j%d_k%d.mat',savePOA,j0,k0),'POAtableJK', 'StatJK', 'GR','j0','k0');
end