function varargout = checkSolFeas(LP, sol, maxInfeas, tol, internal)
%Check the feasibility of a solution given a COBRA model structure or a CPLEX dynamic object and a solution
%infeas = checkSolFeas(LP)
%   returns the maximum infeasibility 'infeas' of the solution in the Cplex 
%   object 'LP'. infeas <= 0 means absolutely feasible.
%infeas = checkSolFeas(LP,[],maxInfeas)
%   returns the maximum infeasibility (maxInfeas = true, default) 
%   or a struct of vectors of infeasibility 'infeas' (maxInfeas = false)
%   of the solution in the Cplex object 'LP' with the fields:
%   'con' for infeasibility of constraints
%   'lb'  for infeasibility of lower bounds
%   'ub'  for infeasibility of upper bounds
%   'ind' for infeasibility of indicator constraints. = -Inf if an
%   indicator is not active.
%infeas = checkSolFeas(LP, sol, maxInfeas)
%   returns a vector 'infeas' of maximum infeasibility for each solution in each 
%   column of 'sol' to the LP problem 'LP', which can be a Cplex object, COBRA
%   model, COBRA LP structure or gurobi LP structure if maxInfeas = true.
%   If maxInfeas = false, return a structure of matrix of infeasibility 
%   (e.g. #constraints x #solutions).
%   'sol' can also be a solution structure resulting from gurobi or solveCobraLP
%[infeas, sol] = checkSolFeas(LP, sol, maxInfeas, tol)
%   returns also the solution structure 'sol' if the input 'sol' is a
%   solution structure, with the field sol.stat = 1 if dev < tol, the
%   feasibility tolerance. 'tol' is defaulted at 1e-8

if nargin < 3 || isempty(maxInfeas)
    maxInfeas = true;
end
if nargin < 5
    internal = false;
    %This arguement is for internal use only. Do not call it.
end
persistent varInd
persistent compleInd
persistent rhsInd
persistent conInd
persistent csInd
if ~internal
    [varInd, compleInd, rhsInd, conInd, csInd] = deal([]);
end
if ~isa(LP, 'Cplex')
    %For COBRA or similar LP problem
    %check field in LP
    if isfield(LP, 'b')
        b0 = LP.b;
    elseif isfield(LP, 'rhs')
        b0 = LP.rhs;
    else
        varargout = {[]};
        return
    end
    if isfield(LP,'A') 
        A = LP.A; %COBRA LP problem or gurobi LP problem
    elseif isfield(LP,'S')
        A = LP.S; %COBRA model
    else
        varargout = {[]};
        return
    end
    if isfield(LP, 'sense')
        cs = LP.csense; %gurobi
        E = '=';
        L = '<';
        G = '>';
    elseif isfield(LP,'csense')
        cs = LP.sense; %COBRA
        E = 'E';
        L = 'L';
        G = 'G';
    else
        %COBRA model, no csense, assume all equal
        cs = char('E' * ones(size(A,1),1));
        E = 'E';
        L = 'L';
        G = 'G';
    end
   %check field in sol
    if isstruct(sol)
        sol0 = sol;
        if isfield(sol,'full') 
            sol = sol.full;
        elseif isfield(sol, 'x')
            sol = sol.x;
        else
            varargout = {[]};
            return
        end
    elseif ismatrix(sol) && size(sol,1) == size(A,2)
        sol = sol;
    else
        varargout = {[]};
        return
    end
    
    b = A * sol;
    b0 = b0 * ones(1, size(sol,2));
    if maxInfeas
        infeas = zeros(1,size(sol,2));
        infeas = max([infeas; max(abs(b(cs == E,:) - b0(cs == E,:)),[],1)],[],1);
        infeas = max([infeas; max(b(cs == L,:) - b0(cs == L,:),[],1)], [], 1);
        infeas = max([infeas; max(b0(cs == G,:) - b(cs == G,:),[],1)], [], 1);
        if isfield(LP,'lb')
            infeas = max([infeas; max(LP.lb*ones(1,size(sol,2)) - sol, [],1)], [], 1);
        end
        if isfield(LP,'ub')
            infeas = max([infeas; max(sol - LP.ub*ones(1,size(sol,2)), [],1)], [], 1);
        end
    else
        infeas.con = zeros(size(A,1),size(sol,2));
        infeas.con(cs == E,:) = abs(b(cs == E,:) - b0(cs == E,:));
        infeas.con(cs == L,:) = b(cs == L,:) - b0(cs == L,:);
        infeas.con(cs == G,:) = b0(cs == G,:) - b(cs == G,:);
        infeas.con(infeas.con < 0) = 0;
        if isfield(LP,'lb')
            infeas.lb = LP.lb*ones(1,size(sol,2)) - sol;
            infeas.lb(infeas.lb < 0) = 0;
        else
            infeas.lb = [];
        end
        if isfield(LP,'ub')
            infeas.ub = sol - LP.ub*ones(1,size(sol,2));
            infeas.ub(infeas.ub < 0) = 0;
        else
            infeas.ub = [];
        end
    end
        
    if isstruct(sol0)
        if nargin < 4
            tol = 1e-8;
        end
        if infeas < tol
            sol0.stat = 1;
        end
        varargout = {infeas, sol0};
    else
        varargout = {infeas};
    end
else
    %For Cplex object
    if nargin < 2 || isempty(sol)
        %check the solution in the Cplex object if no solution input
        if isprop(LP,'Solution') && isfield(LP.Solution,'x') && ~isempty(LP.Solution.x)
            sol = LP.Solution.x;
        else
            varargout = {NaN};
            return
        end
    end
 
    %for checking an array of solutions
    %(Each column in the matrix sol is a solution)
    mSize = 30; %optimized parameter for my computer for models of large size (~1e4)
    %10 ~ 100 are good
    if isfield(LP.Model,'indicator') && ~internal
        [varInd, compleInd, rhsInd] = deal(zeros(numel(LP.Model.indicator),1));
        conInd = zeros(numel(LP.Model.indicator),size(LP.Model.A,2));
        csInd = char('E'*ones(1,numel(LP.Model.indicator)));
        for j = 1:numel(LP.Model.indicator)
            varInd(j) = LP.Model.indicator(j).variable;
            compleInd(j) = LP.Model.indicator(j).complemented;
            conInd(j,:) = LP.Model.indicator(j).a';
            rhsInd(j) = LP.Model.indicator(j).rhs;
            csInd(j) = LP.Model.indicator(j).sense;
        end
        
    end
    b = LP.Model.A * sol;
    
    if maxInfeas
        infeas = zeros(1,size(sol, 2));
        if size(sol, 2) > mSize
            %compute mSize solutions each time
            for j = 1:floor(size(sol, 2) / mSize)
                infeas(((j-1)*mSize + 1):(j*mSize)) = checkSolFeas(LP, ...
                    sol(:,((j-1)*mSize + 1):(j*mSize)),maxInfeas,[],true);
            end
            %compute the remaining solutions
            if mod(size(sol, 2), mSize) > 0
                infeas((floor(size(sol, 2) / mSize)*mSize+1):end) = checkSolFeas(LP, ...
                    sol(:,(floor(size(sol, 2) / mSize)*mSize+1):end),maxInfeas,[],true);
            end
            
        else
            infeas = max([infeas; max(repmat(LP.Model.lhs,1,size(b,2)) - b,[],1)],[],1);
            infeas = max([infeas; max(b - repmat(LP.Model.rhs,1,size(b,2)),[],1)],[],1);
            infeas = max([infeas; max(repmat(LP.Model.lb,1,size(sol,2)) - sol,[],1)],[],1);
            infeas = max([infeas; max(sol - repmat(LP.Model.ub,1,size(sol,2)),[],1)],[],1);
            if isfield(LP.Model,'indicator')
                compleInd2 = compleInd * ones(1, size(sol,2));
                bInd = (conInd * sol) .* ...
                    (sol(varInd,:) > 0.9 & compleInd2 == 0 | ... %active indicator
                    sol(varInd,:) < 1e-2 & compleInd2 == 1);
                b0Ind = (rhsInd * ones(1, size(sol,2))) .* ...
                    (sol(varInd,:) > 0.9 & compleInd2 == 0 | ... %active indicator
                    sol(varInd,:) < 1e-2 & compleInd2 == 1);
                infeas = max([infeas; max(abs(bInd(csInd == 'E',:) - b0Ind(csInd == 'E',:)),[],1)],[],1);
                infeas = max([infeas; max(bInd(csInd == 'L',:) - b0Ind(csInd == 'L',:),[],1)], [], 1);
                infeas = max([infeas; max(b0Ind(csInd == 'G',:) - bInd(csInd == 'G',:),[],1)], [], 1);
            end
        end
    else
        infeas =struct();
        [infeas.con, infeas.lb, infeas.ub, infeas.ind] = deal([]);
        list = {'con','lb','ub','ind'};
        if size(sol, 2) > mSize
            %compute mSize solutions each time
            for j = 1:floor(size(sol, 2) / mSize)
                infeasJ = checkSolFeas(LP, ...
                    sol(:,((j-1)*mSize + 1):(j*mSize)),maxInfeas,[],true);
                for k = list
                    infeas.(k{:}) = [infeas.(k{:}), infeasJ.(k{:})];
                end
            end
            %compute the remaining solutions
            if mod(size(sol, 2), mSize) > 0
                infeasJ = checkSolFeas(LP, ...
                    sol(:,(floor(size(sol, 2) / mSize)*mSize+1):end),maxInfeas,[],true);
                for k = list
                    infeas.(k{:}) = [infeas.(k{:}), infeasJ.(k{:})];
                end
            end
        else
            infeas.con = repmat(LP.Model.lhs,1,size(b,2)) - b;
            infeas.con = max(infeas.con, b - repmat(LP.Model.rhs,1,size(b,2)));
            infeas.con(infeas.con < 0) = 0;
            infeas.lb = repmat(LP.Model.lb,1,size(sol,2)) - sol;
            infeas.lb(infeas.lb < 0) = 0;
            infeas.ub = sol - repmat(LP.Model.ub,1,size(sol,2));
            infeas.ub(infeas.ub < 0) = 0;
            if isfield(LP.Model,'indicator')
                compleInd2 = compleInd * ones(1, size(sol,2));
                bInd = (conInd * sol);
                b0Ind = rhsInd * ones(1, size(sol,2));
                infeas.ind = zeros(size(conInd,1), size(sol,2));
                infeas.ind(csInd == 'E',:) = abs(bInd(csInd == 'E',:) - b0Ind(csInd == 'E',:));
                infeas.ind(csInd == 'L',:) = bInd(csInd == 'L',:) - b0Ind(csInd == 'L',:);
                infeas.ind(csInd == 'G',:) = b0Ind(csInd == 'G',:) - bInd(csInd == 'G',:);
                infeas.ind(infeas.ind < 0) = 0;
                infeas.ind(~((sol(varInd,:) > 0.9 & compleInd2 == 0) | ... %active indicator
                    (sol(varInd,:) < 1e-2 & compleInd2 == 1))) = -inf;
            else
                infeas.ind = [];
            end
        end
    end
    varargout = {infeas};
end

end

