function [D0,infoD0] = SolveD0(Kuc,B0,Ba,opts)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Solve periodic fluctuation matrix D0 in Eq.(9)
% This function is used to replace the original expression:
% D0 = -pinv(B0'*Kuc*B0)*B0'*Kuc*Ba;
%
% The main homogenization process is not changed:
% Da = B0*D0+Ba;
% Ch = (1/V)*Be'*Da'*Kuc*Da*Be;
%
% Method options:
% opts.method = 'auto'        : default, try fast solvers first and use pinv as fallback
% opts.method = 'direct'      : solve by backslash when Kii is well-conditioned
% opts.method = 'lsqminnorm'  : minimum-norm solution, available in MATLAB
% opts.method = 'pinv-single' : single-precision pseudo-inverse, fast but approximate
% opts.method = 'pinv'        : Moore-Penrose pseudo-inverse, robust reference solution
% opts.method = 'gauge'       : fix selected DOFs and solve the reduced system
%
% Note:
% Kii may contain rigid-body modes. Therefore, D0 is not always unique.
% Different solvers may give different D0 matrices, but the final stiffness
% Ch should be almost unchanged if only zero-energy modes are different.
% For strict comparison with the formula in the manuscript, use 'pinv'.

%% Default parameters
if nargin < 4 || isempty(opts)
    opts = struct();
end
if ~isfield(opts,'method') || isempty(opts.method)
    opts.method = 'auto';
end
if ~isfield(opts,'rcondTol') || isempty(opts.rcondTol)
    opts.rcondTol = 1e-12;
end
if ~isfield(opts,'resTol') || isempty(opts.resTol)
    opts.resTol = 1e-8;
end
if ~isfield(opts,'useSingleInAuto') || isempty(opts.useSingleInAuto)
    opts.useSingleInAuto = 1;
end
if ~isfield(opts,'verbose') || isempty(opts.verbose)
    opts.verbose = 0;
end

%% Construct matrices in Eq.(14)
Kuc = sparse(Kuc);
B0 = sparse(B0);
Ba = sparse(Ba);

Kii = B0'*Kuc*B0;
Kia = B0'*Kuc*Ba;

%% Initialization of output information
infoD0.methodRequested = opts.method;
infoD0.methodUsed = '';
infoD0.rcondEst = EstimateRcondKii(Kii);
infoD0.residualRelative = NaN;
infoD0.flag = 0;
infoD0.message = '';
infoD0.tryHistory = {};

methodName = lower(opts.method);

%% Solve D0 according to the selected method
switch methodName
    case 'auto'
        [X,infoD0] = AutoSolveKii(Kii,Kia,opts,infoD0);

    case 'direct'
        [X,infoD0] = DirectSolveKii(Kii,Kia,opts,infoD0);

    case 'lsqminnorm'
        [X,infoD0] = LsqminnormSolveKii(Kii,Kia,opts,infoD0);

    case {'pinv-single','single'}
        [X,infoD0] = SinglePinvSolveKii(Kii,Kia,opts,infoD0);

    case 'pinv'
        [X,infoD0] = PinvSolveKii(Kii,Kia,opts,infoD0);

    case 'gauge'
        [X,infoD0] = GaugeSolveKii(Kii,Kia,opts,infoD0);

    otherwise
        error('SolveD0:UnknownMethod','Unknown D0 solver method.');
end

%% Keep the sign convention of the original code
D0 = -X;

if opts.verbose == 1
    fprintf('D0 solver: %s, relative residual: %.3e\n',...
        infoD0.methodUsed,infoD0.residualRelative);
end
end

function [X,infoD0] = AutoSolveKii(Kii,Kia,opts,infoD0)
%% Automatic solution strategy for Kii*X = Kia
% The automatic mode keeps the original formulation but changes the numerical
% solution strategy. Fast solvers are attempted first. If they fail or their
% residuals are larger than the tolerance, the code falls back to pinv.
%
% Default order:
% 1. Direct sparse solution, only when Kii is not close to singular;
% 2. Minimum-norm least-squares solution, if lsqminnorm is available;
% 3. Single-precision pseudo-inverse, if opts.useSingleInAuto = 1;
% 4. Double-precision pseudo-inverse as the final fallback.

X = [];

%% Try direct sparse solver
if isfinite(infoD0.rcondEst) && infoD0.rcondEst > opts.rcondTol
    [xTrial,infoTrial] = DirectSolveKii(Kii,Kia,opts,infoD0);
    infoD0 = AddTryHistory(infoD0,'direct',infoTrial);
    if infoTrial.flag == 1
        X = xTrial;
        infoD0 = infoTrial;
        infoD0.methodUsed = 'auto-direct';
        return
    end
else
    infoD0.tryHistory{end+1} = 'direct: skipped because Kii is singular or ill-conditioned';
end

%% Try minimum-norm least-squares solver
if exist('lsqminnorm','file') == 2
    [xTrial,infoTrial] = LsqminnormSolveKii(Kii,Kia,opts,infoD0);
    infoD0 = AddTryHistory(infoD0,'lsqminnorm',infoTrial);
    if infoTrial.flag == 1
        X = xTrial;
        infoD0 = infoTrial;
        infoD0.methodUsed = 'auto-lsqminnorm';
        return
    end
else
    infoD0.tryHistory{end+1} = 'lsqminnorm: skipped because this function is unavailable';
end

%% Try single-precision pseudo-inverse solver
if opts.useSingleInAuto == 1
    [xTrial,infoTrial] = SinglePinvSolveKii(Kii,Kia,opts,infoD0);
    infoD0 = AddTryHistory(infoD0,'pinv-single',infoTrial);
    if infoTrial.flag == 1
        X = xTrial;
        infoD0 = infoTrial;
        infoD0.methodUsed = 'auto-pinv-single';
        return
    end
else
    infoD0.tryHistory{end+1} = 'pinv-single: skipped by opts.useSingleInAuto';
end

%% Use double-precision pseudo-inverse as final fallback
[X,infoD0] = PinvSolveKii(Kii,Kia,opts,infoD0);
infoD0.methodUsed = 'auto-pinv';
infoD0 = AddTryHistory(infoD0,'pinv',infoD0);
end

function [X,infoD0] = DirectSolveKii(Kii,Kia,opts,infoD0)
%% Direct sparse solution
% This method is efficient when Kii is nonsingular and well-conditioned.
try
    X = Kii\Kia;
    infoD0.methodUsed = 'direct';
    [flag,residualRelative] = CheckSolveResidual(Kii,Kia,X,opts.resTol);
    infoD0.flag = flag;
    infoD0.residualRelative = residualRelative;
    if flag == 0
        infoD0.message = 'The residual of direct solver is larger than the tolerance.';
    end
catch ME
    X = [];
    infoD0.flag = 0;
    infoD0.residualRelative = Inf;
    infoD0.message = ME.message;
end
end

function [X,infoD0] = LsqminnormSolveKii(Kii,Kia,opts,infoD0)
%% Minimum-norm least-squares solution
% This method is useful for singular or semi-definite Kii in MATLAB.
try
    if exist('lsqminnorm','file') == 2
        X = lsqminnorm(full(Kii),full(Kia));
        infoD0.methodUsed = 'lsqminnorm';
    else
        X = [];
        infoD0.flag = 0;
        infoD0.residualRelative = Inf;
        infoD0.message = 'lsqminnorm is not available.';
        return
    end
    [flag,residualRelative] = CheckSolveResidual(Kii,Kia,X,opts.resTol);
    infoD0.flag = flag;
    infoD0.residualRelative = residualRelative;
    if flag == 0
        infoD0.message = 'The residual of lsqminnorm solver is larger than the tolerance.';
    end
catch ME
    X = [];
    infoD0.flag = 0;
    infoD0.residualRelative = Inf;
    infoD0.message = ME.message;
end
end

function [X,infoD0] = SinglePinvSolveKii(Kii,Kia,opts,infoD0)
%% Single-precision pseudo-inverse solution
% This method may accelerate the computation for some large problems.
% It is accepted only when the residual satisfies opts.resTol.
try
    X = double(pinv(single(full(Kii)))*single(full(Kia)));
    infoD0.methodUsed = 'pinv-single';
    [flag,residualRelative] = CheckSolveResidual(Kii,Kia,X,opts.resTol);
    infoD0.flag = flag;
    infoD0.residualRelative = residualRelative;
    if flag == 0
        infoD0.message = 'The residual of single-precision pinv solver is larger than the tolerance.';
    end
catch ME
    X = [];
    infoD0.flag = 0;
    infoD0.residualRelative = Inf;
    infoD0.message = ME.message;
end
end

function [X,infoD0] = PinvSolveKii(Kii,Kia,opts,infoD0)
%% Moore-Penrose pseudo-inverse solution
% This method is the robust reference solution corresponding to Eq.(9).
try
    X = pinv(full(Kii))*full(Kia);
    infoD0.methodUsed = 'pinv';
    [flag,residualRelative] = CheckSolveResidual(Kii,Kia,X,opts.resTol);
    infoD0.flag = flag;
    infoD0.residualRelative = residualRelative;
    if flag == 0
        infoD0.message = 'The residual of pinv solver is larger than the tolerance.';
    end
catch ME
    X = [];
    infoD0.flag = 0;
    infoD0.residualRelative = Inf;
    infoD0.message = ME.message;
end
end

function [X,infoD0] = GaugeSolveKii(Kii,Kia,opts,infoD0)
%% Gauge-fixed solution
% This method is fast but should be validated against pinv before use.
% The selected fixed DOFs should only remove rigid-body/gauge modes.
try
    nodeNum = size(Kii,1);
    if isfield(opts,'fixedDofs') && ~isempty(opts.fixedDofs)
        fixedDofs = opts.fixedDofs(:)';
    else
        if nodeNum >= 6
            fixedDofs = 1:6;
        elseif nodeNum >= 3
            fixedDofs = 1:3;
        else
            fixedDofs = 1;
        end
    end
    fixedDofs = unique(fixedDofs(fixedDofs>=1 & fixedDofs<=nodeNum));
    freeDofs = setdiff(1:nodeNum,fixedDofs);

    X = zeros(nodeNum,size(Kia,2));
    X(freeDofs,:) = Kii(freeDofs,freeDofs)\Kia(freeDofs,:);

    infoD0.methodUsed = 'gauge';
    infoD0.fixedDofs = fixedDofs;
    infoD0.rcondEst = EstimateRcondKii(Kii(freeDofs,freeDofs));
    [flag,residualRelative] = CheckSolveResidual(Kii,Kia,X,opts.resTol);
    infoD0.flag = flag;
    infoD0.residualRelative = residualRelative;
    if flag == 0
        infoD0.message = 'The residual of gauge-fixed solver is larger than the tolerance.';
    end
catch ME
    X = [];
    infoD0.flag = 0;
    infoD0.residualRelative = Inf;
    infoD0.message = ME.message;
end
end

function rcondEst = EstimateRcondKii(Kii)
%% Estimate reciprocal condition number
try
    if size(Kii,1) <= 3000
        rcondEst = rcond(full(Kii));
    else
        rcondEst = 1/condest(Kii);
    end
catch
    rcondEst = NaN;
end
end

function [flag,residualRelative] = CheckSolveResidual(Kii,Kia,X,resTol)
%% Check the relative residual of Kii*X = Kia
if isempty(X) || any(isnan(X(:))) || any(isinf(X(:)))
    flag = 0;
    residualRelative = Inf;
    return
end
residualRelative = norm(full(Kii*X-Kia),'fro')/max(1,norm(full(Kia),'fro'));
if isfinite(residualRelative) && residualRelative <= resTol
    flag = 1;
else
    flag = 0;
end
end

function infoD0 = AddTryHistory(infoD0,methodName,infoTrial)
%% Record each attempted solver in automatic mode
if isempty(infoTrial.message)
    messageText = '';
else
    messageText = ['; ',infoTrial.message];
end
infoD0.tryHistory{end+1} = [methodName,': residual = ',...
    num2str(infoTrial.residualRelative),messageText];
end
