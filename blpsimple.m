%
% This is an extremely simple implementation of the traditional BLP
% algorithm.  It uses the GMM estimator and the contraction mapping.
%
% This code is written with the goal of being a clear single file
% implementation of the algorithm. It works best with a 2008 or later
% version of MATLAB and makes some use of new/advanced features in Matlab.
% It uses nested functions and structures. It also involves heavy use
% of accumarray.
%
% This is NOT prouction code.
% 1. It does no input checking at all.
% 1a. It doesn't provide much output either.
% 2. You can save a lot of time by precomputing Z_ijt^k = v_ik * x_jt^k but
% at the cost of some memory. (Or in a MEX file).
% 3. You can save some time by defining delta in the main function rather
% than the gmmobjective() function.  The downside is that this may induce
% additional numerical error especially if finite difference derivatives
% are used.
% 4. The gradient loops over theta_2 variables.  This can be done with
% multi-dimensional arrays more quickly but that is a) even harder to
% understand and b) more memory intensive.
% 5. Dataset Z should contain all instruments (not just excluded IV) but X
% variables as well

% It expects a special dataset structure.  To learn how to define the structure visit:
%
% http://www.pantheon.yale.edu/~ctc25/

function [t,f]=blpsimple(data)
    [ns K] = size(data.nu);
    [T] = max(data.subs(:,1));
    [J] = max(data.subs(:,3));
    n = length(data.sjt);
    [ns k] = size(data.nu);
%    data.subs = [repmat(data.mktid,[ns 1]) reshape(repmat(1:ns,[J*T 1]),J*T*ns,1)  repmat(data.prodid,[ns 1])];

    nux=permute(permute(repmat(data.nu,[1 1 n]),[3 2 1]).* repmat(data.x,[1 1 ns]),[1 3 2]);
    
    theta2=ones(K,1);
    % you can define delta here and not in the iter=0 line to reuse last
    % iterations guess-- this speeds things up 2-3x but creates chatter
    %delta = zeros(n,1);
    s0 = accumarray([data.mktid],data.sjt);
    deltastart = data.x*regress(log(data.sjt)-log(s0(data.mktid)),data.x);
    
    % GMM Step 1
    A = (data.Z'*data.Z) \ eye(size(data.Z,2)); 
    op = optimset('Display','iter', 'GradObj','on','DerivativeCheck','on','Algorithm','interior-point','TolFun',1e-8);
tic
    [t,f]=fminunc(@gmmobjective,theta2,op);
 toc
    % GMM Step 2 -- don't reweight if objective near zero!
    if(rcond(data.Z'*xi * xi'*data.Z) < 1e-14)
        return
    else
        A = (data.Z'*xi * xi'*data.Z) \ eye(size(data.Z,2));
        [t,f]=fminunc(@gmmobjective,theta2,op);
    end
    
    function [f,grad]=gmmobjective(theta2)
        % initialize contraction mapping -- if delta_0 is not always the same
        % then we have noisy function evaluation. It should be defined here
        % or above but not in both places.
        iter=0;       deltaold = zeros(n,1); 
        delta=deltastart;
        % contraction mapping step--make sure the tolerance is small!
        while(iter < 1500 & max(abs(delta-deltaold)) > 2e-14),
            iter = iter+1; deltaold=delta;
           delta= delta+log(data.sjt)-log(share(delta,theta2)*data.w);
%             delta= delta.*(data.sjt./(share(delta,theta2)*data.w));
        end
        [theta1,xi]=ivregression((delta),data.x,data.Z,A);
        f=xi'*data.Z * A * data.Z'*xi;
%        disp(iter)
        % Derivative Calculation
        if nargout > 1,
            grad = computegradient((delta),theta2,xi);
        end
    end
    
    % this computes the shares-- there is a big speedup to precomputing
    % x_jtk * nu_ik
    function [pijt]=share(delta,sigma)
        v=data.nu .* repmat(sigma',[ns 1]);
        utils = repmat(exp(delta),[1 ns]).* exp(data.x * v');
        denom = 1+accumarray(data.subs(:,1:2),utils(:));
        pijt = utils./denom(data.mktid,:);
    end

    % this code can be written without loops an much faster-- but tougher on memory!
    function [g] = computegradient(delta, theta2, xi)
        pijt = share(delta,theta2);
        % compute the derivative of s_jt wrt each sigma
        for k=1:K,
            x1=repmat(data.x(:,k),[1 ns]).* pijt;
            x2=accumarray(data.subs(:,1:2),x1(:));
            x3=repmat(data.x(:,k),[1 ns]) - x2(data.mktid,:);
            dsigma(:,k)=(pijt.*x3)*(data.nu(:,k).*data.w);
        end
        % Compute the Jacobian -- invert market by market
        Jacobian = zeros(n,K);
        pmat = accumarray([data.subs(:,[2 3  1])], pijt(:));
        for t=1:T,
            subs = (t-1)*J+1:t*J;
            Jacobian(subs,:)=(diag(data.w'*pmat(:,:,t))-reshape(outerprodall(pmat(:,:,t))*data.w,[J J])) \  dsigma(subs,:);
        end
        g= -2*Jacobian' * data.Z * A *data.Z' *xi; 
    end

    function [beta,resid]=ivregression(Y,X,Z,W)
        if nargin <4
            W = (Z'*Z) \ eye(size(Z,2));
        end
        beta=(X'*Z * W * Z'*X)\(X'*Z * W * Z'*Y);
        if nargout >1
            resid=Y-X * beta;
        end
    end
end