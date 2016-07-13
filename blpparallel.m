%
% This is an extremely simple implementation of the traditional BLP
% algorithm.  It uses the GMM estimator and the contraction mapping.
% The major difference is that this algorithm solves contraction mapping market
% by market and utilizes a separate thread for each market
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
% at the cost of some memory. (Or in a MEX file). This mostly helps when NS
% is big.
% 3. You can save some time by defining delta in the main function. Do not
% do this, it can induce chattering. Most of the time is spent computing
% the gradient anyway.
% 4. The gradient loops over theta_2 variables.  This can be done with
% multi-dimensional arrays more quickly but that is a) even harder to
% understand and b) more memory intensive.
% 5. Dataset Z should contain all instruments (not just excluded IV) but X
% variables as well.
% It expects a special dataset structure.  To learn how to define the structure visit:
%
% http://www.columbia.edu/~cc3264

function [t,f]=blpparallel(data,draws)
    [ns K] = size(draws.nu);
    [T] = max(data.mktid);
    [J] = max(data.prodid);
    n = length(data.sjt);
    subs = [repmat(data.mktid,[ns 1]) reshape(repmat(1:ns,[J*T 1]),J*T*ns,1)  repmat(data.prodid,[ns 1])];
    data.Z=data.IV;
    
    % precompute these to avoid passing the entire structure to each node--
    % this doubles memory usage but eliminates data passing
    for t=1:T,
        mktsubs{t} =  find(data.mktid==t);
        market{t} = data(data.mktid==t,:);
    end
    
    % Do not define delta outside the loop-- it creates chattering    
    % GMM Step 1
    theta2=0.5*ones(K,1);
    A = (data.Z'*data.Z) \ eye(size(data.Z,2)); 
    op = optimset('Display','iter', 'GradObj','on','Algorithm','interior-point','TolFun',1e-8,'Hessian','bfgs','SubproblemAlgorithm','cg');
    tic
    [t,f]=fmincon(@gmmobjective,theta2,-eye(K),zeros(K,1),[],[],[],[],[],op);
    toc
    theta1
    
    function [f,grad]=gmmobjective(theta2)
        delta=[];
        % begin parallel block
        parfor t=1:max(data.mktid),
            
            mkt=market{t};
            iter=0;
            deltas=zeros(J,1);
            deltaold=ones(size(deltas));
            % compute exp(mu) once
            emu=exp(mkt.x*(draws.nu.*repmat(theta2',[ns 1]))');
            % iterate contraction mapping to convergence
            while(iter < 1500 && max(abs(deltas-deltaold)) > 2e-14),
                iter = iter+1; deltaold=deltas;
                utils=repmat(exp(deltas),[1,ns]).*emu;
                pijt=utils./(1+repmat(sum(utils),[J 1]));
                shares = pijt*draws.w;
                deltas= deltas+log(mkt.sjt)-log(shares);
            end
            delta=[delta; deltas];
        end              
        [theta1,xi]=ivregression(delta,data.x,data.Z,A);
        f=xi'*data.Z * A * data.Z'*xi;
        % Derivative Calculation
        if nargout > 1,
            grad = computegradient(delta,theta2,xi);
        end
    end
    
    % this computes the shares: For NS large it helps to precompute X_jt * V_ik
    function [pijt]=share(delta,sigma)
        v=draws.nu .* repmat(sigma',[ns 1]);
        utils = repmat(exp(delta),[1 ns]).* exp(data.x * v');
        denom = 1+accumarray(subs(:,1:2),utils(:));
        pijt = utils./denom(data.mktid,:);
    end

    function [g] = computegradient(delta, theta2, xi)
        pijt = share(delta,theta2);
        % compute the derivative of s_jt wrt each sigma in parallel
        parfor k=1:K,
            x1=repmat(data.x(:,k),[1 ns]).* pijt;
            x2=accumarray(subs(:,1:2),x1(:));
            x3=repmat(data.x(:,k),[1 ns]) - x2(data.mktid,:);
            dsigma(:,k)=(pijt.*x3)*(draws.nu(:,k).*draws.w);
        end
        % 3-D version of choice probabilities
        pmat = accumarray(subs(:,[ 2 3 1]), pijt(:));
        % Compute the Jacobian -- invert market by market
        Jacobian=[];
        parfor tt=1:T,
            JacOut=(diag(draws.w'*pmat(:,:,tt))-reshape(outerprodall(pmat(:,:,tt))*draws.w,[J J])) \ dsigma(mktsubs{tt},:);
            Jacobian=[Jacobian; JacOut];
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