% MODIFIED BY JEREMY ON 23 JULY 18

% Formula for r real
% bias correction implemented
% lambda minimizes MSE now
% RNG initilaized for each worker - results are reproducible

clear variables;
close all;
parpool('reducedform',100);
% Initialize RNG seed (specific to par-loop)
spmd
    rng(0,'combRecursive');
end

% Parameters

n1 = 10;
n0 = 100;
k = 4;
a = 0.10; b = 0.90;
r = 1;
h = 0.10;

v = sqrt(k*((b^(2*r+1)-a^(2*r+1))/((b-a)*(2*r+1)) - (b^(r+1)-a^(r+1))^2/((b-a)*(r+1))^2));

deltaLambda = 0.01;
firstLambda = 0:deltaLambda:1;
maxLambda = 20;

options = optimoptions('quadprog','StepTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',2000,'Display','off');

T = 1000;
Estp = zeros(T,1); Estnp = zeros(T,1); Estm = zeros(T,1); Estmopt = zeros(T,1);
MSEp = zeros(T,1); MSEnp = zeros(T,1); MSEm = zeros(T,1); MSEmopt = zeros(T,1);
Estp_bc = zeros(T,1); Estnp_bc = zeros(T,1); Estm_bc = zeros(T,1); Estmopt_bc = zeros(T,1);
MSEp_bc = zeros(T,1); MSEnp_bc = zeros(T,1); MSEm_bc = zeros(T,1); MSEmopt_bc = zeros(T,1);
Densp = zeros(T,1); Densnp = zeros(T,1); 
maxminDensp = zeros(T,2); maxminDensnp = zeros(T,2); 
lambdavalues = zeros(T,1);
M = 20;
dthr = 0.001;
margin = 25;
dlambda = margin*deltaLambda;
mvalues = zeros(T,1);

parfor t =1:T
    
    tic
    sprintf('Iteration: %d',t)
    
    stream = RandStream.getGlobalStream();
    stream.Substream = t;

    % 0. Simulate Data
    x1 = a+(b-a)*rand(k,n1);
    x0 = (a-h)+(b-a+2*h)*rand(k,n0);
    y1 = sum(x1.^r,1)'/v+randn(n1,1);
    y0 = sum(x0.^r,1)'/v+randn(n0,1);

    H = 2*(x0'*x0);
    
    minMSE = Inf;
    Wp = zeros(n0,n1);
    optlambda = 0;
    Lambda = firstLambda;

    % 1. pen synth w/ optimized lambda
    % here lambda is set to optimize bias
    
    j = 1;
    while j <= length(Lambda)
        lambda = Lambda(j);
        W = zeros(n0,n1);
        for i=1:n1
            x = x1(:,i);
            D = x0 - kron(ones(1,n0),x);
            delta = diag(D'*D); 
            f = lambda*delta-2*x0'*x;
            w = quadprog(H,f,[],[],ones(1,n0),1,zeros(n0,1),ones(n0,1),[],options);
            W(:,i) = w;
        end
        mse = (y1-W'*y0)'*(y1-W'*y0);
        if mse < minMSE
            minMSE = mse;
            optlambda = lambda;
            Wp = W;
        end
        j = j + 1;
        if j > length(Lambda)         
             if   ((round(lambda-optlambda,5)<round(dlambda,5)) && (round(lambda+deltaLambda,5)<round(maxLambda,5)))
                Lambda = [Lambda lambda+(deltaLambda:deltaLambda:dlambda)];
            end
        end
    end
 
    lambdavalues(t) = optlambda;

    % 2. matching
    
    % 2.1 collects the indices of the m closet points in 'matches'
    matches = zeros(M,n1);
    for i=1:n1
         x = x1(:,i);
         D = x0 - kron(ones(1,n0),x);
         delta = diag(D'*D); 
         [sorted,I]=sort(delta);
         matches(:,i) = I(1:M);
    end
    
    % 2.2 computes the corresponding MSE
    minMSEm = Inf;
    Wmopt = zeros(n0,n1);
    optm = 0;
    dim = size(Wmopt); 
    for m = 1:M
        W = zeros(n0,n1);
        W(sub2ind(dim,matches(1:m,:),ones(size(matches(1:m,:),1),1)*(1:dim(2))))=1/m;
        mse = (y1-W'*y0)'*(y1-W'*y0);
        if mse < minMSEm
            minMSEm = mse;
            optm = m;
            Wmopt = W;
        end
    end
    
    mvalues(t) = optm;

    % 3. non pen synth
    Wnp = zeros(n0,n1);
    for i=1:n1
        x = x1(:,i);
        D = x0 - kron(ones(1,n0),x);
        delta = diag(D'*D); 
        w = quadprog(H,-2*x0'*x,[],[],ones(1,n0),1,zeros(n0,1),ones(n0,1),[],options);
        Wnp(:,i) = w;
    end
    

    % 4. one-to-one matching
    Wm = zeros(n0,n1);
    Wm(sub2ind(dim,matches(1,:),ones(size(matches(1,:),1),1)*(1:dim(2))))=1;
    
    % 5. Evaluate performance on subsequent period

    % simulate outcome for second period
    y1 = sum(x1.^r,1)'/v+randn(n1,1);
    y0 = sum(x0.^r,1)'/v+randn(n0,1);
    
    % Bias
    Estp(t) = mean(y1-Wp'*y0);
    Estnp(t) = mean(y1-Wnp'*y0);
    Estm(t) = mean(y1-Wm'*y0);
    Estmopt(t) = mean(y1-Wmopt'*y0);

    % MSE
    MSEp(t) = (y1-Wp'*y0)'*(y1-Wp'*y0)/n1;
    MSEnp(t) = (y1-Wnp'*y0)'*(y1-Wnp'*y0)/n1;
    MSEm(t) = (y1-Wm'*y0)'*(y1-Wm'*y0)/n1;
    MSEmopt(t) = (y1-Wmopt'*y0)'*(y1-Wmopt'*y0)/n1; 

    % Mean sparsity index 
    Densp(t) = mean(sum(Wp>dthr,1));
    Densnp(t) = mean(sum(Wnp>dthr,1));

    % Min and max of sparsity indices
    maxminDensp(t,:) = [min(sum(Wp>dthr)) max(sum(Wp>dthr))];
    maxminDensnp(t,:) = [min(sum(Wnp>dthr)) max(sum(Wnp>dthr))];

    % 6. Bias correction
    f0 = [x0; x0.^2];
    f1 = [x1; x1.^2];
    mu0 = inv(f0*f0')*(f0*y0);

    mu_hat0 = f0'*mu0;
    mu_hat1 = f1'*mu0;
    
    % bc: Bias
    Estp_bc(t) = mean(y1-Wp'*y0 - (mu_hat1-Wp'*mu_hat0));
    Estnp_bc(t) = mean(y1-Wnp'*y0 - (mu_hat1-Wnp'*mu_hat0));
    Estm_bc(t) = mean(y1-Wm'*y0 - (mu_hat1-Wm'*mu_hat0));
    Estmopt_bc(t) = mean(y1-Wmopt'*y0 - (mu_hat1-Wmopt'*mu_hat0));

    % bc: MSE
    MSEp_bc(t) = (y1-Wp'*y0 - (mu_hat1-Wp'*mu_hat0))'*(y1-Wp'*y0 - (mu_hat1-Wp'*mu_hat0))/n1;
    MSEnp_bc(t) = (y1-Wnp'*y0 - (mu_hat1-Wnp'*mu_hat0))'*(y1-Wnp'*y0 - (mu_hat1-Wnp'*mu_hat0))/n1;
    MSEm_bc(t) = (y1-Wm'*y0 - (mu_hat1-Wm'*mu_hat0))'*(y1-Wm'*y0 - (mu_hat1-Wm'*mu_hat0))/n1;
    MSEmopt_bc(t) = (y1-Wmopt'*y0 - (mu_hat1-Wmopt'*mu_hat0))'*(y1-Wmopt'*y0 - (mu_hat1-Wmopt'*mu_hat0))/n1; 
    
    toc
end

% End of loop / saving results
filename = sprintf('n1_%d_n0_%d_k_%d_r_%d_%d%d_%d_T_%d_nested',n1,n0,k,r,100*a,100*b,100*h,T);
save(filename,'n1','n0','k','r','a','b','h','maxLambda','M','MSEp','MSEnp','MSEm','MSEmopt','Estp','Estnp','Estm','Estmopt','Densp','Densnp','MSEp_bc','MSEnp_bc','MSEm_bc','MSEmopt_bc','Estp_bc','Estnp_bc','Estm_bc','Estmopt_bc','maxminDensp','maxminDensnp','lambdavalues','mvalues','T');

% Print to file and screen
Name = {'PenSynth';'NoPenSynth';'Matching';'OptMatching';'PenSynth_bc';'Synth_bc';'Matching_bc';'OptMatching_bc'};
RMSEindiv = sqrt(mean([MSEp MSEnp MSEm MSEmopt MSEp_bc MSEnp_bc MSEm_bc MSEmopt_bc]))';
RMSEatt = sqrt(mean([Estp Estnp Estm Estmopt Estp_bc Estnp_bc Estm_bc Estmopt_bc].^2))';
Bias = abs(mean([Estp Estnp Estm Estmopt Estp_bc Estnp_bc Estm_bc Estmopt_bc]))';
Results = table(num2str(RMSEindiv,'%.4f'),num2str(RMSEatt,'%.4f'),num2str(Bias,'%.4f'),'RowNames',Name);
Results.Properties.VariableNames = {'RMSEindiv' 'RMSEatt' 'Bias'}

txtname = sprintf('n1_%d_n0_%d_k_%d_r_%d_%d%d_%d_T_%d.txt',n1,n0,k,r,100*a,100*b,100*h,T);
writetable(Results,txtname,'Delimiter','\t','WriteRowNames',true);

delete(gcp('nocreate')); 
