% Code from email on 17/07/2018
% Nothing changed except display stuff

clear variables;
close all;
% rng(0);
parpool('local',3);

% Initialize RNG seed (specific to par-loop)
spmd
    rng(0,'combRecursive');
end

%n1 = 100;
%n0 = 1000;
n1 = 10;
n0 = 100;
k = 10;
a = 0.10; b = 0.90;
r = 4;
h = 0.10;

m = 0; 
for t=0:r
    m = m + nchoosek(r,t)*(a^(r-t))*((b-a)^t)/(t+1);
end
m2 = 0; 
for t=0:2*r
    m2 = m2 + nchoosek(2*r,t)*(a^(2*r-t))*((b-a)^t)/(t+1);
end

v = sqrt(k*(m2-m^2));

deltaLambda = 0.01;
firstLambda = 0:deltaLambda:1;
maxLambda = 20;

%options = optimoptions('quadprog','StepTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',2000,'Display','off');
options = optimoptions('quadprog','TolFun',1e-10,'TolX',1e-10,'MaxIter',2000,'Display','off');

T = 100;
Estp = zeros(T,1); Estnp = zeros(T,1); Estm = zeros(T,1); Estmopt = zeros(T,1);
MSEp = zeros(T,1); MSEnp = zeros(T,1); MSEm = zeros(T,1); MSEmopt = zeros(T,1);
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
    
    x1 = a+(b-a)*rand(k,n1);
    x0 = (a-h)+(b-a+2*h)*rand(k,n0);
    x0 = sqrt(x0);
    y1 = sum(x1.^r,1)'/v+randn(n1,1);
    y0 = sum(x0.^r,1)'/v+randn(n0,1);

    H = 2*(x0'*x0);
    
    minMSE = inf;
    Wp = zeros(n0,n1);
    optlambda = 0;
    Lambda = firstLambda;
    
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
    
    matches = zeros(M,n1);
    for i=1:n1
         x = x1(:,i);
         D = x0 - kron(ones(1,n0),x);
         delta = diag(D'*D); 
         [sorted,I]=sort(delta);
         matches(:,i) = I(1:M);
    end
    
    minMSEm = inf;
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
    
    Wm = zeros(n0,n1);
    Wm(sub2ind(dim,matches(1,:),ones(size(matches(1,:),1),1)*(1:dim(2))))=1;
    
    
    y1 = sum(x1.^r,1)'/v+randn(n1,1);
    y0 = sum(x0.^r,1)'/v+randn(n0,1);
    
    Wnp = zeros(n0,n1);
    for i=1:n1
        x = x1(:,i);
        D = x0 - kron(ones(1,n0),x);
        delta = diag(D'*D); 
        w = quadprog(H,-2*x0'*x,[],[],ones(1,n0),1,zeros(n0,1),ones(n0,1),[],options);
        Wnp(:,i) = w;
    end
    
    Estp(t) = mean(y1-Wp'*y0);
    Estnp(t) = mean(y1-Wnp'*y0);
    Estm(t) = mean(y1-Wm'*y0);
    Estmopt(t) = mean(y1-Wmopt'*y0);
    MSEp(t) = (y1-Wp'*y0)'*(y1-Wp'*y0)/n1;
    MSEnp(t) = (y1-Wnp'*y0)'*(y1-Wnp'*y0)/n1;
    MSEm(t) = (y1-Wm'*y0)'*(y1-Wm'*y0)/n1;
    MSEmopt(t) = (y1-Wmopt'*y0)'*(y1-Wmopt'*y0)/n1;  
    Densp(t) = mean(sum(Wp>dthr,1));
    Densnp(t) = mean(sum(Wnp>dthr,1));
    maxminDensp(t,:) = [min(sum(Wp>dthr)) max(sum(Wp>dthr))];
    maxminDensnp(t,:) = [min(sum(Wnp>dthr)) max(sum(Wnp>dthr))];
    
    toc
end

filename = sprintf('n1_%d_n0_%d_k_%d_r_%d_%d%d_%d_T_%d_nested',n1,n0,k,r,100*a,100*b,100*h,T);
save(filename,'n1','n0','k','r','a','b','h','maxLambda','M','MSEp','MSEnp','MSEm','MSEmopt','Estp','Estnp','Estm','Estmopt','Densp','Densnp','maxminDensp','maxminDensnp','lambdavalues','mvalues','T');

sprintf('RMSE individual effects')

sqrt(mean([MSEp MSEnp MSEm MSEmopt]))

sprintf('RMSE average effects')

sqrt(mean([Estp Estnp Estm Estmopt].^2))

sprintf('|Bias|')

abs(mean([Estp Estnp Estm Estmopt]))

% Print to file and screen
Name = {'PenSynth';'NoPenSynth';'Matching';'OptMatching'};
RMSEindiv = sqrt(mean([MSEp MSEnp MSEm MSEmopt]))';
RMSEatt = sqrt(mean([Estp Estnp Estm Estmopt].^2))';
Bias = abs(mean([Estp Estnp Estm Estmopt]))';
Results = table(num2str(RMSEindiv,'%.4f'),num2str(RMSEatt,'%.4f'),num2str(Bias,'%.4f'),'RowNames',Name);
Results.Properties.VariableNames = {'RMSEindiv' 'RMSEatt' 'Bias'}

delete(gcp('nocreate')); 