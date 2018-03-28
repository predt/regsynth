clear variables;
rng(0);

n1 = 100;
n0 = 1000;
k = 4;
a = 0.10; b = 0.90;
r = 8;


x1 = a+(b-a)*rand(k,n1);
x0 = rand(k,n0);

m = 0; 
for t=0:r
    m = m + nchoosek(r,t)*(a^(r-t))*((b-a)^t)/(t+1);
end
m2 = 0; 
for t=0:2*r
    m2 = m2 + nchoosek(2*r,t)*(a^(2*r-t))*((b-a)^t)/(t+1);
end

v = sqrt(k*(m2-m^2));
y1 = sum(x1.^r,1)'/v+randn(n1,1);
y0 = sum(x0.^r,1)'/v+randn(n0,1);

H = 2*(x0'*x0);

options = optimoptions('quadprog','StepTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',2000,'Display','off');

Lambda = 0.01:0.01:1;
MSE = zeros(length(Lambda),1);
for t = 1:length(Lambda)
    sprintf('First loop: %d',t)
    lambda = Lambda(t);
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
    MSE(t) = mse;
end

plot(Lambda,MSE);
[minimum,I]=min(MSE);
lambda = Lambda(I)

T = 1000;
MSEp = zeros(T,1); MSEnp = zeros(T,1); MSEm = zeros(T,1);
for t =1:T
    sprintf('Second loop: %d',t)
    x1 = a+(b-a)*rand(k,n1);
    x0 = rand(k,n0);
    y1 = sum(x1.^r,1)'/v+randn(n1,1);
    y0 = sum(x0.^r,1)'/v+randn(n0,1);

    H = 2*(x0'*x0);
    
    Wp = zeros(n0,n1);
    Wnp = zeros(n0,n1);
    Wm =  zeros(n0,n1);
    for i=1:n1
        x = x1(:,i);
        D = x0 - kron(ones(1,n0),x);
        delta = diag(D'*D); 
        f = lambda*delta-2*x0'*x;
        w = quadprog(H,f,[],[],ones(1,n0),1,zeros(n0,1),ones(n0,1),[],options);
        Wp(:,i) = w;
        w = quadprog(H,-2*x0'*x,[],[],ones(1,n0),1,zeros(n0,1),ones(n0,1),[],options);
        Wnp(:,i) = w;
        [minimum,I]=min(delta);
        Wm(I,i)=1;
    end
    MSEp(t) = (y1-Wp'*y0)'*(y1-Wp'*y0);
    MSEnp(t) = (y1-Wnp'*y0)'*(y1-Wnp'*y0);
    MSEm(t) = (y1-Wm'*y0)'*(y1-Wm'*y0);
end

save('n1_100_n0_1000_k_4_r_8_1090.mat','n1','n0','k','r','a','b','MSEp','MSEnp','MSEm','Lambda','MSE');

%sum(W>.0000000001,1)

%vsample=sqrt(var(sum(x1.^r,2)))
