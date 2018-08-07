%%% Compute Delaunay triangulation and select triangle containing X_1
%%% Jérémy L'Hour
%%% 24/07/2018

clear variables;
close all;
rng(0);

k = 3; % dimension
n1 = 1;
n0 = 100;
X0 = randn(n0,k);
X1 = randn(k,n1);
DT = delaunayn(X0);

minu = -Inf;
t = 0;

while minu < 0 % keep going until X1 is in the simplex
	t = t+1;
    if t >= size(DT,1)
        disp('X1 outside convex hull.');
        break
    end
	P = X0(DT(t,:),:); % k+1 x k matrix of vertices of t
	u = inv([ones(1,k+1);P']) * [1;X1'];
	minu = min(u);
end

% notice u is the solution when lambda is close to zero.

% Alberto's solution (works faster)

DT = delaunayn(X0); DTp = DT';
[rDT,cDT] = size(DT);

A = [X0(DTp(:),:) ones(numel(DT),1)]';
B = mat2cell(sparse(A),[k+1],cDT*[ones(rDT,1)]);
C = sparse(blkdiag(B{:}));

b=C\kron(ones(rDT,1),[X1;ones(1,n1)]);
b_ordered=reshape(b(:),k+1,[]);

sign_sum = reshape(sum(sign(b_ordered)),rDT,[]);
[r,c] = find(sign_sum==k+1);
triangle = zeros(n1,1); triangle(c) = r;  % a zero for units outside convex hull
W = zeros(n1,k+1); W(c,:) = b_ordered(:,sign_sum==k+1)';

Wpure = zeros(n0,n1);
for i = 1:n1
	if(triangle(i)>0)
		idx = DT(triangle(i),:);
		Wpure(idx,i) = W(i,:);
	elseif(triangle(i)==0)
		Wpure(:,i) = Wnp(:,i);
	end
end