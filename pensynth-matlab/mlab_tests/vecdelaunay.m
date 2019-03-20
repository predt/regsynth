clear all;
rng(0);
k = 5;
n1 = 10;
n0 = 100;
X0 = rand(n0,k);

DT = delaunayn(X0); DTp = DT';
[rDT,cDT] = size(DT);

A = [X0(DTp(:),:) ones(numel(DT),1)]';
B = mat2cell(sparse(A),[k+1],cDT*[ones(rDT,1)]);
C = sparse(blkdiag(B{:}));


X1 = rand(k,n1);
b=C\kron(ones(rDT,1),[X1;ones(1,n1)]);

b_ordered=reshape(b(:),k+1,[]);
sign_sum = reshape(sum(sign(b_ordered)),rDT,[]);
[r,c] = find(sign_sum==k+1);
triangle = zeros(n1,1);
triangle(c)=r; 



