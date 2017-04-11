function [beta,tstat]=rols(y,x,cons);

% if cons==1 it adds a constant term to x
   Temp=find(isnan(y)==1);
    if ~isempty(Temp)
        y(Temp)=[];
        x(Temp,:)=[];
    end;
    
    for i=1:min(size(x))
        Temp=find(isnan(x(:,i))==1);
        if ~isempty(Temp)
        y(Temp)=[];
        x(Temp,:)=[];
        end;
    end;
if cons==1
    x=[ones(length(x),1) x];
end;

%length(x)
beta=(x'*x)^-1*(x'*y);
TX=((y-x*beta)*ones(1,min(size(x)))).*x;
        RobustVar=(x'*x)^-1*(TX'*TX)*(x'*x)^-1;
        tstat=beta./sqrt(diag(RobustVar));
        
        