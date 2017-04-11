function [beta,tstat]=fgls(y,x,cons,VarCov);

% if cons==1 it adds a constant term to x

   Temp=find(isnan(y)==1);
    if ~isempty(Temp)
        y(Temp)=[];
        x(Temp,:)=[];
        VarCov(Temp,:)=[];        VarCov(:,Temp)=[];
    end;
    
    for i=1:min(size(x))
        Temp=find(isnan(x(:,i))==1);
        if ~isempty(Temp)
        y(Temp)=[];
        x(Temp,:)=[];
        VarCov(Temp,:)=[];        VarCov(:,Temp)=[];
        end;
    end;
if cons==1
    x=[x ones(length(x),1)];
end;

%length(x)
beta=(x'*x)^-1*(x'*y);
SE=sqrt(diag(inv(x'*x)*(x'*VarCov*x)*inv(x'*x)));
        tstat=beta./SE;
        
        