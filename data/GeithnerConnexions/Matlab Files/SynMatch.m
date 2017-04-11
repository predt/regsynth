% SynMatch function estimates the coefficient on Geithner connection using
% Synthetic Matching method.
function [y,WeightMatrix]=SynMatch(Parameters,Returns,Treatment,Names)

alpha=Parameters.alpha;             % exclude firms with goodness of match alpha*treatment group firms....&&&
NumMatch=Parameters.NumMatch;       % Number of firms used for synthetic matching
NumIter=Parameters.NumIter;         %Number of iterations for the placebo test used for construction of H0 test confidence interval
EstWindow=Parameters.EstWindow;     % Estimation window used for synthetic matching
EvntDate=Parameters.EventDate;      %Date of actual event
PlcbTstWindow=Parameters.TstWindow;  %Placebo Test Window
EvntWindowSize=Parameters.EvntWindowSize;  %This can be 0, 1, ...0 means only event date, 1 means the event date and the date after...
NumReport=Parameters.NumReport;     % Number of firms with the highest weight that needs to be reported
NumMatch=min(sum(Treatment==0)-1,NumMatch);
NumReport=min(NumReport,NumMatch);
% Returns is a matrix that contains the return of all financial
% institutions.

%Treatment is a dummy that is equal to one if the firm is in treatment
%group
%Names is the name of firms

opt=optimset('Diagnostics','off','Display','off','LargeScale','off');
k=0;

for i=1:length(Treatment)   % In this loop we find a synthetic match (from firm in the control group) for each firm
IndCon=0;
if Treatment(i)==1
    IndCon=1;
    k=k+1;
end;

d=Returns(EstWindow,i);
C=[Returns(EstWindow,:) ];             Cx=[Returns];  e=Treatment;     f=1:length(Treatment);  n=Names;

C(:,i)=[];      Cx(:,i)=[];     e(i)=[];  f(i)=[];   n(i)=[];

InTreat=find(e>0);  %eliminating firms in the control group from the pull of potential matches
C(:,InTreat)=[];      Cx(:,InTreat)=[];       f(InTreat)=[];     n(InTreat)=[];

C(isnan(C)==1)=0;        % treating missing values as zeros
Cx(isnan(Cx)==1)=0;
d(isnan(d)==1)=0;       

PreSel=corr(d,C);                   % this is to reduce the dimension of the maximization problem...
PreSel(isnan(PreSel)==1)=0;
[x0 Ind1]=sort(PreSel,'descend');     

C=C(:,Ind1(1:NumMatch));
Cx=Cx(:,Ind1(1:NumMatch));  
x0=x0(1:NumMatch);         % initial guess for optimization weights

A=[];   b=[];   Aeq=ones(1,NumMatch);  beq=1;  lb=zeros(NumMatch,1);     ub=ones(NumMatch,1);
%x0=lb;  x0(2)=.1;    x0(3)=.9;
x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,opt);

SynthRet(:,i)=Cx*x;

if IndCon==1     %saving weights
    [weight Ind2]=sort(x,'descend');
    Weights(k,:)=weight(1:NumReport);
    SynMatchNames(k,:)=n(Ind1(Ind2(1:NumReport)));
end;
end;

AR=Returns-SynthRet;  [m1 n1]=size(AR);
CAR=AR(1:m1-EvntWindowSize-1,:);
for j=2:EvntWindowSize+1
    CAR=CAR+AR(j:m1-EvntWindowSize-2+j,:);
end;

% %%%%Computing SPE and eliminating bad matches 
SPE=sum(AR(EstWindow,:).^2);
SPETr=mean(SPE(Treatment>0))*alpha;
CAR(:,SPE>SPETr)=[];
Weights(SPE(Treatment>0)>SPETr,:)=[];
SynMatchNames(SPE(Treatment>0)>SPETr,:)=[];
Treatment(SPE>SPETr)=[];

[m1 n1]=size(CAR);

CAR(isnan(CAR)==1)=0;

Ind=find(Treatment>0); 
Weights2=(std(CAR(EstWindow,Ind))).^-1;

WMatrix=ones(m1,1)*Weights2;
    
EstCoef=sum(CAR(:,Ind).*WMatrix,2)/sum(Weights2);

% Constructing confidence interval using placebo test
 CAR(:,Ind)=[];     

[m1 n1]=size(CAR);
for i=1:NumIter
    Sample=randperm(n1);     Sample=Sample(1:k);
    Weights2=(std(CAR(EstWindow,Sample))).^-1;
%    Weights2=max(Weights2,10*median(Weights2));
    WMatrix=ones(m1,1)*Weights2;
    PlaceboCoef(:,i)=sum(CAR(:,Sample).*WMatrix,2)/sum(Weights2);
end;
    
    
for i=1:length(EstCoef)
    a1=sort((PlaceboCoef(i,:)),'descend');
    pvalue10H(i)=a1(floor(.05*NumIter));    pvalue10L(i)=a1(floor(.95*NumIter));
    pvalue5H(i)=a1(floor(.025*NumIter));     pvalue5L(i)=a1(floor(.975*NumIter));
    pvalue1H(i)=a1(floor(.005*NumIter));      pvalue1L(i)=a1(floor(.995*NumIter));
    
end;
 
a= length(find((EstCoef(PlcbTstWindow)>pvalue10H(PlcbTstWindow)')|(EstCoef(PlcbTstWindow)<pvalue10L(PlcbTstWindow)')));
b= length(find((EstCoef(PlcbTstWindow)>pvalue5H(PlcbTstWindow)')|(EstCoef(PlcbTstWindow)<pvalue5L(PlcbTstWindow)')));
c= length(find((EstCoef(PlcbTstWindow)>pvalue1H(PlcbTstWindow)')|(EstCoef(PlcbTstWindow)<pvalue1L(PlcbTstWindow)')));

Signi=-1*(EstCoef(EvntDate)<pvalue1L(EvntDate))-1*(EstCoef(EvntDate)<pvalue5L(EvntDate))-1*(EstCoef(EvntDate)<pvalue10L(EvntDate))+1*(EstCoef(EvntDate)>pvalue10H(EvntDate))+1*(EstCoef(EvntDate)>=pvalue5H(EvntDate))+1*(EstCoef(EvntDate)>=pvalue1H(EvntDate));

y=[EstCoef(EvntDate) pvalue1L(EvntDate) pvalue5L(EvntDate) pvalue10L(EvntDate)  pvalue10H(EvntDate)  pvalue5H(EvntDate)   pvalue1H(EvntDate)   a  b c length(Treatment) length(find(Treatment>0)) Signi]';
%
k=length(find(Treatment>0));

tmpnames=reshape(SynMatchNames,k*NumReport,1); 
tmpnames=tabulate(tmpnames);
tmpnames=tmpnames(:,1); 
tmpnames=sort(tmpnames);
WeightMatrix=cell(length(tmpnames)+1,k+1);
WeightMatrix(2:length(tmpnames)+1,1)=tmpnames;
WeightMatrix(1,2:k+1)=Names(Treatment>0)';
[m1 n1]=size(SynMatchNames);
for i=1:m1
    for j=1:n1
        Ind=find(strcmp(tmpnames,SynMatchNames(i,j))>0);
        WeightMatrix(Ind+1,i+1)={Weights(i,j)};
    end;
end;
