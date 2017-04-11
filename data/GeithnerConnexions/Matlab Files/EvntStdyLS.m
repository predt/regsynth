% This function is used for event study using least square estimation
% method and correcting for the VarCov matrix of standard errors.It also
% perform placebo test for non-event days and reports the number of false
% positives at different significance levels. 
function yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment)

EstWindow=Parameters.EstWindow;     % Estimation window used for estimating market model as well as the VarCov matrix 
EvntDate=Parameters.EventDate;      %Date of actual event
PlcbTstWindow=Parameters.TstWindow;  %Placebo Test Window
EvntWindowSize=Parameters.EvntWindowSize;  %This is the size of the event study window. This can be 0, 1, 2, ...0 means only the event date, 1 means the event date and the date after...

%Returns is individual firms stock return
%MarRet is the market return
%Controls are controls in the least square regression
%Treatment is a measure of exposure to treatment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now Computing Abnormal Returns for entire sample
[m n]=size(Returns);
for i=1:n
    x=MarRet(EstWindow);
    y=Returns(EstWindow,i);
    [beta,tst]=rols(y,x,1);
    AbnRet(:,i)=Returns(:,i)-beta(1)-beta(2)*MarRet;
end;
[m1 n1]=size(AbnRet);
CAR=AbnRet(1:m1-EvntWindowSize-1,:);
for j=2:EvntWindowSize+1
    CAR=CAR+AbnRet(j:m1-EvntWindowSize-2+j,:);
end;

%Estimating VarCov matrix
k=0;
for Date=EstWindow(1):EstWindow(length(EstWindow))
    k=k+1;
    y=CAR(Date,:)';
    [beta,tstat]=rols(y,Controls,1);
    Err(:,k)=y-[ones(length(y),1) Controls]*beta;
end;

Err(isnan(Err)==1)=0;
VarCov=cov(Err');

[m1 n1]=size(CAR);
for Date=1:m1
    y=CAR(Date,:)';
    [beta,tstat]=fgls(y,[Treatment Controls],1,VarCov);
    Bet(Date)=beta(1);
    Tstat(Date)=tstat(1);
end;

SigniLevel=-1*(Tstat<-tinv(.995,length(y)))-1*(Tstat<-tinv(.975,length(y)))-1*(Tstat<-tinv(.95,length(y)))+1*(Tstat>tinv(.95,length(y)))+1*(Tstat>tinv(.975,length(y)))+1*(Tstat>tinv(.995,length(y)));
a= length(find(abs(SigniLevel(PlcbTstWindow))>2));
b= length(find(abs(SigniLevel(PlcbTstWindow))>1));
c= length(find(abs(SigniLevel(PlcbTstWindow))>0));

yout=[Bet(EvntDate); SigniLevel(EvntDate); Tstat(EvntDate) ; c;b;a; length(y); length(find(Treatment>0)) ];
