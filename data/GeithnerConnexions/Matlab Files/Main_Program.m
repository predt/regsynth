% This file is the main file that generates Table 5, 6, 7, 8(B) , 
% A10(B) , 
%A7 is based on SynMatchWeights output from SynMatch function

% The file uses EvntStdyLS and SynMatch functions for least square estimation of value of connections as well as the Synthetic matching estimator. 
tic
clear
addpath '\\ulysse\users\JL.HOUR\1A_These\A. Research\RegSynthProject\regsynth\data\GeithnerConnexions\Matlab Files'

load Data;       % Data file contains MarRet: daily market return , Re: daily return of financial firms, num: characteristics of financial firms, VarNames: name of the variables in num, ticker: name and ticker of the financial firms (it has the same order as variables "num" and "Re") , Statadate: stata date corresponding to dates of "Re" and "MarRet"
% note that in MarRet and Re the returns for nomination date are returns
% for 3 to 4 pm.
% This is the order of variables in num: {'geithner2num','ny','geithnerschednum','geithner07num','index_ds','name','ta2008','ta2008_log','roe2008','tdtc2008','mcap2008';}

%%%data cleaning
Ind=find((isnan(num(:,8))==1)|(isnan(num(:,9))==1)|(isnan(num(:,10))==1));   %eliminating firms with no data for 'ta2008_log','roe2008','tdtc2008'
num(Ind,:)=[];  Re(:,Ind)=[];  ticker(Ind,:)=[]; 

%treating missing information
Ind3=find(isnan(Re)==1);
Re(Ind3)=0;     

% Connection measure
ConnMeasure=3;  %1: Shared Board 2: NY Connection 3: Geithner Schedule 4: Geithner Schedule 2007 



Parameters.EstWindow=(GeiNomDat-281:GeiNomDat-32); % Window of 250 days ending 30 days prior to Geithner nomination
Parameters.NumMatch=20;   % This is the number of firms used in Synthetic matching. To avoid the dimensionality issue, we chose the N firms from the control group that are most correlated to the financial firm in treatment group and find the best synthetic match for that firm from a portfolio of those N firms.  
Parameters.alpha=inf;     % This parameters govern the correction for bad matches. inf means no correction. X means we get rid of firms for which the standard deviation of match errors is X times larger than the average standard deviation of matches for firms in Treatment group.
Parameters.NumIter=5000; % number of iterations for estimation of confidence intervals.
Parameters.EventDate=GeiNomDat-1;
Parameters.TstWindow=[340:353  357:389    395:447];   % This is the window that we use for falsification test. It should not include any event day.
Parameters.EvntWindowSize=0;   % The size of the window used for event study: 0 means only the event day, 1 mean returns on event day and the day after, ...
Parameters.NumReport=10;      % This is the number of weights in the synthetic matching to be reported

%Correlations calculations
a=find(num(:,5)==140);  % Citi Group 
b=find(num(:,5)==56);   % Bank of America

CorrCiti=(corr(Re(Parameters.EstWindow,a),Re(Parameters.EstWindow,:)))';
CorrBAC=(corr(Re(Parameters.EstWindow,b),Re(Parameters.EstWindow,:)))';

CorrCitiTr=sort(CorrCiti,'descend'); CorrCitiTr=CorrCitiTr(58);  
CorrBACTr=sort(CorrBAC,'descend'); CorrBACTr=CorrBACTr(58);

ticker=ticker(:,[2 1]);


%%%%% TABLE 5 %%%%%%%%%

%% Panel A

Parameters.EvntWindowSize=0; 
Parameters.alpha=inf;
%Col 1 - Col 3
Ind2=[];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl5PnlAClmn1=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlAClmn2=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlAClmn3=y([1 3 6 8:13]);

%Col 4 - Col 6
Parameters.alpha=inf; 
Ind2=find(num(:,3)>0 & num(:,3)<=2);

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl5PnlAClmn4=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlAClmn5=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlAClmn6=y([1 3 6 8:13]);

%Col 7 - Col 9
Parameters.alpha=inf; 
Ind2=find(num(:,3)>2);

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl5PnlAClmn7=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlAClmn8=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlAClmn9=y([1 3 6 8:13]);

Tbl5PnlA=[Tbl5PnlAClmn1 Tbl5PnlAClmn2 Tbl5PnlAClmn3 Tbl5PnlAClmn4 Tbl5PnlAClmn5 Tbl5PnlAClmn6 Tbl5PnlAClmn7 Tbl5PnlAClmn8 Tbl5PnlAClmn9];

%% Panel B
%Col 1 - Col 3

Parameters.EvntWindowSize=1; 
Parameters.alpha=inf;
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl5PnlBClmn1=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlBClmn2=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlBClmn3=y([1 3 6 8:13]);

%Col 4 - Col 6
Parameters.alpha=inf; 
Ind2=[find(CorrCiti>=CorrCitiTr); find(num(:,3)>0 & num(:,3)<=2)];

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl5PnlBClmn4=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlBClmn5=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlBClmn6=y([1 3 6 8:13]);

%Col 7 - Col 9
Parameters.alpha=inf; 
Ind2=[find(CorrCiti>=CorrCitiTr); find(num(:,3)>2)];

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl5PnlBClmn7=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlBClmn8=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlBClmn9=y([1 3 6 8:13]);

Tbl5PnlB=[Tbl5PnlBClmn1 Tbl5PnlBClmn2 Tbl5PnlBClmn3 Tbl5PnlBClmn4 Tbl5PnlBClmn5 Tbl5PnlBClmn6 Tbl5PnlBClmn7 Tbl5PnlBClmn8 Tbl5PnlBClmn9];

%% Panel C

%Col 1 - Col 3
Parameters.EvntWindowSize=10; 
Parameters.alpha=inf;
Ind2=find(CorrCiti>=CorrCitiTr);

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl5PnlCClmn1=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlCClmn2=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlCClmn3=y([1 3 6 8:13]);

%Col 4 - Col 6
Parameters.alpha=inf; 
Ind2=[find(CorrCiti>=CorrCitiTr); find(num(:,3)>0 & num(:,3)<=2)];

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl5PnlCClmn4=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlCClmn5=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlCClmn6=y([1 3 6 8:13]);

%Col 7 - Col 9
Parameters.alpha=inf; 
Ind2=[find(CorrCiti>=CorrCitiTr); find(num(:,3)>2)];

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl5PnlCClmn7=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlCClmn8=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl5PnlCClmn9=y([1 3 6 8:13]);

Tbl5PnlC=[Tbl5PnlCClmn1 Tbl5PnlCClmn2 Tbl5PnlCClmn3 Tbl5PnlCClmn4 Tbl5PnlCClmn5 Tbl5PnlCClmn6 Tbl5PnlCClmn7 Tbl5PnlCClmn8 Tbl5PnlCClmn9];

Tbl5=[Tbl5PnlA; zeros(1,9); Tbl5PnlB; zeros(1,9); Tbl5PnlC];


%%%% Table 6
%% Panel A

Parameters.EstWindow=297:320; % Financial Crisis Window
Parameters.alpha=inf; 
Parameters.EvntWindowSize=1; 

Ind2=find(CorrCiti>=CorrCitiTr);

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl6PnlAClmn1=yout([1:3 7 8]);

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlAClmn2=y([1 3 6 11:12]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlAClmn3=y([1 3 6 11:12]);

%Col 4 - Col 6
Parameters.alpha=inf; 
Ind2=[find(CorrCiti>=CorrCitiTr); find(num(:,3)>0 & num(:,3)<=2)];

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl6PnlAClmn4=yout([1:3 7 8]);

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlAClmn5=y([1 3 6 11:12]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlAClmn6=y([1 3 6 11:12]);

%Col 7 - Col 9
Parameters.alpha=inf; 
Ind2=[find(CorrCiti>=CorrCitiTr); find(num(:,3)>2)];

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl6PnlAClmn7=yout([1:3 7 8]);

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlAClmn8=y([1 3 6 11:12]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlAClmn9=y([1 3 6 11:12]);

Tbl6PnlA=[Tbl6PnlAClmn1 Tbl6PnlAClmn2 Tbl6PnlAClmn3 Tbl6PnlAClmn4 Tbl6PnlAClmn5 Tbl6PnlAClmn6 Tbl6PnlAClmn7 Tbl6PnlAClmn8 Tbl6PnlAClmn9];

%% Panel B
%Col 1 - Col 3
Parameters.EstWindow=(GeiNomDat-281:GeiNomDat-32); 
Parameters.alpha=inf; 
Parameters.EvntWindowSize=1; 

ConnMeasure=1;   % Shared Board
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl6PnlBClmn1=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlBClmn2=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlBClmn3=y([1 3 6 8:13]);

%Col 4 - Col 6
Parameters.alpha=inf; 
Ind2=[find(CorrCiti>=CorrCitiTr); find(num(:,ConnMeasure)>0 & num(:,ConnMeasure)<=1)];

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl6PnlBClmn4=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlBClmn5=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlBClmn6=y([1 3 6 8:13]);

%Col 7 - Col 9
Parameters.alpha=inf; 
Ind2=[find(CorrCiti>=CorrCitiTr); find(num(:,ConnMeasure)>1)];

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl6PnlBClmn7=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlBClmn8=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlBClmn9=y([1 3 6 8:13]);

Tbl6PnlB=[Tbl6PnlBClmn1 Tbl6PnlBClmn2 Tbl6PnlBClmn3 Tbl6PnlBClmn4 Tbl6PnlBClmn5 Tbl6PnlBClmn6 Tbl6PnlBClmn7 Tbl6PnlBClmn8 Tbl6PnlBClmn9];
Tbl6PnlB=Tbl6PnlB([1:3 7:8],:);

%% Panel C
%Col 1 - Col 3
Parameters.EstWindow=(GeiNomDat-281:GeiNomDat-32); 
Parameters.alpha=inf; 
Parameters.EvntWindowSize=1; 

ConnMeasure=2;  % NY measure
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl6PnlCClmn1=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlCClmn2=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlCClmn3=y([1 3 6 8:13]);

Tbl6PnlC=[Tbl6PnlCClmn1 Tbl6PnlCClmn2 Tbl6PnlCClmn3];
Tbl6PnlC=Tbl6PnlC([1:3 7:8],:);
%% Panel D

%Col 1 - Col 3
Parameters.EstWindow=(GeiNomDat-281:GeiNomDat-32); 
Parameters.alpha=inf; 
Parameters.EvntWindowSize=1; 

ConnMeasure=4;   % Shared Board
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl6PnlDClmn1=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlDClmn2=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlDClmn3=y([1 3 6 8:13]);

%Col 4 - Col 6
Parameters.alpha=inf; 
Ind2=[find(CorrCiti>=CorrCitiTr); find(num(:,ConnMeasure)>0 & num(:,ConnMeasure)<=1)];

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl6PnlDClmn4=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlDClmn5=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlDClmn6=y([1 3 6 8:13]);

%Col 7 - Col 9
Parameters.alpha=inf; 
Ind2=[find(CorrCiti>=CorrCitiTr); find(num(:,ConnMeasure)>1)];

Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;
Controls=[Tnum(:,[8,9,10]) Tnum(:,[8,9,10]).^2 Tnum(:,[8,9,10]).^3 ];

yout=EvntStdyLS(Parameters,Returns,MarRet,Controls,Treatment);
Tbl6PnlDClmn7=[yout; 0];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlDClmn8=y([1 3 6 8:13]);

Parameters.alpha=3; 
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
Tbl6PnlDClmn9=y([1 3 6 8:13]);

Tbl6PnlD=[Tbl6PnlDClmn1 Tbl6PnlDClmn2 Tbl6PnlDClmn3 Tbl6PnlDClmn4 Tbl6PnlDClmn5 Tbl6PnlDClmn6 Tbl6PnlDClmn7 Tbl6PnlDClmn8 Tbl6PnlDClmn9];
Tbl6PnlD=Tbl6PnlD([1:3 7:8],:);

%%%% Table A9  -- Pretrend

%% Panel A
load DailyReturns.mat;   % Replacing the retun on nomination date: in the previous file it is 3-4pm return for the nomination date
Re(:,Ind)=[];

Parameters.EstWindow=(GeiNomDat-281:GeiNomDat-32); % Window of 250 days ending 30 days prior to Geithner nomination
Parameters.alpha=inf; 
Parameters.EventDate=(GeiNomDat-1)-1;
Parameters.EvntWindowSize=1; 

ConnMeasure=3;
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlAClmn1=y([1 3 6 8:13]);

Ind2=[ find(num(:,ConnMeasure)==1|num(:,ConnMeasure)==2); find(CorrCiti>=CorrCitiTr)];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlAClmn2=y([1 3 6 8:13]);

Ind2=[find(num(:,ConnMeasure)>2); find(CorrCiti>=CorrCitiTr)];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlAClmn3=y([1 3 6 8:13]);

ConnMeasure=1;
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlAClmn4=y([1 3 6 8:13]);

ConnMeasure=2;
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlAClmn5=y([1 3 6 8:13]);

TblA9PnlA=[TblA9PnlAClmn1 TblA9PnlAClmn2 TblA9PnlAClmn3 TblA9PnlAClmn4 TblA9PnlAClmn5];
TblA9PnlA=TblA9PnlA([1:3 7:8],:);

%% Panel B
Parameters.EventDate=(GeiNomDat-1)-5;
Parameters.EvntWindowSize=5; 

ConnMeasure=3;
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlBClmn1=y([1 3 6 8:13]);

Ind2=[ find(num(:,ConnMeasure)==1|num(:,ConnMeasure)==2); find(CorrCiti>=CorrCitiTr)];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlBClmn2=y([1 3 6 8:13]);

Ind2=[find(num(:,ConnMeasure)>2); find(CorrCiti>=CorrCitiTr)];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlBClmn3=y([1 3 6 8:13]);

ConnMeasure=1;
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlBClmn4=y([1 3 6 8:13]);

ConnMeasure=2;
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlBClmn5=y([1 3 6 8:13]);

TblA9PnlB=[TblA9PnlBClmn1 TblA9PnlBClmn2 TblA9PnlBClmn3 TblA9PnlBClmn4 TblA9PnlBClmn5];
TblA9PnlB=TblA9PnlB([1:3 7:8],:);

%% Panel C
Parameters.EventDate=(GeiNomDat-1)-10;
Parameters.EvntWindowSize=10; 

ConnMeasure=3;
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlCClmn1=y([1 3 6 8:13]);

Ind2=[ find(num(:,ConnMeasure)==1|num(:,ConnMeasure)==2); find(CorrCiti>=CorrCitiTr)];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlCClmn2=y([1 3 6 8:13]);

Ind2=[find(num(:,ConnMeasure)>2); find(CorrCiti>=CorrCitiTr)];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlCClmn3=y([1 3 6 8:13]);

ConnMeasure=1;
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlCClmn4=y([1 3 6 8:13]);

ConnMeasure=2;
Ind2=find(CorrCiti>=CorrCitiTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA9PnlCClmn5=y([1 3 6 8:13]);

TblA9PnlC=[TblA9PnlCClmn1 TblA9PnlCClmn2 TblA9PnlCClmn3 TblA9PnlCClmn4 TblA9PnlCClmn5];
TblA9PnlC=TblA9PnlC([1:3 7:8],:);


%%
%%%% Table A13   Tax Problem
%%
Parameters.EstWindow=(GeiNomDat-281:GeiNomDat-32); % Window of 250 days ending 30 days prior to Geithner nomination
Parameters.alpha=inf; 
Parameters.EventDate=389;
Parameters.EvntWindowSize=1; 

ConnMeasure=3;
Ind2=find(CorrCiti>=CorrCitiTr | CorrBAC>=CorrBACTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn1=y([1 3 6 8:13]);

Ind2=[ find(num(:,ConnMeasure)==1|num(:,ConnMeasure)==2);find(CorrCiti>=CorrCitiTr | CorrBAC>=CorrBACTr)];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn2=y([1 3 6 8:13]);

Ind2=[find(num(:,ConnMeasure)>2); find(CorrCiti>=CorrCitiTr | CorrBAC>=CorrBACTr)];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn3=y([1 3 6 8:13]);


Ind2=[];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn4=y([1 3 6 8:13]);

Ind2=find(num(:,ConnMeasure)==1|num(:,ConnMeasure)==2);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn5=y([1 3 6 8:13]);

Ind2=find(num(:,ConnMeasure)>2); 
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn6=y([1 3 6 8:13]);


Parameters.EvntWindowSize=3; 

Ind2=find(CorrCiti>=CorrCitiTr | CorrBAC>=CorrBACTr);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn7=y([1 3 6 8:13]);

Ind2=[ find(num(:,ConnMeasure)==1|num(:,ConnMeasure)==2);find(CorrCiti>=CorrCitiTr | CorrBAC>=CorrBACTr)];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn8=y([1 3 6 8:13]);

Ind2=[find(num(:,ConnMeasure)>2); find(CorrCiti>=CorrCitiTr | CorrBAC>=CorrBACTr)];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn9=y([1 3 6 8:13]);


Ind2=[];
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn10=y([1 3 6 8:13]);

Ind2=find(num(:,ConnMeasure)==1|num(:,ConnMeasure)==2);
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn11=y([1 3 6 8:13]);

Ind2=find(num(:,ConnMeasure)>2); 
Tnum=num; Tticker=ticker(:,2); TRe=Re;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; 
Returns=TRe; Treatment=Tnum(:,ConnMeasure)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA131PnlBClmn12=y([1 3 6 8:13]);

TblA131PnlB=[TblA131PnlBClmn1 TblA131PnlBClmn2 TblA131PnlBClmn3 TblA131PnlBClmn4 TblA131PnlBClmn5 TblA131PnlBClmn6 TblA131PnlBClmn7 TblA131PnlBClmn8 TblA131PnlBClmn9 TblA131PnlBClmn10 TblA131PnlBClmn11 TblA131PnlBClmn12];
TblA131PnlB=TblA131PnlB([1:3 7:8],:);

%%
%% Other Connections
%Table A11 PanelB
% Holder  11/19/2008
load OtherConnection;
Parameters.EventDate=353-1;
Parameters.EstWindow=(Parameters.EventDate-280:Parameters.EventDate-31); % Window of 250 days ending 30 days prior to Geithner nomination

Parameters.EvntWindowSize=1; 
Parameters.alpha=inf;
Ind2=[find(CorrCiti>=.99) 359];  %359 is MXB
%Ind2=[];
Tnum=num; Tticker=ticker(:,2); TRe=Re;   TOtherConnections=OtherConnections;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; TOtherConnections(Ind2,:)=[];

Returns=TRe; Treatment=TOtherConnections(:,8)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA11PnlBClmn1=y([1 3 6 8:13]);

%Vilsack 12/17/2008
Parameters.EventDate=373-1;
Parameters.EstWindow=(Parameters.EventDate-280:Parameters.EventDate-31); % Window of 250 days ending 30 days prior to Geithner nomination

Parameters.EvntWindowSize=1; 
Parameters.alpha=inf;
Ind2=find(CorrCiti>=.99);
%Ind2=[];
Tnum=[num; num_fabk_ebsb(1,:)]; Tticker=ticker(:,2); TRe=[Re Re_FABK_EBSB(:,1)];  TOtherConnections=OtherConnections;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; TOtherConnections(Ind2,:)=[];

Returns=TRe; Treatment=[TOtherConnections(:,9)>0;0]; Names=[Tticker; 'FABK'];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA11PnlBClmn2=y([1 3 6 8:13]);

%Locke 2/24/2009
Parameters.EventDate=417-1;
Parameters.EstWindow=(Parameters.EventDate-280:Parameters.EventDate-31); % Window of 250 days ending 30 days prior to Geithner nomination

Parameters.EvntWindowSize=1; 
Parameters.alpha=inf;
Ind2=find(CorrCiti>=.99);
%Ind2=[];
Tnum=[num; num_fabk_ebsb]; Tticker=ticker(:,2); TRe=[Re Re_FABK_EBSB];  TOtherConnections=OtherConnections;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; TOtherConnections(Ind2,:)=[];

Returns=TRe; Treatment=[TOtherConnections(:,10)>0;0;0]; Names=[Tticker; 'FABK' ; 'EBSB'];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA11PnlBClmn3=y([1 3 6 8:13]);

%Chu 12/11/2008
Parameters.EventDate=368-1;
Parameters.EstWindow=(Parameters.EventDate-280:Parameters.EventDate-31); % Window of 250 days ending 30 days prior to Geithner nomination

Parameters.EvntWindowSize=1; 
Parameters.alpha=inf;
Ind2=find(CorrCiti>=.99);
%Ind2=[];
Tnum=[num; num_fabk_ebsb(1,:)]; Tticker=ticker(:,2); TRe=[Re Re_FABK_EBSB(:,1)];  TOtherConnections=OtherConnections;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; TOtherConnections(Ind2,:)=[];

Returns=TRe; Treatment=[TOtherConnections(:,11)>0;0]; Names=[Tticker; 'FABK'];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA11PnlBClmn4=y([1 3 6 8:13]);

%Duncan 12/16/2008
Parameters.EventDate=371-1;
Parameters.EstWindow=(Parameters.EventDate-280:Parameters.EventDate-31); % Window of 250 days ending 30 days prior to Geithner nomination

Parameters.EvntWindowSize=1; 
Parameters.alpha=inf;
Ind2=find(CorrCiti>=.99);
%Ind2=[];
Tnum=[num; num_fabk_ebsb(1,:)]; Tticker=ticker(:,2); TRe=[Re Re_FABK_EBSB(:,1)];  TOtherConnections=OtherConnections;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; TOtherConnections(Ind2,:)=[];

Returns=TRe; Treatment=[TOtherConnections(:,12)>0;0]; Names=[Tticker; 'FABK'];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA11PnlBClmn5=y([1 3 6 8:13]);

%Shinseki 12/8/2008
Parameters.EventDate=365-1;
Parameters.EstWindow=(Parameters.EventDate-280:Parameters.EventDate-31); % Window of 250 days ending 30 days prior to Geithner nomination

Parameters.EvntWindowSize=1; 
Parameters.alpha=inf;
Ind2=find(CorrCiti>=.99);
%Ind2=[];
Tnum=[num; num_fabk_ebsb(1,:)]; Tticker=ticker(:,2); TRe=[Re Re_FABK_EBSB(:,1)];  TOtherConnections=OtherConnections;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; TOtherConnections(Ind2,:)=[];

Returns=TRe; Treatment=[TOtherConnections(:,13)>0;0]; Names=[Tticker; 'FABK'];

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA11PnlBClmn6=y([1 3 6 8:13]);

%Napolitano 11/20/2008
Parameters.EventDate=354-1;
Parameters.EstWindow=(Parameters.EventDate-280:Parameters.EventDate-31); % Window of 250 days ending 30 days prior to Geithner nomination

Parameters.EvntWindowSize=1; 
Parameters.alpha=inf;
Ind2=find(CorrCiti>=.99);
%Ind2=[];
Tnum=num; Tticker=ticker(:,2); TRe=Re;  TOtherConnections=OtherConnections;
Tnum(Ind2,:)=[];    TRe(:,Ind2)=[];  Tticker(Ind2)=[]; TOtherConnections(Ind2,:)=[];
Returns=TRe; Treatment=TOtherConnections(:,14)>0; Names=Tticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Treatment,Names);
TblA11PnlBClmn7=y([1 3 6 8:13]);

TblA11PnlB=[TblA11PnlBClmn1 TblA11PnlBClmn2 TblA11PnlBClmn3 TblA11PnlBClmn4 TblA11PnlBClmn5 TblA11PnlBClmn6 TblA11PnlBClmn7];

%% Table 7 Panel B--------- CDS Results
%%
load CDSData


Parameters.EstWindow=101:201; % Window of 250 days ending 30 days prior to Geithner nomination
Parameters.alpha=inf; 
Parameters.NumIter=5000; % number of iterations for placebo test
Parameters.EventDate=233;
Parameters.TstWindow=[150:232 235:250];
Parameters.EvntWindowSize=0; 
Parameters.NumReport=10;
Parameters.NumMatch=20;
[m,n]=size(Re); 
% Cleaning data for missing values. % We take out each firm for which more than 50% of CDS data is
% missing
th=.5;      Ind2=[]; 
for j=1:n
    Ind3=find(isnan(Re(:,j))==1);
    if length(Ind3)>th*m
        Ind2=[Ind2 j];
    elseif isnan(Re(233,j))
        Ind2=[Ind2 j];
    else
        Re(Ind3)=0;
    end;
end;

Re(:,Ind2)=[];  Re(isnan(Re)==1)=0;     ticker(Ind2)=[];
Conn(Ind2,:)=[];
%Ind2=find(std(Re)==0);  Re(:,Ind2)=[];  ticker(Ind2)=[]; Conn(Ind2,:)=[];

%Citigroup Included
ConnMeasure=3;
Returns=Re;  Names=ticker;

[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn1=y([1 3 6 8:13]);

ConnMeasure=1;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn2=y([1 3 6 8:13]);


ConnMeasure=2;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn3=y([1 3 6 8:13]);

Parameters.EvntWindowSize=9;

ConnMeasure=3;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn7=y([1 3 6 8:13]);

ConnMeasure=1;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn8=y([1 3 6 8:13]);


ConnMeasure=2;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn9=y([1 3 6 8:13]);

%Citigroup Excluded
Parameters.EvntWindowSize=0;
I=find(strcmp(ticker,'C')==1);
Returns(:,I)=[];
Conn(I,:)=[];
Names(I)=[];

ConnMeasure=3;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn4=y([1 3 6 8:13]);

ConnMeasure=1;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn5=y([1 3 6 8:13]);


ConnMeasure=2;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn6=y([1 3 6 8:13]);

Parameters.EvntWindowSize=9;

ConnMeasure=3;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn10=y([1 3 6 8:13]);

ConnMeasure=1;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn11=y([1 3 6 8:13]);


ConnMeasure=2;
[y,SynMatchWeights]=SynMatch(Parameters,Returns,Conn(:,ConnMeasure)>0,Names);
Tbl7PnlBClmn12=y([1 3 6 8:13]);


Tbl7PnlB=[Tbl7PnlBClmn1 Tbl7PnlBClmn2 Tbl7PnlBClmn3 Tbl7PnlBClmn4 Tbl7PnlBClmn5 Tbl7PnlBClmn6 Tbl7PnlBClmn7 Tbl7PnlBClmn8 Tbl7PnlBClmn9 Tbl7PnlBClmn10 Tbl7PnlBClmn11 Tbl7PnlBClmn12];
Tbl7PnlB=Tbl7PnlB([1:3 7:8],:);


toc