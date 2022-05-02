% =========================================================================
% ----------------- RIVER DISCHARGE IN SOUTH AMERICA: ---------------------
% ------ AGREEMENT BETWEEN PROJECTED CHANGE AND RECENT ALTERATION ---------
% =========================================================================
% AUTHORS:
% 1 - PEDRO TORRES MIRANDA: pedrotorresm121@gmail.com
% 2 - RODRIGO CAUDURO DIAS DE PAIVA: rodrigo.paiva@ufrgs.br
% 3 - CLÉBER HENRIQUE DE ARAÚJO GAMA: cleber.hag@gmail.com
% 4 - JOÃO PAULO LYRA FIALHO BRÊDA: joaopaulolfb@gmail.com
% =========================================================================

nMini=33749;
n1=20; %1980-1999
n2=20; %2000-2020
%---------------------- RIVER DISCHARGE 1980-2020 --------------------------
InitialDate=datenum('1-1-1979','mm-dd-yyyy');
EndDate=datenum('12-31-2019','mm-dd-yyyy');
NT=EndDate-InitialDate+1;

fid=fopen('QTUDO.MGB','r');
    flows=fread(fid,[nMini,NT],'single');
    flows=flows';
fclose(fid);

InitialDate2 = datenum('1-1-1980','mm-dd-yyyy');

dates=(InitialDate2:EndDate)';

descarteFlows = InitialDate2-InitialDate;

flows(1:descarteFlows,:)=[];
%--------------------------------------------------------------------------

%------------------- ANUAL MEAN DISCHARGES 1980-2020 ----------------------
datesVEC=datevec(dates);
anos=datesVEC(:,1);
ano=(min(anos):max(anos))';

for i=1:length(ano)
    i
    ni=find(anos==ano(i),1, 'first');
    nf=find(anos==ano(i),1, 'last');
    qmin(i,:)=min(flows(ni:nf,:));
    qmean(i,:)=mean(flows(ni:nf,:));
    qmax(i,:)=max(flows(ni:nf,:));
end

for i=1:length(ano)
    ndias=find(anos==ano(i));
    dias_ano(i,1)=length(ndias);
end

qmeanXdias=qmean.*dias_ano;
qmean1=sum(qmeanXdias(1:n1,:))./sum(dias_ano(1:n1,:));
qmean2=sum(qmeanXdias(n1+1:end,:)./sum(dias_ano(n1+1:end,:)));

qmin1=mean(qmin(1:n1,:));
qmin2=mean(qmin(n1+1:end,:));
% qmean1=mean(qmean(1:n1,:));
% qmean2=mean(qmean(n1+1:end,:));
qmax1=mean(qmax(1:n1,:));
qmax2=mean(qmax(n1+1:end,:));
%---------------------------------------------------------------------------

%---------- DISCHARGE ALTERATION BETWEEN 1980-1999 E 2000-2020 -------------

changeAbsMin = (qmin2-qmin1)';
changeAbsMean = (qmean2-qmean1)';
changeAbsMax = (qmax2-qmax1)';

changePerMin = 100*changeAbsMin./abs(qmin1');
changePerMean = 100*changeAbsMean./abs(qmean1');
changePerMax = 100*changeAbsMax./abs(qmax1');

signChangeMin = sign(changePerMin);
signChangeMean = sign(changePerMean);
signChangeMax = sign(changePerMax);
%-----------------------------------------------------------------------------

%--------------- SIGNIFICANT CHANGE: T-TEST AND MANN KENDALL -----------------
% SERIE'S SIZE (n) AND LEVEL OF SIGNIFICANCE (a)
n=min(size(qmax));
a = .05;
z = norminv(1-a/2);
varS=(n*(n-1)*(2*n+5))/18;% var(s) mann-kendall not considering repeated values

qmin1=qmin(1:n1,:);
qmin2=qmin(n1+1:end,:);
qmean1=qmean(1:n1,:);
qmean2=qmean(n1+1:end,:);
qmax1=qmax(1:n1,:);
qmax2=qmax(n1+1:end,:);

%STUDENT'S T-TEST
for i=1:nMini
    ttestqmin(i,1)=ttest2(qmin1(:,i),qmin2(:,i)); %default p-value = 5%
    ttestqmean(i,1)=ttest2(qmean1(:,i),qmean2(:,i)); %default p-value = 5%
    ttestqmax(i,1)=ttest2(qmax1(:,i),qmax2(:,i)); %default p-value = 5%
end

%SIGN ATTRIBUTION TO T-TEST
ttestqmin=ttestqmin.*signChangeMin;
ttestqmean=ttestqmean.*signChangeMean;
ttestqmax=ttestqmax.*signChangeMax;

%loads series with autocorrelation correction for MK's application
load('autocorrelacao', 'qminMK', 'qmeanMK', 'qmaxMK')

for k=1:nMini
    k
    %MATRIX (Xj-Xi)
    for i=1:n
        for j=1:n
            if j<=i
                Smin(j,i)=0;
                Smean(j,i)=0;
                Smax(j,i)=0;
            else
                Smin(j,i)=sign(qminMK(j,k)-qminMK(i,k));
                Smean(j,i)=sign(qmeanMK(j,k)-qmeanMK(i,k));
                Smax(j,i)=sign(qmaxMK(j,k)-qmaxMK(i,k));
            end
        end
    end

    %S PARAMETER DEFINITION
    sumMin(1,k)=sum(sum(Smin));
    sumMean(1,k)=sum(sum(Smean));
    sumMax(1,k)=sum(sum(Smax));

    %Zs PARAMETER DEFINITION
    %qmin
    if sumMin(1,k) > 0
        Zmin(1,k) = (sumMin(1,k)-1)/(varS^(1/2));
    elseif sumMin(1,k) < 0
        Zmin(1,k) = (sumMin(1,k)+1)/(varS^(1/2));
    else
        Zmin(1,k) = 0;
    end
    
    %qmean
    if sumMean(1,k) > 0
        Zmean(1,k) = (sumMean(1,k)-1)/(varS^(1/2));
    elseif sumMean(1,k) < 0
        Zmean(1,k) = (sumMean(1,k)+1)/(varS^(1/2));
    else
        Zmean(1,k) = 0;
    end

    %qmax
    if sumMax(1,k) > 0
        Zmax(1,k) = (sumMax(1,k)-1)/(varS^(1/2));
    elseif sumMax(1,k) < 0
        Zmax(1,k) = (sumMax(1,k)+1)/(varS^(1/2));
    else
        Zmax(1,k) = 0;
    end
    
    %TREND
    %qmin
    if Zmin(1,k) > 0 && abs(Zmin(1,k)) > z
        MKtendmin(k,1) = 1;
    elseif Zmin(1,k) < 0 && abs(Zmin(1,k)) > z
        MKtendmin(k,1) = -1;
    else
        MKtendmin(k,1) = 0;
    end
    
    %qmean
    if Zmean(1,k) > 0 && abs(Zmean(1,k)) > z
        MKtendmean(k,1) = 1;
    elseif Zmean(1,k) < 0 && abs(Zmean(1,k)) > z
        MKtendmean(k,1) = -1;
    else
        MKtendmean(k,1) = 0;
    end
    
    %qmax
    if Zmax(1,k) > 0 && abs(Zmax(1,k)) > z
        MKtendmax(k,1) = 1;
    elseif Zmax(1,k) < 0 && abs(Zmax(1,k)) > z
        MKtendmax(k,1) = -1;
    else
        MKtendmax(k,1) = 0;
    end
    
    %SEN'S SLOPE
    l=1;
    for i=1:n
        for j=1:n
            if i==j || i<j
            else
                Min(l,k)=(qmin(i,k)-qmin(j,k))/(i-j);
                Mean(l,k)=(qmean(i,k)-qmean(j,k))/(i-j);
                Max(l,k)=(qmax(i,k)-qmax(j,k))/(i-j);
                l=l+1;
            end
        end
    end
end

%SEN'S SLOPE
sen_qmin=(median(Min))';
sen_qmean=(median(Mean))';
sen_qmax=(median(Max))';

%ALTERATION BETWEEN 1980-2020 BY SEN'S SLOPE
alt_sen_qmin=sen_qmin.*n;
alt_sen_qmean=sen_qmean.*n;
alt_sen_qmax=sen_qmax.*n;

%-------------------- LINEAR REGRESSION ---------------------------------
for k=1:nMini
    k
    %x=i
    for i=1:n
        x(i,1)=i;
        xqmin(i,k)=i*qmin(i,k);
        xqmean(i,k)=i*qmean(i,k);
        xqmax(i,k)=i*qmax(i,k);
        x2(i,1)=i^2;
    end
end

%sum
sx=sum(x);
sxqmin=sum(xqmin);
sxqmean=sum(xqmean);
sxqmax=sum(xqmax);
sx2=sum(x2);
sqmin=sum(qmin);
sqmean=sum(qmean);
sqmax=sum(qmax);

%y = m.x + b
%m
m_min=(sxqmin-sx*sqmin/n)/(sx2-sx^2/n);
m_mean=(sxqmean-sx*sqmean/n)/(sx2-sx^2/n);
m_max=(sxqmax-sx*sqmax/n)/(sx2-sx^2/n);

%b
b_min=sqmin./n-m_min*sx./n;
b_mean=sqmean./n-m_mean*sx./n;
b_max=sqmax./n-m_max*sx./n;

for i=1:n
    Lqmin(i,:)=m_min*i+b_min;
    Lqmean(i,:)=m_mean*i+b_mean;
    Lqmax(i,:)=m_max*i+b_max;
end

alt_lin_qmin=(100*(Lqmin(end,:)-Lqmin(1,:))./Lqmin(1,:))';
alt_lin_qmean=(100*(Lqmean(end,:)-Lqmean(1,:))./Lqmean(1,:))';
alt_lin_qmax=(100*(Lqmax(end,:)-Lqmax(1,:))./Lqmax(1,:))';

%------------------- T-TEST + MANN-KENDALL (EXCLUSIONARY) --------------------
for i=1:nMini
    if MKtendmin(i)==ttestqmin(i)
        ttestMKqmin(i,1)=ttestqmin(i);
    else
        ttestMKqmin(i,1)=0;
    end
    
    if MKtendmean(i)==ttestqmean(i)
        ttestMKqmean(i,1)=ttestqmean(i);
    else
        ttestMKqmean(i,1)=0;
    end
    
    if MKtendmax(i)==ttestqmax(i)
        ttestMKqmax(i,1)=ttestqmax(i);
    else
        ttestMKqmax(i,1)=0;
    end
end

%------------------- T-TEST + MANN-KENDALL (INCLUSIVE) --------------------
for i=1:nMini
    if MKtendmin(i)==ttestqmin(i) && ttestqmin(i)==0 
        combTtestMKqmin(i,1)=0;
    elseif MKtendmin(i)~= 0
        combTtestMKqmin(i,1)=MKtendmin(i);
    elseif ttestqmin(i)~= 0
        combTtestMKqmin(i,1)=ttestqmin(i);
    end

    if MKtendmean(i)==ttestqmean(i) && ttestqmean(i)==0 
        combTtestMKqmean(i,1)=0;
    elseif MKtendmean(i)~= 0
        combTtestMKqmean(i,1)=MKtendmean(i);
    elseif ttestqmean(i)~= 0
        combTtestMKqmean(i,1)=ttestqmean(i);
    end
    
    if MKtendmax(i)==ttestqmax(i) && ttestqmax(i)==0 
        combTtestMKqmax(i,1)=0;
    elseif MKtendmax(i)~= 0
        combTtestMKqmax(i,1)=MKtendmax(i);
    elseif ttestqmax(i)~= 0
        combTtestMKqmax(i,1)=ttestqmax(i);
    end
end

%NEUTRAL VALUES
%significant neutral values
nSignMin=find(combTtestMKqmin==0);
nSignMean=find(combTtestMKqmean==0);
nSignMax=find(combTtestMKqmax==0);

%mean neutral values (|q_alt| < +/- 10%)
nAvgMin=find(abs(changePerMin)<=10);
nAvgMean=find(abs(changePerMean)<=10);
nAvgMax=find(abs(changePerMax)<=10);

%mean AND significant neutral values
nSAvgMin=find(abs(changePerMin(nSignMin))<=10);
nSAvgMean=find(abs(changePerMean(nSignMean))<=10);
nSAvgMax=find(abs(changePerMax(nSignMax))<=10);

%(mean neutral)/(significant neutral)
pNAvgMin = length(nAvgMin)/length(nSignMin)*100;
pNAvgMean = length(nAvgMean)/length(nSignMean)*100;
pNAvgMax = length(nAvgMax)/length(nSignMax)*100;

%(mean AND significant neutral)/(significant neutral)
pNSAvgMin = length(nSAvgMin)/length(nSignMin)*100;
pNSAvgMean = length(nSAvgMean)/length(nSignMean)*100;
pNSAvgMax = length(nSAvgMax)/length(nSignMax)*100;

%(mean AND significant neutral)/(mean neutral)
pNSAvg_AvgMin = length(nSAvgMin)/length(nAvgMin)*100;
pNSAvg_AvgMean = length(nSAvgMean)/length(nAvgMean)*100;
pNSAvg_AvgMax = length(nSAvgMax)/length(nAvgMax)*100;



neutralMin=median(abs(changePerMin(nSignMin)));
neutralMean=median(abs(changePerMean(nSignMean)));
neutralMax=median(abs(changePerMax(nSignMax)));

neutral2Min=mean(abs(changePerMin(nSignMin)));
neutral2Mean=mean(abs(changePerMean(nSignMean)));
neutral2Max=mean(abs(changePerMax(nSignMax)));
%-----------------------

%NUMBER OF MINIBASINS WITH SIGNIFICANT CHANGE (MK + T-TEST)
ministtestmin=100*sum(abs(ttestqmin))/nMini;
ministtestmean=100*sum(abs(ttestqmean))/nMini;
ministtestmax=100*sum(abs(ttestqmax))/nMini;
minisMKmin=100*sum(abs(MKtendmin))/nMini;
minisMKmean=100*sum(abs(MKtendmean))/nMini;
minisMKmax=100*sum(abs(MKtendmax))/nMini;

%ALTERATION OF STANDARD DEVIATION BETWEEN 1980-1999 AND 2000-2019
std_qmin1=(std(qmin(1:n1,:)))';
std_qmean1=(std(qmean(1:n1,:)))';
std_qmax1=(std(qmax(1:n1,:)))';

std_qmin2=(std(qmin(n1+1:end,:)))';
std_qmean2=(std(qmean(n1+1:end,:)))';
std_qmax2=(std(qmax(n1+1:end,:)))';

dif_std_qmin=100*(std_qmin2-std_qmin1)./std_qmin1;
dif_std_qmean=100*(std_qmean2-std_qmean1)./std_qmean1;
dif_std_qmax=100*(std_qmax2-std_qmax1)./std_qmax1;

save vazoes1980_2020