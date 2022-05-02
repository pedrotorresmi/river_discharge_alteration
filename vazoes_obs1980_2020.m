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

%------------------------ MEAN ALTERATION MGB-SA -----------------------------------
load('vazoes1980_2020','changePerMin','changePerMean','changePerMax')

%------------------------- OBSERVATION DATA ---------------------------------
load('Nova_selecao_postos_larissa.mat', 'selecao_ANA');
ANA = struct2table(selecao_ANA);
sortedANA = sortrows(ANA, 'mini');
ANA = table2struct(sortedANA);

%------------------- INITIAL DEFINITIONS ---------------
nT1=datenum('1-Jan-1980');
nT2=datenum('31-Dec-1999');
nT3=datenum('1-Jan-2000');
nT4=datenum('31-Dec-2019');

n1=20; %1980-2000
n2=20; %2000-2020

q_postos=ones(nT4-nT1+1,length(ANA))*NaN;
%----------------------------------------------------------

%------------------ MATRIX OF ANA STATIONS' DISCHARGES BETWEEN 1980-2020 -----------------
for i=1:length(ANA)
    i
    data=[];
    k=1;
    for j=1:length(ANA(i).Datas)        
        data(k)=datenum(ANA(i).Datas(j));
        k=k+1;
    end
    data=data-nT1+1;
    q=ANA(i).Vazao;
    ii=find(data<1|data>nT4-nT1+1);
    data(ii)=[];
    q(ii)=[];    
    q_postos(data,i)=q;    
end

% load('vazoes_obs1980_2020','q_postos'); %carrega a matriz qpostos
%------------------------------------------------------------

%----------------------------- ANA STATIONS' SELECTION -----------------------------------
data=[];
data=(nT1:nT4)';
anos=datevec(data);
anos=anos(:,1);
ano=unique(anos);

%PERCENTUAL OF ANUAL DATA IN EACH STATION
for i=1:length(ano)
    x=[];
    x=find(anos==ano(i));
    dias=length(x);
    nan(1,:)=sum(isnan(q_postos(x,:)));
    dados_ano(i,:)=(1-nan./dias)*100;
end
%------

%DISPOSAL OF STATIONS WITH LESS THAN 75% YEARS WITH DATA (YEAR WITH DATA = year with less than 20% of data gaps)
k=1;
for i=1:min(size(q_postos))
    nan=[];
    nan=find(dados_ano(:,i)<80);
    dados_ano(nan,i)=NaN;
    falha1=100*sum(isnan(dados_ano(1:n1,i)))/20;
    falha2=100*sum(isnan(dados_ano(n1+1:end,i)))/20;
    if falha1 > 25 || falha2 > 25
        descarte(k)=i;
        k=k+1;
    else
    end
end
%-----

q_postos_ok=q_postos;
q_postos_ok(:,descarte)=[]; %matrix of filtered stations
q1=q_postos_ok(1:nT2-nT1+1,:); %matrix with daily river discharge 1980-2000
q2=q_postos_ok(nT3-nT1+1:nT4-nT1+1,:); %matrix with daily river discharge 2000-2020
%------------------------------------------------------------------------------------------

%----------------- MINIMUM, MEAN AND MAXIMUM DISCHARGE ALTERATION -------------------------
%SELECTED STATIONS' MINIBASINS
for i=1:length(ANA)
    minis(i)=ANA(i).mini;
end
minis_ok=minis';
minis_ok(descarte)=[];
%---------

%DATA PERCENTAGE FROM SELECTED STATIONS
dados_ano(:,descarte)=[];
%---------

%MEAN VALUE OF MINIMUM, MEAN AND MAXIMUM DISCHARGES FOR 1980-2000 AND FOR 2000-2020
%MEAN DISCHARGE
for i=1:length(minis_ok)
    j=[];
    j=find(isnan(q1(:,i)));
    media=q1(:,i);
    media(j)=[];
    qmean1(i,1)=mean(media);
    j=[];
    j=find(isnan(q2(:,i)));
    media=q2(:,i);
    media(j)=[];
    qmean2(i,1)=mean(media);
end
%--------

%MINIMUM AND MAXIMUM DISCHARGES
for i=1:length(ano)
    i
    ni=[];
    nf=[];
    ni=find(anos==ano(i),1, 'first');
    nf=find(anos==ano(i),1, 'last');
%     dados_ano(i,:)=100-100*sum(isnan(q_postos_ok(ni:nf,:)))/(nf-ni+1);
    qmax(i,:)=max(q_postos_ok(ni:nf,:));
    qmin(i,:)=min(q_postos_ok(ni:nf,:));
end

anos100max = ones(length(ano),length(minis_ok))*NaN;
anos100min = ones(length(ano),length(minis_ok))*NaN;

for i=1:length(minis_ok)
    j=[];
    j=find(dados_ano(:,i)==100);
    anos100max(j,i)=qmax(j,i);
    anos100min(j,i)=qmin(j,i);
end

for i=1:length(minis_ok)
    j=[];
    j=find(isnan(anos100max(1:n1,i)));
    max=anos100max(1:n1,i);
    max(j)=[];
    qmax1(i,1)=mean(max);
    
    j=[];
    j=find(isnan(anos100max(n1+1:end,i)));
    max=anos100max(n1+1:end,i);
    max(j)=[];
    qmax2(i,1)=mean(max);
    
    j=[];
    j=find(isnan(anos100min(1:n1,i)));
    min=anos100min(1:n1,i);
    min(j)=[];
    qmin1(i,1)=mean(min);
    
    j=[];
    j=find(isnan(anos100min(n1+1:end,i)));
    min=anos100min(n1+1:end,i);
    min(j)=[];
    qmin2(i,1)=mean(min);
end
%-------

%MEAN ALTERATION
changeAbsMinObs = [minis_ok qmin2-qmin1];
changeAbsMeanObs = [minis_ok qmean2-qmean1];
changeAbsMaxObs = [minis_ok qmax2-qmax1];

changePerMinObs = [minis_ok 100*changeAbsMinObs(:,2)./abs(qmin1)];
changePerMeanObs = [minis_ok 100*changeAbsMeanObs(:,2)./abs(qmean1)];
changePerMaxObs = [minis_ok 100*changeAbsMaxObs(:,2)./abs(qmax1)];
%---------------------------------------------

%------------- AGREEMENT BETWEEN ANA(OBSERVATION) AND MGB-SA (SIMULATION) ------------------
%MGB-SA ALTERATION IN SELECTED STATIONS
changePerMin=changePerMin(minis_ok);
changePerMean=changePerMean(minis_ok);
changePerMax=changePerMax(minis_ok);
%-----------

%ALTERATION AGREEMENT QMIN
for i=1:length(minis_ok)
    if (changePerMinObs(i,2)<-10 && changePerMin(i)<-10)
        agrMin(i) = -1;
    else
        if (changePerMinObs(i,2)>10 && changePerMin(i)>10)
            agrMin(i) = 1;
        else
            if abs(changePerMinObs(i,2))<=10 && abs(changePerMin(i))<=10
                agrMin(i) = 11;
            else
                if (changePerMinObs(i,2)<-10 && abs(changePerMin(i))<=10 && changePerMin(i)<0)||(changePerMin(i)<-10 && abs(changePerMinObs(i,2))<=10 && changePerMinObs(i,2)<0)
                    agrMin(i) = -2;
                else
                    if (changePerMinObs(i,2)>10 && abs(changePerMin(i))<=10 && changePerMin(i)>0)||(changePerMin(i)>10 && abs(changePerMinObs(i,2))<=10 && changePerMinObs(i,2)>0)
                        agrMin(i) = 2;
                    else
                        if (changePerMinObs(i,2)>10 && changePerMin(i)<-10) || (changePerMinObs(i,2)<-10 && changePerMin(i)>10)
                            agrMin(i) = 0;
                        else
                            agrMin(i) = 3;
                        end
                    end
                end
            end
        end
    end
end
agrMin = agrMin';
%--------------

%ALTERATION AGREEMENT QMEAN
for i=1:length(minis_ok)
    if (changePerMeanObs(i,2)<-10 && changePerMean(i)<-10)
        agrMean(i) = -1;
    else
        if (changePerMeanObs(i,2)>10 && changePerMean(i)>10)
            agrMean(i) = 1;
        else
            if abs(changePerMeanObs(i,2))<=10 && abs(changePerMean(i))<=10
                agrMean(i) = 11;
            else
                if (changePerMeanObs(i,2)<-10 && abs(changePerMean(i))<=10 && changePerMean(i)<0)||(changePerMean(i)<-10 && abs(changePerMeanObs(i,2))<=10 && changePerMeanObs(i,2)<0)
                    agrMean(i) = -2;
                else
                    if (changePerMeanObs(i,2)>10 && abs(changePerMean(i))<=10 && changePerMean(i)>0)||(changePerMean(i)>10 && abs(changePerMeanObs(i,2))<=10 && changePerMeanObs(i,2)>0)
                        agrMean(i) = 2;
                    else
                        if (changePerMeanObs(i,2)>10 && changePerMean(i)<-10) || (changePerMeanObs(i,2)<-10 && changePerMean(i)>10)
                            agrMean(i) = 0;
                        else
                            agrMean(i) = 3;
                        end
                    end
                end
            end
        end
    end
end
agrMean = agrMean';
%--------------

%ALTERATION AGREEMENT QMAX
for i=1:length(minis_ok)
    if (changePerMaxObs(i,2)<-10 && changePerMax(i)<-10)
        agrMax(i) = -1;
    else
        if (changePerMaxObs(i,2)>10 && changePerMax(i)>10)
            agrMax(i) = 1;
        else
            if abs(changePerMaxObs(i,2))<=10 && abs(changePerMax(i))<=10
                agrMax(i) = 11;
            else
                if (changePerMaxObs(i,2)<-10 && abs(changePerMax(i))<=10 && changePerMax(i)<0)||(changePerMax(i)<-10 && abs(changePerMaxObs(i,2))<=10 && changePerMaxObs(i,2)<0)
                    agrMax(i) = -2;
                else
                    if (changePerMaxObs(i,2)>10 && abs(changePerMax(i))<=10 && changePerMax(i)>0)||(changePerMax(i)>10 && abs(changePerMaxObs(i,2))<=10 && changePerMaxObs(i,2)>0)
                        agrMax(i) = 2;
                    else
                        if (changePerMaxObs(i,2)>10 && changePerMax(i)<-10) || (changePerMaxObs(i,2)<-10 && changePerMax(i)>10)
                            agrMax(i) = 0;
                        else
                            agrMax(i) = 3;
                        end
                    end
                end
            end
        end
    end
end
agrMax = agrMax';
%--------------------------------------------------------------------------------
% %PORCENTAGEM CONCORDÂNCIA
% %concordam
% concordanciaMin(1,1) = length(find(agrMin==1|agrMin==-1|agrMin==11))/length(minis_ok)*100;
% concordanciaMean(1,1) = length(find(agrMean==1|agrMean==-1|agrMean==11))/length(minis_ok)*100;
% concordanciaMax(1,1) = length(find(agrMax==1|agrMax==-1|agrMax==11))/length(minis_ok)*100;
% 
% %concordam parcialmente
% concordanciaMin(2,1) = length(find(agrMin==2|agrMin==-2))/length(minis_ok)*100;
% concordanciaMean(2,1) = length(find(agrMean==2|agrMean==-2))/length(minis_ok)*100;
% concordanciaMax(2,1) = length(find(agrMax==2|agrMax==-2))/length(minis_ok)*100;
% 
% %discordam parcialmente
% concordanciaMin(3,1) = length(find(agrMin==3))/length(minis_ok)*100;
% concordanciaMean(3,1) = length(find(agrMean==3))/length(minis_ok)*100;
% concordanciaMax(3,1) = length(find(agrMax==3))/length(minis_ok)*100;
% 
% %discordam
% concordanciaMin(4,1) = length(find(agrMin==0))/length(minis_ok)*100;
% concordanciaMean(4,1) = length(find(agrMean==0))/length(minis_ok)*100;
% concordanciaMax(4,1) = length(find(agrMax==0))/length(minis_ok)*100;
% %--------------------------------------------------------------------

save vazoes_obs1980_2020