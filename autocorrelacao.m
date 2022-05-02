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

%DISCHARGE AUTOCORRELATION CORRECTION THROUGH TREND FREE PRE-WHITENING (YUE ET AL., 2002)

load('vazoes1980_2020', 'qmin', 'qmean', 'qmax', 'sen_qmin', 'sen_qmean', 'sen_qmax')
nT = 40;

dif_zero_min = find(sen_qmin~=0);
dif_zero_mean = find(sen_qmean~=0);
dif_zero_max = find(sen_qmax~=0);

for i=1:nT
    qmin_adjust(i,:) = qmin(i,dif_zero_min) - sen_qmin(dif_zero_min)'.*i;
    qmean_adjust(i,:) = qmean(i,dif_zero_mean) - sen_qmean(dif_zero_mean)'.*i;
    qmax_adjust(i,:) = qmax(i,dif_zero_max) - sen_qmax(dif_zero_max)'.*i;
end

for i=1:length(qmin_adjust)
    [acf, lags, bounds] = autocorr(qmin_adjust(:,i), [], 1);
    acf_min(:,i) = acf;
    bound_min(:,i)=bounds(1);
end

acf = []; lags = []; bounds = [];

for i=1:length(qmean_adjust)
    [acf, lags, bounds] = autocorr(qmean_adjust(:,i), [], 1);
    acf_mean(:,i) = acf;
    bound_mean(:,i)=bounds(1);
end

acf = []; lags = []; bounds = [];

for i=1:length(qmax_adjust)
    [acf, lags, bounds] = autocorr(qmax_adjust(:,i), [], 1);
    acf_max(:,i) = acf;
    bound_max(:,i)=bounds(1);
end

autocorr_min = ones(21,length(qmin_adjust)).*0;
autocorr_mean = ones(21,length(qmean_adjust)).*0;
autocorr_max = ones(21,length(qmax_adjust)).*0;

for i=1:length(qmin_adjust)
    a=[];
    a = find(abs(acf_min(2:end,i))>bound_min(:,i));
    autocorr_min(a,i)=1;
end

for i=1:length(qmean_adjust)
    a=[];
    a = find(abs(acf_mean(2:end,i))>bound_mean(:,i));
    autocorr_mean(a,i)=1;
end

for i=1:length(qmax_adjust)
    a=[];
    a = find(abs(acf_max(2:end,i))>bound_max(:,i));
    autocorr_max(a,i)=1;
end

autocorr_min = autocorr_min(1,:).*acf_min(1,:);
autocorr_mean = autocorr_mean(1,:).*acf_mean(1,:);
autocorr_max = autocorr_max(1,:).*acf_max(1,:);

minis_acf_min = dif_zero_min(find(autocorr_min~=0));
minis_acf_mean = dif_zero_mean(find(autocorr_mean~=0));
minis_acf_max = dif_zero_max(find(autocorr_max~=0));

%Y't = Yt - phi.Yt-1
qmin_adjust2 = qmin_adjust;
qmean_adjust2 = qmean_adjust;
qmax_adjust2 = qmax_adjust;

for i=1:nT-1
    qmin_adjust2(i+1,:) = qmin_adjust(i+1,:) - autocorr_min.*qmin_adjust(i,:);
    qmean_adjust2(i+1,:) = qmean_adjust(i+1,:) - autocorr_mean.*qmean_adjust(i,:);
    qmax_adjust2(i+1,:) = qmax_adjust(i+1,:) - autocorr_max.*qmax_adjust(i,:);
end

%Y''t = Y't + b.t
for i=1:nT
    qmin_adjust3(i,:) = qmin_adjust2(i,:) + sen_qmin(dif_zero_min)'.*i;
    qmean_adjust3(i,:) = qmean_adjust2(i,:) + sen_qmean(dif_zero_mean)'.*i;
    qmax_adjust3(i,:) = qmax_adjust2(i,:) + sen_qmax(dif_zero_max)'.*i;
end

qminMK = qmin;
qmeanMK = qmean;
qmaxMK = qmax;

for i=1:length(minis_acf_min)
    qminMK(:,minis_acf_min(i)) = qmin_adjust3(:,i);
end

for i=1:length(minis_acf_mean)
    qmeanMK(:,minis_acf_mean(i)) = qmean_adjust3(:,i);
end

for i=1:length(minis_acf_max)
    qmaxMK(:,minis_acf_max(i)) = qmax_adjust3(:,i);
end

save autocorrelacao