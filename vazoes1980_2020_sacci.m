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

nMini = 33749;
load('vazoes1980_2020','changePerMean')
%---------------------------- SACCI'S MINIBASINS ------------------------------
rios = shaperead('vazao_geral2');

for i=1:length(rios)
    minis(i,1)=rios(i).Mini;
end
minis = sort(minis);
minis_unique = unique(minis); %eliminar valores repetidos e botar em ordem crescente

%------------------ COMPATIBILIZAR COM O SHP DO SACCI ---------------------
j = 1;
for i=minis_unique(1):nMini
    if minis_unique(j) == i
        ChangePerSacci(j,1) = changePerMean(i);
        j = j+1;
    else       
    end
end

for i=1:length(rios)
    for j=1:length(minis_unique)
        if minis(i)==minis_unique(j)
            ChangePerSacciShp(i,1)= ChangePerSacci(j);
        else
            j=j+1;
        end
    end
end

rios = struct2table(rios);
rios = sortrows(rios, 'Mini');
rios = table2struct(rios);

% ------------------------- AGREEMENT - RCP 4.5---------------------------
for i=1:length(rios)
    if rios(i).q_agree45==0
        agr_m45(i)=5;
    else
       if ChangePerSacciShp(i)<-10 && rios(i).q_mean45<-10
           agr_m45(i) = -1;
       else
           if ChangePerSacciShp(i)>10 && rios(i).q_mean45>10
               agr_m45(i) = 1;
           else
               if abs(ChangePerSacciShp(i))<=10 && abs(rios(i).q_mean45)<=10
                  agr_m45(i) = 11;
               else
                   if (ChangePerSacciShp(i)<-10 && abs(rios(i).q_mean45)<=10 && rios(i).q_mean45<0)||(rios(i).q_mean45<-10 && abs(ChangePerSacciShp(i))<=10 && ChangePerSacciShp(i)<0)
                       agr_m45(i) = -2;
                   else
                       if (ChangePerSacciShp(i)>10 && abs(rios(i).q_mean45)<=10 && rios(i).q_mean45>0)||(rios(i).q_mean45>10 && abs(ChangePerSacciShp(i))<=10 && ChangePerSacciShp(i)>0)
                           agr_m45(i) = 2;
                       else
                           if (ChangePerSacciShp(i)>10 && rios(i).q_mean45<-10) || (ChangePerSacciShp(i)<-10 && rios(i).q_mean45>10)
                               agr_m45(i) = 0;
                           else
                               agr_m45(i) = 3;
                           end
                       end
                   end
               end
           end
       end
    end
end
agr_m45 = table(agr_m45');
agr_m45.Properties.VariableNames={'agr_m45'};

% ------------------------- AGREEMENT - RCP 8.5---------------------------
for i=1:13240
    if rios(i).q_agree85==0
        agr_m85(i)=5;
    else
       if ChangePerSacciShp(i)<-10 && rios(i).q_mean85<-10
           agr_m85(i) = -1;
       else
           if ChangePerSacciShp(i)>10 && rios(i).q_mean85>10
               agr_m85(i) = 1;
           else
               if abs(ChangePerSacciShp(i))<=10 && abs(rios(i).q_mean85)<=10
                  agr_m85(i) = 11;
               else
                   if (ChangePerSacciShp(i)<-10 && abs(rios(i).q_mean85)<=10 && rios(i).q_mean85<0)||(rios(i).q_mean85<-10 && abs(ChangePerSacciShp(i))<=10 && ChangePerSacciShp(i)<0)
                       agr_m85(i) = -2;
                   else
                       if (ChangePerSacciShp(i)>10 && abs(rios(i).q_mean85)<=10 && rios(i).q_mean85>0)||(rios(i).q_mean85>10 && abs(ChangePerSacciShp(i))<=10 && ChangePerSacciShp(i)>0)
                           agr_m85(i) = 2;
                       else
                           if (ChangePerSacciShp(i)>10 && rios(i).q_mean85<-10) || (ChangePerSacciShp(i)<-10 && rios(i).q_mean85>10)
                               agr_m85(i) = 0;
                           else
                               agr_m85(i) = 3;
                           end
                       end
                   end
               end
           end
       end
    end
end
agr_m85 = table(agr_m85');
agr_m85.Properties.VariableNames={'agr_m85'};

%------------------------ ADICIONAR AO SHAPE ------------------------------
ChangePerSacciShp = table(ChangePerSacciShp);
ChangePerSacciShp.Properties.VariableNames={'changeMean'};
rios = struct2table(rios);
rios = [rios, ChangePerSacciShp, agr_m45, agr_m85];
rios=table2struct(rios);

save vazoes1980_2020_sacci

shapewrite(rios, 'sacci_geral');
