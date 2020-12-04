%% Metodo 1 con filtro fir (frecuencias normales)
clc
clear all
input_folder = 'D:\Documentos N\Tesis\BCICIV_2a_gdf'; %Nombre de la carpeta de donde se van a obtener los documentos
files = dir(fullfile(input_folder, '*.gdf')); % lectura de todos los documentos con una especifica extension, en nuestro caso gdf
file_paths = fullfile({files.folder}, {files.name});% it will specify images names with full path and extension
  
%% Diseño de Filtros
D_Filt1 = designfilt('bandpassfir','FilterOrder',1024, ....
    'StopbandFrequency1',(3),'PassbandFrequency1',(4), ...
    'PassbandFrequency2',(8),'StopbandFrequency2',(9), ...
    'SampleRate',250);   
D_Filt2 = designfilt('bandpassfir','FilterOrder',1024, ....
    'StopbandFrequency1',(7),'PassbandFrequency1',(8),...
    'PassbandFrequency2',(12),'StopbandFrequency2',(13),...
    'SampleRate',250);
D_Filt3 = designfilt('bandpassfir','FilterOrder',1024, ....
    'StopbandFrequency1',(11),'PassbandFrequency1',(12),...
    'PassbandFrequency2',(16),'StopbandFrequency2',(17),...
    'SampleRate',250);  
D_Filt4 = designfilt('bandpassfir','FilterOrder',1024, ....
    'StopbandFrequency1',(15),'PassbandFrequency1',(16),...
    'PassbandFrequency2',(20),'StopbandFrequency2',(21),...
    'SampleRate',250);   
D_Filt5 = designfilt('bandpassfir','FilterOrder',1024, ....
    'StopbandFrequency1',(19),'PassbandFrequency1',(20),...
    'PassbandFrequency2',(24),'StopbandFrequency2',(25),...
    'SampleRate',250);   
D_Filt6 = designfilt('bandpassfir','FilterOrder',1024, ....
    'StopbandFrequency1',(23),'PassbandFrequency1',(24),...
    'PassbandFrequency2',(28),'StopbandFrequency2',(29),...
    'SampleRate',250);    
D_Filt7 = designfilt('bandpassfir','FilterOrder',1024, ....
    'StopbandFrequency1',(27),'PassbandFrequency1',(28),...
    'PassbandFrequency2',(32),'StopbandFrequency2',(33),...
    'SampleRate',250);    
D_Filt8 = designfilt('bandpassfir','FilterOrder',1024, ....
    'StopbandFrequency1',(31),'PassbandFrequency1',(32),...
    'PassbandFrequency2',(36),'StopbandFrequency2',(37),...
    'SampleRate',250);    
D_Filt9 = designfilt('bandpassfir','FilterOrder',1024, ....
    'StopbandFrequency1',(35),'PassbandFrequency1',(36),...
    'PassbandFrequency2',(40),'StopbandFrequency2',(41),...
    'SampleRate',250);
%%
for i = 1 : 9
         
    [SIG,HEAD] = sload(file_paths{i}, '');%Cargar los archivos
     SIG(isnan(SIG))=0;
     SIGEEG=SIG(:,1:22);%Tomar los 22 canales del sujeto debido a que los ultmos 3 no son necesarios (señales EOG)
     % Encontrar movimientos 
     typ_l = find(HEAD.EVENT.TYP==769);%Encontrar el movimiento de mano izquierda (769)en la matriz head.event.typ
     pos_l=(HEAD.EVENT.POS(typ_l)); %Posiciones del evento 
     typ_r = find(HEAD.EVENT.TYP==770);%Encontrar el movimiento de mano derecha (770) en la matriz head.event.typ
     pos_r = (HEAD.EVENT.POS(typ_r));
     typ_f = find(HEAD.EVENT.TYP==771);
     pos_f = (HEAD.EVENT.POS(typ_f));
     typ_t = find(HEAD.EVENT.TYP==772);
     pos_t = (HEAD.EVENT.POS(typ_t));  
    
     %Obtener una matriz 3D con las señales filtradas, son 9 paginas debido
     %a las 9 bandas que se obtienen
     A= filtfilt(D_Filt1,SIGEEG);%Aplico el filtro 1
     A(:,:,2)= filtfilt(D_Filt2,SIGEEG);%Aplicacion del filtro 2
     A(:,:,3)= filtfilt(D_Filt3,SIGEEG);%Aplicacion del filtro 3
     A(:,:,4)= filtfilt(D_Filt4,SIGEEG);%Aplicacion del filtro 4
     A(:,:,5)= filtfilt(D_Filt5,SIGEEG);%Aplicacion del filtro 5
     A(:,:,6)= filtfilt(D_Filt6,SIGEEG);%Aplicacion del filtro 6
     A(:,:,7)= filtfilt(D_Filt7,SIGEEG);%Aplicacion del filtro 7
     A(:,:,8)= filtfilt(D_Filt8,SIGEEG);%Aplicacion del filtro 8
     A(:,:,9)= filtfilt(D_Filt9,SIGEEG);%Aplicacion del filtro 9

% Visualizacion de los filtros 

% fvtool(D_Filt1,'Analysis','freq');
% fvtool(D_Filt2,'Analysis','freq');
% fvtool(D_Filt3,'Analysis','freq');
% fvtool(D_Filt4,'Analysis','freq');
% fvtool(D_Filt5,'Analysis','freq');
% fvtool(D_Filt6,'Analysis','freq');
% fvtool(D_Filt7,'Analysis','freq');
% fvtool(D_Filt8,'Analysis','freq');
% fvtool(D_Filt9,'Analysis','freq');
% 
% fvtool(D_Filt1,'Analysis','freq');xlim([1 50])
% fvtool(D_Filt2,'Analysis','freq');xlim([1 50])
% fvtool(D_Filt3,'Analysis','freq');xlim([1 50])
% fvtool(D_Filt4,'Analysis','freq');xlim([1 50])
% fvtool(D_Filt5,'Analysis','freq');xlim([1 50])
% fvtool(D_Filt6,'Analysis','freq');xlim([1 50])
% fvtool(D_Filt7,'Analysis','freq');xlim([1 50])
% fvtool(D_Filt8,'Analysis','freq');xlim([1 50])
% fvtool(D_Filt9,'Analysis','freq');xlim([1 50])
% 
% fvtool(D_Filt1,'Analysis','freq');xlim([1 50]);hold on;plot(linspace(-3,-3),'r','linewidth',2);hold on;plot(linspace(-10,-10),'m','linewidth',2);
% fvtool(D_Filt2,'Analysis','freq');xlim([1 50]);hold on;plot(linspace(-3,-3),'r','linewidth',2);hold on;plot(linspace(-10,-10),'m','linewidth',2);
% fvtool(D_Filt3,'Analysis','freq');xlim([1 50]);hold on;plot(linspace(-3,-3),'r','linewidth',2);hold on;plot(linspace(-10,-10),'m','linewidth',2);
% fvtool(D_Filt4,'Analysis','freq');xlim([1 50]);hold on;plot(linspace(-3,-3),'r','linewidth',2);hold on;plot(linspace(-10,-10),'m','linewidth',2);
% fvtool(D_Filt5,'Analysis','freq');xlim([1 50]);hold on;plot(linspace(-3,-3),'r','linewidth',2);hold on;plot(linspace(-10,-10),'m','linewidth',2);
% fvtool(D_Filt6,'Analysis','freq');xlim([1 50]);hold on;plot(linspace(-3,-3),'r','linewidth',2);hold on;plot(linspace(-10,-10),'m','linewidth',2);
% fvtool(D_Filt7,'Analysis','freq');xlim([1 50]);hold on;plot(linspace(-3,-3),'r','linewidth',2);hold on;plot(linspace(-10,-10),'m','linewidth',2);
% fvtool(D_Filt8,'Analysis','freq');xlim([1 50]);hold on;plot(linspace(-3,-3),'r','linewidth',2);hold on;plot(linspace(-10,-10),'m','linewidth',2);
% fvtool(D_Filt9,'Analysis','freq');xlim([1 50]);hold on;plot(linspace(-3,-3),'r','linewidth',2);hold on;plot(linspace(-10,-10),'m','linewidth',2);

% Visualizacion señal original y filtradas
% 
% sig_o=SIG(93872:93872+1252,1);
% sig_1=A(93872:93872+1252,1,1);
% sig_2=A(93872:93872+1252,1,2);
% sig_3=A(93872:93872+1252,1,3);
% sig_4=A(93872:93872+1252,1,4);
% sig_5=A(93872:93872+1252,1,5);
% sig_6=A(93872:93872+1252,1,6);
% sig_7=A(93872:93872+1252,1,7);
% sig_8=A(93872:93872+1252,1,8);
% sig_9=A(93872:93872+1252,1,9);
% 
% plot(sig_o);hold on;
% plot(sig_1+50);hold on;
% plot(sig_2+100); hold on;
% plot(sig_3+150); hold on;
% plot(sig_4+200); hold on;
% plot(sig_5+250); hold on;
% plot(sig_6+300); hold on;
% plot(sig_7+350); hold on;
% plot(sig_8+400); hold on;
% plot(sig_9+450); 
% xlim([0 1300]);
% str = ['Señal Original y Señales Filtradas Sujeto ',num2str(i)];
% title(str);
% xlabel ('Tiempo(ms)');
% % ylabel ('Amplitud(mV)');
% xticks([0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300]);
% yticks([0 50 100 150 200 250 300 350 400 450]);
% yticklabels({'Señal Original','Señal Filtro (4-8)Hz',...
%     'Señal Filtro (8-12)Hz','Señal Filtro (12-16)Hz',...
%     'Señal Filtro (16-20)Hz','Señal Filtro (20-24)Hz',...
%     'Señal Filtro (20-28)Hz','Señal Filtro (28-32)Hz',...
%     'Señal Filtro (32-36)Hz','Señal Filtro (36-40)Hz'});
% ax = gca;
% ax.XGrid = 'on';
% ax.YGrid = 'off';

% subplot(2,1,1)
% plot(A(93872:93872+1252,1,1))
% hold on
% plot(SIG(93872:93872+1252,1))
% legend('Señal Filtrada (4-8)Hz', 'Señal original');
% str = ['Señal Original y Señal Filtrada (4-8)Hz Sujeto',num2str(i)];
% title(str)
% xlabel ('Tiempo(ms)');
% ylabel ('Amplitud(mV)');
% subplot(2,1,2)
% plot(A(93872:93872+1252,1,2))
% hold on
% plot(SIG(93872:93872+1252,1))
% legend('Señal Filtrada (8-12)Hz', 'Señal original');
% str = ['Señal Original y Señal Filtrada (8-12)Hz Sujeto',num2str(i)];
% title(str)
% xlabel ('Tiempo(ms)');
% ylabel ('Amplitud(mV)');

% CSP   
for j=1:72
        for k=1:9
        p1=(A(pos_l(j):pos_l(j)+1252,:,(k)));%movimientos de manos izquierda
        p2=(A(pos_r(j):pos_r(j)+1252,:,(k)));%movimientos de mano derecha
        p3=(A(pos_f(j):pos_f(j)+1252,:,(k)));%movimientos de pies
        p4=(A(pos_t(j):pos_t(j)+1252,:,(k)));%movimientos de lengua 
        S1 = cov(p1);% S1~[C x C] Covarianza de la matriz que contiene los movimientos de la mano izquierda
        S2 = cov(p2);% S1~[C x C] Covarianza de la matriz que contiene los movimientos de la mano derehca
        S3 = cov(p3);% S1~[C x C] Covarianza de la matriz que contiene los movimientos de los pies
        S4 = cov(p4);% S1~[C x C] Covarianza de la matriz que contiene los movimientos de la lengua
        [W1,L1] = eig(S1,S1+S2+S3+S4);   % Los filtros espaciales son las columnas de W
        [W2,L2] = eig(S2,S1+S2+S3+S4);
        [W3,L3] = eig(S3,S1+S2+S3+S4);
        [W4,L4] = eig(S4,S1+S2+S3+S4);
        lambda1 = diag(L1); % Eigenvalues
        lambda2 = diag(L2); % Eigenvalues
        lambda3 = diag(L3); % Eigenvalues
        lambda4 = diag(L4); % Eigenvalues
        s1=sort(lambda1,'descend');%Organizar los eigenvalues de mayor a menor del movimiento de mano izquierda
        s2=sort(lambda2,'descend');%Organizar los eigenvalues de mayor a menor del movimiento de mano derecha
        s3=sort(lambda3,'descend');%Organizar los eigenvalues de mayor a menor del movimiento de pies
        s4=sort(lambda4,'descend');%Organizar los eigenvalues de mayor a menor del movimiento de lengua
        %Eigenvalues maximos y minimos de cada matriz de movimiento 
        max_1(:,k)=s1(1); max2_1(:,k)=s1(2); min_1(:,k)=s1(22); min2_1(:,k)=s1(21);
        max_2(:,k)=s2(1); max2_2(:,k)=s2(2); min_2(:,k)=s2(22); min2_2(:,k)=s2(21);
        max_3(:,k)=s3(1); max2_3(:,k)=s3(2); min_3(:,k)=s3(22); min2_3(:,k)=s3(21);
        max_4(:,k)=s4(1); max2_4(:,k)=s4(2); min_4(:,k)=s4(22); min2_4(:,k)=s4(21);
        end        
     %Concatenar los maximos y minimos   
     V1(j,:)=horzcat(max_1,max2_1,min_1,min2_1);% maximos y minimos de mano izquierda
     V2(j,:)=horzcat(max_2,max2_2,min_2,min2_2);% maximos y minimos de mano derecha
     V3(j,:)=horzcat(max_3,max2_3,min_3,min2_3);% maximos y minimos de pies
     V4(j,:)=horzcat(max_4,max2_4,min_4,min2_4);% maximos y minimos de lengua
end
%Concatenar las matrices de maximos y minimos 
VTotal=horzcat(V1,V2,V3,V4);
VTransTotal=VTotal';
VTTotal=[V1;V2;V3;V4];

%Matriz de labels
Label1(1:72,1)=1;%Clase 1 es la clase de la mano izquierda
Label2(1:72,1)=2;%Clase 2 es la clase de la mano derecha
Label3(1:72,1)=3;%Clase 3 es la clase de los pies
Label4(1:72,1)=4;%Clase 4 es la clase de la lengua
labels=[Label1;Label2;Label3;Label4];

% Cross validation (train: 70%, test: 30%) CSP
cvy = cvpartition(labels,'HoldOut',0.3);
Parititon = cvy.test;
% Separate to training and test data
%Train data
TrainV = VTTotal(~Parititon,:);
TrainL = labels(~Parititon,:);
%test data
TestV = VTTotal(Parititon,:);
TestL  =labels(Parititon,:);

%CV movimiento 1 
CV1 = cvpartition(Label1,'HoldOut',0.3);
Parititon1 = CV1.test;
TrainV1 = V1(~Parititon1,:);
TrainL1 = Label1(~Parititon1,:);
TestV1 = V1(Parititon1,:);
TestL1  =Label1(Parititon1,:);

%CV movimiento 2
cvy2 = cvpartition(Label2,'HoldOut',0.3);
Parititon2 = cvy2.test;
TrainV2 = V2(~Parititon2,:);
TrainL2 = Label2(~Parititon2,:);
TestV2 = V2(Parititon2,:);
TestL2  =Label2(Parititon2,:);

%CV movimiento 3
cvy3 = cvpartition(Label3,'HoldOut',0.3);
Parititon3 = cvy3.test;
TrainV3 = V3(~Parititon3,:);
TrainL3 = Label3(~Parititon3,:);
TestV3 = V3(Parititon3,:);
TestL3  =Label3(Parititon3,:);

%CV movimiento 4
cvy4 = cvpartition(Label4,'HoldOut',0.3);
Parititon4 = cvy4.test;
TrainV4 = V4(~Parititon4,:);
TrainL4 = Label4(~Parititon4,:);
TestV4 = V4(Parititon4,:);
TestL4  =Label4(Parititon4,:);

%Concatenar todas las particiones y obtener las matrices de training y de
%test
TrainVT=[TrainV1;TrainV2;TrainV3;TrainV4];
TrainLT=[TrainL1;TrainL2;TrainL3;TrainL4];
TestVT=[TestV1;TestV2;TestV3;TestV4];
TestLT=[TestL1;TestL2;TestL3;TestL4];

% MIBIF
Nclass1=numel(pos_l); %Numero de repeticiones de la clase 1 (mano izq)
Nclass2=numel(pos_r); %Numero de repeticiones de la clase 2 (Mano der)
Nclass3=numel(pos_f); %Numero de repeticiones de la clase 3 (Pies)
Nclass4=numel(pos_t); %Numero de repeticiones de la clase 4 (Lengua)
Ntotal=numel(HEAD.Classlabel);%Numero total de repeticiones de todas las clases 
Pw=[Nclass1 Nclass2 Nclass3 Nclass4]/Ntotal;%Probabilidad de cada clase
Hw=-sum(Pw.*log2(Pw)); %Entropia

x1=TrainV1';%matriz de entrenamiento clase 1
x2=TrainV2';%matriz de entrenamiento clase 2
x3=TrainV3';%matriz de entrenamiento clase 3
x4=TrainV4';%matriz de entrenamiento clase 4

for n=1:144
    for m=1:51
        [f1,xi1,bw1]=ksdensity(x1(:,m),VTransTotal(n,m)); %Pfi,j_w de la clase 1  
        [f2,xi2,bw2]=ksdensity(x2(:,m),VTransTotal(n,m)); %Pfi,j_w de la clase 2 
        [f3,xi3,bw3]=ksdensity(x3(:,m),VTransTotal(n,m)); %Pfi,j_w de la clase 3 
        [f4,xi4,bw4]=ksdensity(x4(:,m),VTransTotal(n,m)); %Pfi,j_w de la clase 4  
        Pfji=sum(f1*Pw+f2*Pw+f3*Pw+f4*Pw);%sumatoria Pfj,i
        Pw1_fji=(f1*Pw)/Pfji;%Probabilidad Pw_fj,i usando regla de bayes de la clase 1
        Pw2_fji=(f2*Pw)/Pfji;%Probabilidad Pw_fj,i usando regla de bayes de la clase 2
        Pw3_fji=(f3*Pw)/Pfji;%Probabilidad Pw_fj,i usando regla de bayes de la clase 3
        Pw4_fji=(f4*Pw)/Pfji;%Probabilidad Pw_fj,i usando regla de bayes de la clase 4
        Hw_fj=-sum(Pw1_fji.*log2(Pw1_fji)+Pw2_fji.*log2(Pw2_fji)+...
    Pw3_fji.*log2(Pw3_fji)+Pw4_fji.*log2(Pw4_fji)); %Entropia condicional
MI(1:36,1:4)=(Hw-Hw_fj); %Mutual information
    end

end
% %4 valores maximos y sus respectivas posiciones de la matriz de MI 
[max1, pos1] = max(MI(:,1)); MI(pos1,1)= -Inf;[max2, pos2] = max(MI(:,1));MI(pos2,1)= -Inf;
[max3, pos3] = max(MI(:,1));MI(pos3,1)= -Inf;[max4, pos4] = max(MI(:,1));MI(pos4,1)= -Inf;
[max5, pos5] = max(MI(:,2));MI(pos5,2)= -Inf;[max6, pos6] = max(MI(:,2));MI(pos6,2)= -Inf;
[max7, pos7] = max(MI(:,2));MI(pos7,2)= -Inf;[max8, pos8] = max(MI(:,2));MI(pos8,2)= -Inf;
[max9, pos9] = max(MI(:,3));MI(pos9,3)= -Inf;[max10, pos10] = max(MI(:,3));MI(pos10,3)= -Inf;
[max11, pos11] = max(MI(:,3));MI(pos11,3)= -Inf;[max12, pos12] = max(MI(:,3));MI(pos12,3)= -Inf;
[max13, pos13] = max(MI(:,4));MI(pos13,4)= -Inf;[max14, pos14] = max(MI(:,4));MI(pos14,4)= -Inf;
[max15, pos15] = max(MI(:,4));MI(pos15,4)= -Inf;[max16, pos16] = max(MI(:,4));MI(pos16,4)= -Inf;

%Obtener las matrices respectivas en las posiciones anteriormente halladas
%matrices de training para el mutual information
Trainv1=[TrainV1(:,pos1),TrainV1(:,pos2),TrainV1(:,pos3),TrainV1(:,pos4)];
Trainv2=[TrainV2(:,pos4),TrainV2(:,pos6),TrainV2(:,pos7),TrainV2(:,pos8)];
Trainv3=[TrainV3(:,pos9),TrainV3(:,pos10),TrainV3(:,pos11),TrainV3(:,pos12)];
Trainv4=[TrainV4(:,pos13),TrainV4(:,pos14),TrainV4(:,pos15),TrainV4(:,pos16)];
Train2=[Trainv1;Trainv2;Trainv3;Trainv4];
%matrices de test para el mutual information
Testv2_1=[TestV1(:,pos1),TestV1(:,pos2),TestV1(:,pos3),TestV1(:,pos4)];
Testv2_2=[TestV2(:,pos4),TestV2(:,pos6),TestV2(:,pos7),TestV2(:,pos8)];
Testv2_3=[TestV3(:,pos9),TestV3(:,pos10),TestV3(:,pos11),TestV3(:,pos12)];
Testv2_4=[TestV4(:,pos13),TestV4(:,pos14),TestV4(:,pos15),TestV4(:,pos16)];
Test2v2_2=[Testv2_1;Testv2_2;Testv2_3;Testv2_4];
    
% LDA (MI)
ClassificationLDAMI = fitcdiscr(Train2,TrainLT,'DiscrimType', 'linear', 'Gamma', 0,...
    'FillCoeffs', 'off','ClassNames', [1; 2; 3; 4]);% crear el clasificador
ResultLDAMI = predict(ClassificationLDAMI,Test2v2_2);%predecir con el clasificador 
CVLDAMI = crossval(ClassificationLDAMI);% validación cruzada
genError1MI = kfoldLoss(CVLDAMI);% Error de generlización
CLDAMI(:,:,i)=confusionmat(TestLT,ResultLDAMI);%Matriz de confusion de todos los sujetos 
AccuracyLDAMI(i)=(1-genError1MI)*100;%Prediccion de todos los sujetos 
%sensibilidad
TPR1LDAMI = CLDAMI(1,1,i)/sum(CLDAMI(1,:,i));%True positive rate clase 1
TPR2LDAMI = CLDAMI(2,2,i)/sum(CLDAMI(2,:,i));%True positive rate clase 2
TPR3LDAMI = CLDAMI(3,3,i)/sum(CLDAMI(3,:,i));%True positive rate clase 3
TPR4LDAMI = CLDAMI(4,4,i)/sum(CLDAMI(4,:,i));%True positive rate clase 4
TPRTLDAMI (i)= (TPR1LDAMI+TPR2LDAMI+TPR3LDAMI+TPR4LDAMI)/4*100;%True positive rate total
%especificidad
 TNR11_LDAMI=(CLDAMI(2,2,i)+CLDAMI(2,3,i)+CLDAMI(2,4,i)+CLDAMI(3,2,i)+...
     CLDAMI(3,3,i)+CLDAMI(3,4,i)+CLDAMI(4,2,i)+CLDAMI(4,3,i)+CLDAMI(4,4,i));
 TNR22_LDAMI=(CLDAMI(1,1,i)+CLDAMI(1,3,i)+CLDAMI(1,4,i)+CLDAMI(3,1,i)+...
     CLDAMI(3,3,i)+CLDAMI(3,4,i)+CLDAMI(4,1,i)+CLDAMI(4,3,i)+CLDAMI(4,4,i));
 TNR33_LDAMI=((CLDAMI(1,1,i)+CLDAMI(1,2,i)+CLDAMI(1,4,i)+CLDAMI(2,1,i)+...
     CLDAMI(2,2,i)+CLDAMI(2,4,i)+CLDAMI(4,1,i)+CLDAMI(4,2,i)+CLDAMI(4,4,i)));
 TNR44_LDAMI=((CLDAMI(1,1,i)+CLDAMI(1,2,i)+CLDAMI(1,3,i)+CLDAMI(2,1,i)+...
     CLDAMI(2,2,i)+CLDAMI(2,3,i)+CLDAMI(3,1,i)+CLDAMI(3,2,i)+CLDAMI(3,3,i)));
TNR1LDAMI = TNR11_LDAMI/(TNR11_LDAMI+CLDAMI(2,1,i)+CLDAMI(3,1,i)+CLDAMI(4,1,i));%True negative rate clase 1
TNR2LDAMI = TNR22_LDAMI/(TNR22_LDAMI+CLDAMI(1,2,i)+CLDAMI(3,2,i)+CLDAMI(4,2,i));%True negative rate clase 2
TNR3LDAMI = TNR33_LDAMI/(TNR33_LDAMI+CLDAMI(1,3,i)+CLDAMI(2,3,i)+CLDAMI(4,3,i));%True negative rate clase 3
TNR4LDAMI = TNR44_LDAMI/(TNR44_LDAMI+CLDAMI(1,4,i)+CLDAMI(2,4,i)+CLDAMI(3,4,i));%True negative rate clase 4
TNRTLDAMI (i)= (TNR1LDAMI+TNR2LDAMI+TNR3LDAMI+TNR4LDAMI)/4*100;%True negative rate total
%Kappa LDAMI
N = length(TestLT);
CM_LDAMI=(CLDAMI(1,1,i)+CLDAMI(2,2,i)+CLDAMI(3,3,i)+CLDAMI(4,4,i));
C1corr_LDAMI=sum(CLDAMI(1,:,i));
C2corr_LDAMI=sum(CLDAMI(2,:,i));
C3corr_LDAMI=sum(CLDAMI(3,:,i));
C4corr_LDAMI=sum(CLDAMI(4,:,i));
C1pre_LDAMI=sum(CLDAMI(:,1,i));
C2pre_LDAMI=sum(CLDAMI(:,2,i));
C3pre_LDAMI=sum(CLDAMI(:,3,i));
C4pre_LDAMI=sum(CLDAMI(:,4,i));
TotalsumLDAMI(i)=(C1corr_LDAMI*C1pre_LDAMI+C2corr_LDAMI*C2pre_LDAMI+C3corr_LDAMI*C3pre_LDAMI+C4corr_LDAMI*C4pre_LDAMI);
Kappa_LDAMI(i)=(N*CM_LDAMI-(TotalsumLDAMI(i)))/(N^2-TotalsumLDAMI(i));

% ONE VS ALL SVM (MI)
template2 = templateSVM('KernelFunction', 'gaussian', 'PolynomialOrder', [], 'KernelScale', 'auto', ...
    'BoxConstraint', 1, 'Standardize', true);
classificationSVMMI = fitcecoc (Train2, TrainLT, 'Learners', template2, 'Coding', 'onevsall', ...
    'ClassNames', [1; 2; 3; 4]);
CVSVMMI = crossval(classificationSVMMI,'kfold',10); % validación cruzada
genError2MI = kfoldLoss(CVSVMMI); % Error de generlización
%Ones vs all prediccion
predictSVM = predict(classificationSVMMI,Test2v2_2);
% Confusion matrix generation
CSVMMI(:,:,i)= confusionmat(TestLT,predictSVM);
AccuracyMI(i)=(1-genError2MI)*100;

%sensibilidad
TPR1SVMMI = CSVMMI(1,1,i)/sum(CSVMMI(1,:,i));%True positive rate clase 1
TPR2SVMMI = CSVMMI(2,2,i)/sum(CSVMMI(2,:,i));
TPR3SVMMI = CSVMMI(3,3,i)/sum(CSVMMI(3,:,i));
TPR4SVMMI = CSVMMI(4,4,i)/sum(CSVMMI(4,:,i));
TPRTSVMMI (i)= (TPR1SVMMI+TPR2SVMMI+TPR3SVMMI+TPR4SVMMI)/4*100;
%especificidad
 TNR11_SVMMI=(CSVMMI(2,2,i)+CSVMMI(2,3,i)+CSVMMI(2,4,i)+CSVMMI(3,2,i)+...
     CSVMMI(3,3,i)+CSVMMI(3,4,i)+CSVMMI(4,2,i)+CSVMMI(4,3,i)+CSVMMI(4,4,i));
 TNR22_SVMMI=(CSVMMI(1,1,i)+CSVMMI(1,3,i)+CSVMMI(1,4,i)+CSVMMI(3,1,i)+...
     CSVMMI(3,3,i)+CSVMMI(3,4,i)+CSVMMI(4,1,i)+CSVMMI(4,3,i)+CSVMMI(4,4,i));
 TNR33_SVMMI=((CSVMMI(1,1,i)+CSVMMI(1,2,i)+CSVMMI(1,4,i)+CSVMMI(2,1,i)+...
     CSVMMI(2,2,i)+CSVMMI(2,4,i)+CSVMMI(4,1,i)+CSVMMI(4,2,i)+CSVMMI(4,4,i)));
 TNR44_SVMMI=((CSVMMI(1,1,i)+CSVMMI(1,2,i)+CSVMMI(1,3,i)+CSVMMI(2,1,i)+...
     CSVMMI(2,2,i)+CSVMMI(2,3,i)+CSVMMI(3,1,i)+CSVMMI(3,2,i)+CSVMMI(3,3,i)));
TNR1SVMMI = TNR11_SVMMI/(TNR11_SVMMI+CSVMMI(2,1,i)+CSVMMI(3,1,i)+CSVMMI(4,1,i));
TNR2SVMMI = TNR22_SVMMI/(TNR22_SVMMI+CSVMMI(1,2,i)+CSVMMI(3,2,i)+CSVMMI(4,2,i));
TNR3SVMMI = TNR33_SVMMI/(TNR33_SVMMI+CSVMMI(1,3,i)+CSVMMI(2,3,i)+CSVMMI(4,3,i));
TNR4SVMMI = TNR44_SVMMI/(TNR44_SVMMI+CSVMMI(1,4,i)+CSVMMI(2,4,i)+CSVMMI(3,4,i));
TNRTSVMMI (i)= (TNR1SVMMI+TNR2SVMMI+TNR3SVMMI+TNR4SVMMI)/4*100;
%Kappa SVMMI
N = length(TestLT);
CM_SVMMI=(CSVMMI(1,1,i)+CSVMMI(2,2,i)+CSVMMI(3,3,i)+CSVMMI(4,4,i));
C1corr_SVMMI=sum(CSVMMI(1,:,i));
C2corr_SVMMI=sum(CSVMMI(2,:,i));
C3corr_SVMMI=sum(CSVMMI(3,:,i));
C4corr_SVMMI=sum(CSVMMI(4,:,i));
C1pre_SVMMI=sum(CSVMMI(:,1,i));
C2pre_SVMMI=sum(CSVMMI(:,2,i));
C3pre_SVMMI=sum(CSVMMI(:,3,i));
C4pre_SVMMI=sum(CSVMMI(:,4,i));
TotalsumSVMMI(i)=(C1corr_SVMMI*C1pre_SVMMI+C2corr_SVMMI*C2pre_SVMMI+...
    C3corr_SVMMI*C3pre_SVMMI+C4corr_SVMMI*C4pre_SVMMI);
Kappa_SVMMI(i)=(N*CM_SVMMI-(TotalsumSVMMI(i)))/(N^2-TotalsumSVMMI(i));


% LDA (CSP)
ClassificationLDA = fitcdiscr(TrainVT,TrainLT,'DiscrimType', 'linear', 'Gamma', 0,...
    'FillCoeffs', 'off','ClassNames', [1; 2; 3; 4]);
ResultLDA = predict(ClassificationLDA,TestVT);
CVLDA = crossval(ClassificationLDA);
genError1 = kfoldLoss(CVLDA);
CLDA(:,:,i)=confusionmat(TestLT,ResultLDA);
AccuracyLDA(i)=(1-genError1)*100;
%sensibilidad
TPR1LDA = CLDA(1,1,i)/sum(CLDA(1,:,i));%True positive rate clase 1
TPR2LDA = CLDA(2,2,i)/sum(CLDA(2,:,i));
TPR3LDA = CLDA(3,3,i)/sum(CLDA(3,:,i));
TPR4LDA = CLDA(4,4,i)/sum(CLDA(4,:,i));
TPRTLDA (i)= (TPR1LDA+TPR2LDA+TPR3LDA+TPR4LDA)/4*100;
%especificidad
 TNR11LDA=(CLDA(2,2,i)+CLDA(2,3,i)+CLDA(2,4,i)+CLDA(3,2,i)+...
     CLDA(3,3,i)+CLDA(3,4,i)+CLDA(4,2,i)+CLDA(4,3,i)+CLDA(4,4,i));
 TNR22LDA=(CLDA(1,1,i)+CLDA(1,3,i)+CLDA(1,4,i)+CLDA(3,1,i)+...
     CLDA(3,3,i)+CLDA(3,4,i)+CLDA(4,1,i)+CLDA(4,3,i)+CLDA(4,4,i));
  TNR33LDA=((CLDA(1,1,i)+CLDA(1,2,i)+CLDA(1,4,i)+CLDA(2,1,i)+...
     CLDA(2,2,i)+CLDA(2,4,i)+CLDA(4,1,i)+CLDA(4,2,i)+CLDA(4,4,i)));
 TNR44LDA=((CLDA(1,1,i)+CLDA(1,2,i)+CLDA(1,3,i)+CLDA(2,1,i)+...
     CLDA(2,2,i)+CLDA(2,3,i)+CLDA(3,1,i)+CLDA(3,2,i)+CLDA(3,3,i)));
TNR1LDA = TNR11LDA/(TNR11LDA+CLDA(2,1,i)+CLDA(3,1,i)+CLDA(4,1,i));
TNR2LDA = TNR22LDA/(TNR22LDA+CLDA(1,2,i)+CLDA(3,2,i)+CLDA(4,2,i));
TNR3LDA = TNR33LDA/(TNR33LDA+CLDA(1,3,i)+CLDA(2,3,i)+CLDA(4,3,i));
TNR4LDA = TNR44LDA/(TNR44LDA+CLDA(1,4,i)+CLDA(2,4,i)+CLDA(3,4,i));
TNRTLDA (i)= (TNR1LDA+TNR2LDA+TNR3LDA+TNR4LDA)/4*100;
%Kappa LDA
N = length(TestLT);
CM_lda=(CLDA(1,1,i)+CLDA(2,2,i)+CLDA(3,3,i)+CLDA(4,4,i));
C1corr_lda=sum(CLDA(1,:,i));
C2corr_lda=sum(CLDA(2,:,i));
C3corr_lda=sum(CLDA(3,:,i));
C4corr_lda=sum(CLDA(4,:,i));
C1pre_lda=sum(CLDA(:,1,i));
C2pre_lda=sum(CLDA(:,2,i));
C3pre_lda=sum(CLDA(:,3,i));
C4pre_lda=sum(CLDA(:,4,i));
TotalsumLDA(i)=(C1corr_lda*C1pre_lda+C2corr_lda*C2pre_lda+C3corr_lda*C3pre_lda+C4corr_lda*C4pre_lda);
Kappa_lda(i)=(N*CM_lda-(TotalsumLDA(i)))/(N^2-TotalsumLDA(i))


% ONE VS ALL SVM (CSP)
template = templateSVM('KernelFunction', 'gaussian', 'PolynomialOrder', [], 'KernelScale', 'auto', ...
    'BoxConstraint', 1, 'Standardize', true);
classificationSVM = fitcecoc (TrainVT, TrainLT, 'Learners', template, 'Coding', 'onevsall', ...
    'ClassNames', [1; 2; 3; 4]);
CVSVM = crossval(classificationSVM,'kfold',10); % validación cruzada
genError2 = kfoldLoss(CVSVM); % Error de generlización
%Ones vs all prediccion
predictSVM = predict(classificationSVM,TestVT);
% Confusion matrix generation
CSVM(:,:,i)= confusionmat(TestLT,predictSVM);
Accuracy(i)=(1-genError2)*100;

%sensibilidad
TPR1SVM = CSVM(1,1,i)/sum(CSVM(1,:,i));%True positive rate clase 1
TPR2SVM = CSVM(2,2,i)/sum(CSVM(2,:,i));
TPR3SVM = CSVM(3,3,i)/sum(CSVM(3,:,i));
TPR4SVM = CSVM(4,4,i)/sum(CSVM(4,:,i));
TPRTSVM (i)= (TPR1SVM+TPR2SVM+TPR3SVM+TPR4SVM)/4*100;
%especificidad
 TNR11SVM=(CSVM(2,2,i)+CSVM(2,3,i)+CSVM(2,4,i)+CSVM(3,2,i)+...
     CSVM(3,3,i)+CSVM(3,4,i)+CSVM(4,2,i)+CSVM(4,3,i)+CSVM(4,4,i));
 TNR22SVM=(CSVM(1,1,i)+CSVM(1,3,i)+CSVM(1,4,i)+CSVM(3,1,i)+...
     CSVM(3,3,i)+CSVM(3,4,i)+CSVM(4,1,i)+CSVM(4,3,i)+CSVM(4,4,i));
  TNR33SVM=((CSVM(1,1,i)+CSVM(1,2,i)+CSVM(1,4,i)+CSVM(2,1,i)+...
     CSVM(2,2,i)+CSVM(2,4,i)+CSVM(4,1,i)+CSVM(4,2,i)+CSVM(4,4,i)));
 TNR44SVM=((CSVM(1,1,i)+CSVM(1,2,i)+CSVM(1,3,i)+CSVM(2,1,i)+...
     CSVM(2,2,i)+CSVM(2,3,i)+CSVM(3,1,i)+CSVM(3,2,i)+CSVM(3,3,i)));
TNR1SVM = TNR11SVM/(TNR11SVM+CSVM(2,1,i)+CSVM(3,1,i)+CSVM(4,1,i));
TNR2SVM = TNR22SVM/(TNR22SVM+CSVM(1,2,i)+CSVM(3,2,i)+CSVM(4,2,i));
TNR3SVM = TNR33SVM/(TNR33SVM+CSVM(1,3,i)+CSVM(2,3,i)+CSVM(4,3,i));
TNR4SVM = TNR44SVM/(TNR44SVM+CSVM(1,4,i)+CSVM(2,4,i)+CSVM(3,4,i));
TNRTSVM (i)= (TNR1SVM+TNR2SVM+TNR3SVM+TNR4SVM)/4*100;
%Kappa SVM
N = length(TestLT);
CM_SVM=(CSVM(1,1,i)+CSVM(2,2,i)+CSVM(3,3,i)+CSVM(4,4,i));
C1corr_SVM=sum(CSVM(1,:,i));
C2corr_SVM=sum(CSVM(2,:,i));
C3corr_SVM=sum(CSVM(3,:,i));
C4corr_SVM=sum(CSVM(4,:,i));
C1pre_SVM=sum(CSVM(:,1,i));
C2pre_SVM=sum(CSVM(:,2,i));
C3pre_SVM=sum(CSVM(:,3,i));
C4pre_SVM=sum(CSVM(:,4,i));
TotalsumSVM(i)=(C1corr_SVM*C1pre_SVM+C2corr_SVM*C2pre_SVM+C3corr_SVM*C3pre_SVM+C4corr_SVM*C4pre_SVM);
Kappa_SVM(i)=(N*CM_SVM-(TotalsumSVM(i)))/(N^2-TotalsumSVM(i));
end

%% matrices de confusion
 %visualizacion de las matrices de confusion de la implementacion csp
 for i=1:9
    figure;
    subplot(1,2,1)
    confusionchart(CSVM(:,:,i),'Title','Matriz de confusión SVM CSP','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
    subplot(1,2,2)
    confusionchart(CLDA(:,:,i),'Title','Matriz de confusión LDA CSP','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
 end
 %visualizacion de las matrices de confusion de la implementacion
 %informacion mutua 
for i=1:9
    figure;
    subplot(1,2,1)
    confusionchart(CSVMMI(:,:,i),'Title','Matriz de confusión SVM MI','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
    subplot(1,2,2)
    confusionchart(CLDAMI(:,:,i),'Title','Matriz de confusión LDA MI','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
end
 
%% Tablas metodos vs sujetos 
Total1=[AccuracyLDAMI,sum(AccuracyLDAMI)/9;TNRTLDAMI,sum(TNRTLDAMI)/9;...
    TPRTLDAMI,sum(TPRTLDAMI)/9;Kappa_LDAMI,sum(Kappa_LDAMI)/9];
Total2=[AccuracyMI,sum(AccuracyMI)/9;TNRTSVMMI,sum(TNRTSVMMI)/9;...
    TPRTSVMMI,sum(TPRTSVMMI)/9;Kappa_SVMMI,sum(Kappa_SVMMI)/9];
Total1=[AccuracyLDA,sum(AccuracyLDA)/9;TNRTLDA,sum(TNRTLDA)/9;...
    TPRTLDA,sum(TPRTLDA)/9;Kappa_lda,sum(Kappa_lda)/9];