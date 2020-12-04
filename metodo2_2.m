input_folder = 'D:\Documentos N\Tesis\BCICIV_2a_gdf'; %Nombre de la carpeta de donde se van a obtener los documentos
files = dir(fullfile(input_folder, '*.gdf')); % lectura de todos los documentos con una especifica extension, en nuestro caso gdf
file_paths = fullfile({files.folder}, {files.name});% it will specify images names with full path and extension

D_Filt1 = designfilt('bandpassfir','FilterOrder',1024, ....
    'StopbandFrequency1',(6),'PassbandFrequency1',(8), ...
    'PassbandFrequency2',(40),'StopbandFrequency2',(41), ...
    'SampleRate',250);

%Inicializacion de marices necesarias para el codigo
cross1(1:2505,1:87)=0;cross2(1:2505,1:87)=0;
cross3(1:2505,1:87)=0;cross4(1:2505,1:87)=0;
Accuracy(1:9)=0;TPRT(1:9)=0;TNRT(1:9)=0;Kappa_svm(1:9)=0;
TotalsumSVM(1:9)=0;
Accuracy2(1:9)=0;TPRT_2(1:9)=0;TNRT_2(1:9)=0;Kappa_svm_2(1:9)=0;
TotalsumSVM_2(1:9)=0;CSVM2(2,2,9)=0;
Accuracy3(1:9)=0;TPRT_3(1:9)=0;TNRT_3(1:9)=0;Kappa_svm_3(1:9)=0;
TotalsumSVM_3(1:9)=0;CSVM3(2,2,9)=0;
Accuracy4(1:9)=0;TPRT_4(1:9)=0;TNRT_4(1:9)=0;Kappa_svm_4(1:9)=0;
TotalsumSVM_4(1:9)=0;CSVM4(2,2,9)=0;
TotalAccuracySVM(1:9)=0;TotalSensibilidadSVM(1:9)=0;
TotalEspecificidadSVM(1:9)=0;TotalKappaSVM(1:9)=0;
AccuracyLDA(1:9)=0;TPRTLDA(1:9)=0;TNRTLDA(1:9)=0;Kappa_lda(1:9)=0;
TotalsumLDA(1:9)=0;CLDA(2,2,9)=0;
AccuracyLDA2(1:9)=0;TPRTLDA_2(1:9)=0;TNRTLDA_2(1:9)=0;Kappa_lda_2(1:9)=0;
TotalsumLDA_2(1:9)=0;CLDA2(2,2,9)=0;
AccuracyLDA3(1:9)=0;TPRTLDA_3(1:9)=0;TNRTLDA_3(1:9)=0;Kappa_lda_3(1:9)=0;
TotalsumLDA_3(1:9)=0;CLDA3(2,2,9)=0;
AccuracyLDA4(1:9)=0;TPRTLDA_4(1:9)=0;TNRTLDA_4(1:9)=0;Kappa_lda_4(1:9)=0;
TotalsumLDA_4(1:9)=0;CLDA4(2,2,9)=0;Mean1(1:72,1:87)=0;
Median1(1:72,1:87)=0;Mode1(1:72,1:87)=0;
Stan_dev1(1:72,1:87)=0;Max1(1:72,1:87)=0;Min1(1:72,1:87)=0;
Mean2(1:72,1:87)=0;Median2(1:72,1:87)=0;Mode2(1:72,1:87)=0;
Stan_dev2(1:72,1:87)=0;Max2(1:72,1:87)=0;Min2(1:72,1:87)=0;
Mean3(1:72,1:87)=0;Median3(1:72,1:87)=0;Mode3(1:72,1:87)=0;
Stan_dev3(1:72,1:87)=0;Max3(1:72,1:87)=0;Min3(1:72,1:87)=0;
Mean4(1:72,1:87)=0;Median4(1:72,1:87)=0;Mode4(1:72,1:87)=0;
Stan_dev4(1:72,1:87)=0;Max4(1:72,1:87)=0;Min4(1:72,1:87)=0;

for i = 1 
    
    [SIG,HEAD] = sload(file_paths{i}, '');
     SIG(isnan(SIG))=0;
     SIGEEG=SIG(:,1:22);
     A=filtfilt(D_Filt1,SIGEEG);  

     typ_l = (HEAD.EVENT.TYP==769);%Encontrar el movimiento de mano derecha (769)en la matriz head.event.typ
     pos_l=(HEAD.EVENT.POS(typ_l)); %
     typ_r = (HEAD.EVENT.TYP==770);
     pos_r = (HEAD.EVENT.POS(typ_r));
     typ_f = (HEAD.EVENT.TYP==771);
     pos_f = (HEAD.EVENT.POS(typ_f));
     typ_t = (HEAD.EVENT.TYP==772);
     pos_t = (HEAD.EVENT.POS(typ_t));  
    
for j=1:72
   p1=(A(pos_l(j):pos_l(j)+1252,(1:22)));%matriz de repeticiones del movimiento de mano izq
   p2=(A(pos_r(j):pos_r(j)+1252,(1:22)));%matriz de repeticiones del movimiento de mano der
   p3=(A(pos_f(j):pos_f(j)+1252,(1:22)));%matriz de repeticiones del movimiento de pies
   p4=(A(pos_t(j):pos_t(j)+1252,(1:22)));%matriz de repeticiones del movimiento de lengua
   p5=[p1(:,2:22),p2,p3,p4];
   p6=[p2(:,2:22),p1,p3,p4];
   p7=[p3(:,2:22),p1,p2,p4];
   p8=[p4(:,2:22),p1,p2,p3]; 

    cross1=fliplr(xcorr2(p1(:,1),p5(:,:)));%correlacion cruzada entre el fp1 y el resto de movimiento 
    cross2=fliplr(xcorr2(p2(:,1),p6(:,:)));%correlacion cruzada entre el segundo movimiento y el resto de movimiento
    cross3=fliplr(xcorr2(p3(:,1),p7(:,:)));%correlacion cruzada entre el tercer movimiento y el resto de movimiento
    cross4=fliplr(xcorr2(p4(:,1),p8(:,:)));%correlacion cruzada entre el cuarto movimiento y el resto de movimiento

    %Obtener las caracteristicas siguientes de cada correlacion cruzada 
    %moda, mediana, media, desviacion estandar, valor maximo y minimo
    mode_1=sort(cross1);mode_1=mode_1(1,:);Mode1(j,:)=mode_1;%obtencion de la moda de la correlacion cruzada
    Mean1(j,:)=mean(cross1);Median1(j,:)=median(cross1);
    Stan_dev1(j,:)=std(cross1);Max1(j,:)=max(cross1);Min1(j,:)=min(cross1);
    
    mode_2=sort(cross2);mode_2=mode_2(1,:);Mode2(j,:)=mode_2;%obtencion de la moda de la correlacion cruzada
    Mean2(j,:)=mean(cross2);Median2(j,:)=median(cross2);
    Stan_dev2(j,:)=std(cross2);Max2(j,:)=max(cross2);Min2(j,:)=min(cross2);
    
    mode_3=sort(cross3);mode_3=mode_3(1,:);Mode3(j,:)=mode_3;%obtencion de la moda de la correlacion cruzada
    Mean3(j,:)=mean(cross3);Median3(j,:)=median(cross3);
    Stan_dev3(j,:)=std(cross3);Max3(j,:)=max(cross3);Min3(j,:)=min(cross3);
    
    mode_4=sort(cross4);mode_4=mode_4(1,:);Mode4(j,:)=mode_4;%obtencion de la moda de la correlacion cruzada
    Mean4(j,:)=mean(cross4);Median4(j,:)=median(cross4);
    Stan_dev4(j,:)=std(cross4);Max4(j,:)=max(cross4);Min4(j,:)=min(cross4);
    
end

%matriz de caracteristicas teniendo en cuenta las caracteristicas que se
%obtuvieron en el apartado anterior
features1=[Mean1(:,1:21),Median1(:,1:21),Mode1(:,1:21),Stan_dev1(:,1:21),...
    Max1(:,1:21),Min1(:,1:21),Median1(:,5:10);Mean1(:,22:43),Median1(:,22:43),Mode1(:,22:43),Stan_dev1(:,22:43),...
    Max1(:,22:43),Min1(:,22:43);Mean1(:,44:65),Median1(:,44:65),Mode1(:,44:65),Stan_dev1(:,44:65),...
    Max1(:,44:65),Min1(:,44:65);Mean1(:,66:87),Median1(:,66:87),Mode1(:,66:87),Stan_dev1(:,66:87),...
    Max1(:,66:87),Min1(:,66:87)];
features2=[Mean2(:,1:21),Median2(:,1:21),Mode2(:,1:21),Stan_dev2(:,1:21),...
    Max2(:,1:21),Min2(:,1:21),Median2(:,5:10);Mean2(:,22:43),Median2(:,22:43),Mode2(:,22:43),Stan_dev2(:,22:43),...
    Max2(:,22:43),Min2(:,22:43);Mean2(:,44:65),Median2(:,44:65),Mode2(:,44:65),Stan_dev2(:,44:65),...
    Max2(:,44:65),Min2(:,44:65);Mean2(:,66:87),Median2(:,66:87),Mode2(:,66:87),Stan_dev2(:,66:87),...
    Max2(:,66:87),Min2(:,66:87)];
features3=[Mean3(:,1:21),Median3(:,1:21),Mode3(:,1:21),Stan_dev3(:,1:21),...
    Max3(:,1:21),Min3(:,1:21),Median3(:,5:10);Mean3(:,22:43),Median3(:,22:43),Mode3(:,22:43),Stan_dev3(:,22:43),...
    Max3(:,22:43),Min3(:,22:43);Mean3(:,44:65),Median3(:,44:65),Mode3(:,44:65),Stan_dev3(:,44:65),...
    Max3(:,44:65),Min3(:,44:65);Mean3(:,66:87),Median3(:,66:87),Mode3(:,66:87),Stan_dev3(:,66:87),...
    Max3(:,66:87),Min3(:,66:87)];
features4=[Mean4(:,1:21),Median4(:,1:21),Mode4(:,1:21),Stan_dev4(:,1:21),...
    Max4(:,1:21),Min4(:,1:21),Median4(:,5:10);Mean4(:,22:43),Median4(:,22:43),Mode4(:,22:43),Stan_dev4(:,22:43),...
    Max4(:,22:43),Min4(:,22:43);Mean4(:,44:65),Median4(:,44:65),Mode4(:,44:65),Stan_dev4(:,44:65),...
    Max4(:,44:65),Min4(:,44:65);Mean4(:,66:87),Median4(:,66:87),Mode4(:,66:87),Stan_dev4(:,66:87),...
    Max4(:,66:87),Min4(:,66:87)];

%Matriz de etiquetas
Label1(1:72,1)=1;
Label2(1:216,1)=2;
labels=[Label1;Label2];

% Cross validation (train: 70%, test: 30%) CSP
cvy = cvpartition(labels,'Holdout',0.3);
Parititon = cvy.test;
% Separate to training and test data
% %Train data
TrainF1 = features1(~Parititon,:);
TrainF2 = features2(~Parititon,:);
TrainF3 = features3(~Parititon,:);
TrainF4 = features4(~Parititon,:);
TrainL = labels(~Parititon,:);
%test data
TestF1 = features1(Parititon,:);
TestF2 = features2(Parititon,:);
TestF3 = features3(Parititon,:);
TestF4 = features4(Parititon,:);
TestL  = labels(Parititon,:);
% 
% SVM Clase A vs todas
template = templateSVM('KernelFunction', 'gaussian', 'PolynomialOrder', [], 'KernelScale', 'auto', ...
    'BoxConstraint', 1, 'Standardize', true);
classificationSVM = fitcecoc (TrainF1,TrainL,'Learners', template, 'Coding', 'onevsall', ...
    'ClassNames', [1;2]);
CVSVM = crossval(classificationSVM,'kfold',10); % validación cruzada
res=kfoldPredict(CVSVM);
genError2 = kfoldLoss(CVSVM); % Error de generlización
%Ones vs all prediccion
predictSVM = predict(classificationSVM,TestF1);
%Confusion matrix generation
CSVM(:,:,i)= confusionmat(TestL,predictSVM);
Accuracy(i)=(1-genError2)*100;
%sensibilidad
TPRT(i)=(CSVM(1,1,i)/(CSVM(1,1,i)+CSVM(1,2,i)))*100;
%especificidad
TNRT(i)=(CSVM(2,2,i)/(CSVM(2,2,i)+CSVM(2,1,i)))*100;
%Kappa
N = length(TestL);
CM_svm=(CSVM(1,1,i)+CSVM(2,2,i));
C1corr_svm=sum(CSVM(1,:,i));
C2corr_svm=sum(CSVM(2,:,i));
C1pre_svm=sum(CSVM(:,1,i));
C2pre_svm=sum(CSVM(:,2,i));
TotalsumSVM(i)=(C1corr_svm*C1pre_svm+C2corr_svm*C2pre_svm);
Kappa_svm(i)=(N*CM_svm-(TotalsumSVM(i)))/(N^2-TotalsumSVM(i));

% SVM Clase B vs todas
template2 = templateSVM('KernelFunction', 'gaussian', 'PolynomialOrder', [], 'KernelScale', 'auto', ...
    'BoxConstraint', 1, 'Standardize', true);
classificationSVM2 = fitcecoc (TrainF2,TrainL,'Learners', template, 'Coding', 'onevsall', ...
    'ClassNames', [1;2]);
CVSVM2 = crossval(classificationSVM2,'kfold',10); % validación cruzada
genError22 = kfoldLoss(CVSVM2); % Error de generlización
%Ones vs all prediccion
predictSVM2 = predict(classificationSVM2,TestF2);
%Confusion matrix generation
CSVM2(:,:,i)= confusionmat(TestL,predictSVM2);
Accuracy2(i)=(1-genError22)*100;

%sensibilidad
TPRT_2(i)=(CSVM2(1,1,i)/(CSVM2(1,1,i)+CSVM2(1,2,i)))*100;
%especificidad
TNRT_2(i)=(CSVM2(2,2,i)/(CSVM2(2,2,i)+CSVM2(2,1,i)))*100;

%Kappa
CM_svm_2=(CSVM2(1,1,i)+CSVM2(2,2,i));
C1corr_svm_2=sum(CSVM2(1,:,i));
C2corr_svm_2=sum(CSVM2(2,:,i));
C1pre_svm_2=sum(CSVM2(:,1,i));
C2pre_svm_2=sum(CSVM2(:,2,i));
TotalsumSVM_2(i)=(C1corr_svm_2*C1pre_svm_2+C2corr_svm_2*C2pre_svm_2);
Kappa_svm_2(i)=(N*CM_svm_2-(TotalsumSVM_2(i)))/(N^2-TotalsumSVM_2(i));

% SVM Clase C vs todas
template3 = templateSVM('KernelFunction', 'gaussian', 'PolynomialOrder', [], 'KernelScale', 'auto', ...
    'BoxConstraint', 1, 'Standardize', true);
classificationSVM3 = fitcecoc (TrainF3,TrainL,'Learners', template3, 'Coding', 'onevsall', ...
    'ClassNames', [1;2]);
CVSVM3 = crossval(classificationSVM3,'kfold',10); % validación cruzada
genError23 = kfoldLoss(CVSVM3); % Error de generlización
%Ones vs all prediccion
predictSVM3 = predict(classificationSVM3,TestF3);
%Confusion matrix generation
[C2, ~] = confusionmat(TestL,predictSVM3);
CSVM3(:,:,i)= confusionmat(TestL,predictSVM3);
Accuracy3(i)=(1-genError23)*100;

%sensibilidad
 TPRT_3(i)=(CSVM3(1,1,i)/(CSVM3(1,1,i)+CSVM3(1,2,i)))*100;
%especificidad
TNRT_3(i)=(CSVM3(2,2,i)/(CSVM3(2,2,i)+CSVM3(2,1,i)))*100;
%Kappa
CM_svm_3=(CSVM3(1,1,i)+CSVM3(2,2,i));
C1corr_svm_3=sum(CSVM3(1,:,i));
C2corr_svm_3=sum(CSVM3(2,:,i));
C1pre_svm_3=sum(CSVM3(:,1,i));
C2pre_svm_3=sum(CSVM3(:,2,i));
TotalsumSVM_3(i)=(C1corr_svm_3*C1pre_svm_3+C2corr_svm_3*C2pre_svm_3);
Kappa_svm_3(i)=(N*CM_svm_3-(TotalsumSVM_3(i)))/(N^2-TotalsumSVM_3(i));

% SVM Clase D vs todas
template4 = templateSVM('KernelFunction', 'gaussian', 'PolynomialOrder', [], 'KernelScale', 'auto', ...
    'BoxConstraint', 1, 'Standardize', true);
classificationSVM4 = fitcecoc (TrainF4,TrainL,'Learners', template4, 'Coding', 'onevsone', ...
    'ClassNames', [1;2]);
CVSVM4 = crossval(classificationSVM4,'kfold',10); % validación cruzada
genError24 = kfoldLoss(CVSVM4); % Error de generlización
%Ones vs all prediccion
predictSVM4 = predict(classificationSVM4,TestF4);
%Confusion matrix generation
CSVM4(:,:,i)= confusionmat(TestL,predictSVM4);
Accuracy4(i)=(1-genError24)*100;

%sensibilidad
TPRT_4(i)=(CSVM4(1,1,i)/(CSVM4(1,1,i)+CSVM(1,2,i)))*100;
% %especificidad
TNRT_4(i)=(CSVM4(2,2,i)/(CSVM4(2,2,i)+CSVM4(2,1,i)))*100;

%Kappa
CM_svm_4=(CSVM4(1,1,i)+CSVM4(2,2,i));
C1corr_svm_4=sum(CSVM4(1,:,i));
C2corr_svm_4=sum(CSVM4(2,:,i));
C1pre_svm_4=sum(CSVM4(:,1,i));
C2pre_svm_4=sum(CSVM4(:,2,i));
TotalsumSVM_4(i)=(C1corr_svm_4*C1pre_svm_4+C2corr_svm_4*C2pre_svm_4);
Kappa_svm_4(i)=(N*CM_svm_4-(TotalsumSVM_4(i)))/(N^2-TotalsumSVM_4(i));

TotalAccuracySVM(i)=(Accuracy(i)+Accuracy2(i)+Accuracy3(i)+Accuracy4(i))/4;
TotalSensibilidadSVM(i)=(TPRT(i)+TPRT_2(i)+TPRT_3(i)+TPRT_4(i))/4;
TotalEspecificidadSVM(i)=(TNRT(i)+TNRT_2(i)+TNRT_3(i)+TNRT_4(i))/4;
TotalKappaSVM(i)=(Kappa_svm(i)+Kappa_svm_2(i)+Kappa_svm_3(i)+Kappa_svm_4(i))/4;

%LDA A versus todas
ClassificationLDA = fitcdiscr(TrainF1,TrainL,'DiscrimType', 'linear', 'Gamma', 0,...
    'FillCoeffs', 'off','ClassNames', [1; 2]);
ResultLDA = predict(ClassificationLDA,TestF1);
CVLDA = crossval(ClassificationLDA);
genError1 = kfoldLoss(CVLDA);
AccuracyLDA(i)=(1-genError1)*100;
CLDA(:,:,i)= confusionmat(TestL,ResultLDA);

% %sensibilidad
TPRTLDA(i)=(CLDA(1,1,i)/(CLDA(1,1,i)+CLDA(1,2,i)))*100;
% %especificidad
TNRTLDA(i)=(CLDA(2,2,i)/(CLDA(2,2,i)+CLDA(2,1,i)))*100;
%Kappa
CM_lda=(CLDA(1,1,i)+CLDA(2,2,i));
C1corr_lda=sum(CLDA(1,:,i));
C2corr_lda=sum(CLDA(2,:,i));
C1pre_lda=sum(CLDA(:,1,i));
C2pre_lda=sum(CLDA(:,2,i));
TotalsumLDA(i)=(C1corr_lda*C1pre_lda+C2corr_lda*C2pre_lda);
Kappa_lda(i)=(N*CM_lda-(TotalsumLDA(i)))/(N^2-TotalsumLDA(i));

%LDA b versus todas
ClassificationLDA2 = fitcdiscr(TrainF2,TrainL,'DiscrimType', 'linear', 'Gamma', 0,...
    'FillCoeffs', 'off','ClassNames', [1; 2]);
ResultLDA2 = predict(ClassificationLDA2,TestF2);
CVLDA2 = crossval(ClassificationLDA2);
genError12 = kfoldLoss(CVLDA2);
AccuracyLDA2(i)=(1-genError12)*100;
CLDA2(:,:,i)= confusionmat(TestL,ResultLDA2);

% %sensibilidad
TPRTLDA_2(i)=(CLDA2(1,1,i)/(CLDA(1,1,i)+CLDA(1,2,i)))*100;
% % %especificidad
TNRTLDA_2(i)=(CLDA2(2,2,i)/(CLDA2(2,2,i)+CLDA2(2,1,i)))*100;
%Kappa
CM_lda_2=(CLDA2(1,1,i)+CLDA2(2,2,i));
C1corr_lda_2=sum(CLDA2(1,:,i));
C2corr_lda_2=sum(CLDA2(2,:,i));
C1pre_lda_2=sum(CLDA2(:,1,i));
C2pre_lda_2=sum(CLDA2(:,2,i));
TotalsumLDA_2(i)=(C1corr_lda_2*C1pre_lda_2+C2corr_lda_2*C2pre_lda_2);
Kappa_lda_2(i)=(N*CM_lda_2-(TotalsumLDA_2(i)))/(N^2-TotalsumLDA_2(i));

%LDA C versus todas
ClassificationLDA3 = fitcdiscr(TrainF3,TrainL,'DiscrimType', 'linear', 'Gamma', 0,...
    'FillCoeffs', 'off','ClassNames', [1; 2]);
ResultLDA3 = predict(ClassificationLDA3,TestF3);
CVLDA3 = crossval(ClassificationLDA3);
genError13 = kfoldLoss(CVLDA3);
CLDA3(:,:,i)=confusionmat(TestL,ResultLDA3);
AccuracyLDA3(i)=(1-genError13)*100;

%sensibilidad
TPRTLDA_3(i)=(CLDA3(1,1,i)/(CLDA3(1,1,i)+CLDA3(1,2,i)))*100;
% % %especificidad
TNRTLDA_3(i)=(CLDA3(2,2,i)/(CLDA3(2,2,i)+CLDA3(2,1,i)))*100;
%Kappa
CM_lda_3=(CLDA3(1,1,i)+CLDA3(2,2,i));
C1corr_lda_3=sum(CLDA3(1,:,i));
C2corr_lda_3=sum(CLDA3(2,:,i));
C1pre_lda_3=sum(CLDA3(:,1,i));
C2pre_lda_3=sum(CLDA3(:,2,i));
TotalsumLDA_3(i)=(C1corr_lda_3*C1pre_lda_3+C2corr_lda_3*C2pre_lda_3);
Kappa_lda_3(i)=(N*CM_lda_3-(TotalsumLDA_3(i)))/(N^2-TotalsumLDA_3(i));

%LDA D versus todas
ClassificationLDA4 = fitcdiscr(TrainF4,TrainL,'DiscrimType', 'linear', 'Gamma', 0,...
    'FillCoeffs', 'off','ClassNames', [1; 2]);
ResultLDA4 = predict(ClassificationLDA4,TestF4);
CVLDA4 = crossval(ClassificationLDA4);
genError14 = kfoldLoss(CVLDA4);
AccuracyLDA4(i)=(1-genError14)*100;
CLDA4(:,:,i)= confusionmat(TestL,ResultLDA4);

%sensibilidad
TPRTLDA_4(i)=(CLDA4(1,1,i)/(CLDA4(1,1,i)+CLDA4(1,2,i)))*100;
% % %especificidad
TNRTLDA_4(i)=(CLDA4(2,2,i)/(CLDA4(2,2,i)+CLDA4(2,1,i)))*100;
%Kappa
CM_lda_4=(CLDA4(1,1,i)+CLDA4(2,2,i));
C1corr_lda_4=sum(CLDA4(1,:,i));
C2corr_lda_4=sum(CLDA4(2,:,i));
C1pre_lda_4=sum(CLDA4(:,1,i));
C2pre_lda_4=sum(CLDA4(:,2,i));
TotalsumLDA_4(i)=(C1corr_lda_4*C1pre_lda_4+C2corr_lda_4*C2pre_lda_4);
Kappa_lda_4(i)=(N*CM_lda_4-(TotalsumLDA_4(i)))/(N^2-TotalsumLDA_4(i));
end
%%
% m = max(max(abs(p1)));
% canales = size(p1,2);
% 
% plot(p1 + (ones(size(p1,1),1)*(1:canales) - 1) * m * 1.5);
% 
% Dx=200;
% y1=min(p1(:,1));
% y2=max(p1(:,22)+22*80);
% for t=1:20:1052
% for k=1:22
% plot(p1(:,k)+k*80)
% hold on;
% end
% xlim([t t+Dx]);
% ylim([y1 y2]);
% drawnow;
% end

%% Tablas de los metodos vs sujetos

 AvsTodas=[Accuracy,sum(Accuracy)/9;TNRT,sum(TNRT)/9;TPRT,sum(TPRT)/9;...
     Kappa_svm,sum(Kappa_svm)/9];
 AvsTodasLDA=[AccuracyLDA,sum(AccuracyLDA)/9;TNRTLDA,sum(TNRTLDA)/9;...
     TPRTLDA,sum(TPRTLDA)/9;Kappa_lda,sum(Kappa_lda)/9];
 BvsTodas=[Accuracy2,sum(Accuracy2)/9;TNRT_2,sum(TNRT_2)/9;...
     TPRT_2,sum(TPRT_2)/9;Kappa_svm_2,sum(Kappa_svm_2)/9];
 BvsTodasLDA=[AccuracyLDA2,sum(AccuracyLDA2)/9;TNRTLDA_2,...
     sum(TNRTLDA_2)/9;TPRTLDA_2,sum(TPRTLDA_2)/9;Kappa_lda_2,sum(Kappa_lda_2)/9];
  CvsTodas=[Accuracy3,sum(Accuracy3)/9;TNRT_3,sum(TNRT_3)/9;TPRT_3,...
      sum(TPRT_3)/9;Kappa_svm_3,sum(Kappa_svm_3)/9];
  CvsTodasLDA=[AccuracyLDA3,sum(AccuracyLDA3)/9;TNRTLDA_3,sum(TNRTLDA_3)/9;...
      TPRTLDA_3,sum(TPRTLDA_3)/9;Kappa_lda_3,sum(Kappa_lda_3)/9];
    DvsTodas=[Accuracy4,sum(Accuracy4)/9;TNRT_4,sum(TNRT_4)/9;TPRT_4,...
      sum(TPRT_4)/9;Kappa_svm_4,sum(Kappa_svm_4)/9];
  DvsTodasLDA=[AccuracyLDA4,sum(AccuracyLDA4)/9;TNRTLDA_4,sum(TNRTLDA_4)/9;...
      TPRTLDA_4,sum(TPRTLDA_4)/9;Kappa_lda_4,sum(Kappa_lda_4)/9];
  
  %% Visualizacion de las matrices de confusion dentro de una sola figura 
  
   for i=1:9
    figure;
    subplot(2,2,1)
    confusionchart(CSVM(:,:,i),'Title','Matriz de confusión SVM A vs Todas','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
    subplot(2,2,2)
    confusionchart(CSVM2(:,:,i),'Title','Matriz de confusión SVM B vs Todas','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
    subplot(2,2,3)
    confusionchart(CSVM3(:,:,i),'Title','Matriz de confusión SVM C vs Todas','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
    subplot(2,2,4)
    confusionchart(CSVM4(:,:,i),'Title','Matriz de confusión SVM D vs Todas','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
   end
 
      for i=1:9
    figure;
    subplot(2,2,1)
    confusionchart(CLDA(:,:,i),'Title','Matriz de confusión LDA A vs Todas','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
    subplot(2,2,2)
    confusionchart(CLDA2(:,:,i),'Title','Matriz de confusión LDA B vs Todas','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
    subplot(2,2,3)
    confusionchart(CLDA3(:,:,i),'Title','Matriz de confusión LDA C vs Todas','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
    subplot(2,2,4)
    confusionchart(CLDA4(:,:,i),'Title','Matriz de confusión LDA D vs Todas','RowSummary','row-normalized',...
    'ColumnSummary','column-normalized','FontSize',12);
      end
 
%%
Totalpredict=[predictSVM,predictSVM2,predictSVM3,predictSVM4];