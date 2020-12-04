%% Lectura de señales
[ST1, HT1] = sload('A01T.gdf'); %LECTURA signals train and header train subject 1
[ST2, HT2] = sload('A02T.gdf'); %LECTURA signals train and header train subject 2
[ST3, HT3] = sload('A03T.gdf'); %LECTURA signals train and header train subject 3
[ST4, HT4] = sload('A04T.gdf'); %LECTURA signals train and header train subject 4
[ST5, HT5] = sload('A05T.gdf'); %LECTURA signals train and header train subject 5
[ST6, HT6] = sload('A06T.gdf'); %LECTURA signals train and header train subject 6
[ST7, HT7] = sload('A07T.gdf'); %LECTURA signals train and header train subject 7
[ST8, HT8] = sload('A08T.gdf'); %LECTURA signals train and header train subject 8
[ST9, HT9] = sload('A09T.gdf'); %LECTURA signals train and header train subject 9

%%
pos1=HT1.EVENT.POS; %TIPO DE EVENTO DEL SUJETO 1
typ1=HT1.EVENT.TYP; %POSICION DEL EVENTO DEL SUJETO 1
dur1=HT1.EVENT.DUR; %DURACION DEL EVENTO DEL SUJETO 1
pos2=HE1.EVENT.POS;
typ2=HE1.EVENT.TYP;
dur2=HE1.EVENT.DUR;

%% VISUALIZACION DE SEÑALES
figure;
a =ST1(93872:94372,8); 
plot(a)
% hold on;
% plot(ST1(93872:94372,10))
% hold on
% plot(ST1(93872:94372,12))
% hold off
% 
% figure;
% plot(ST9(,8))
% hold on;
% plot(ST9(:,10))
% hold on
% plot(ST9(:,12))
% hold off
        
%% MAPAS DE COLOR CANALES C3,C4 Y CZ POR MOVIMIENTO (DE TODOS LOS 9 SUJETOS) 
clc
clear all
input_folder = 'D:\Documentos N\Tesis\BCICIV_2a_gdf'; %Nombre de la carpeta de donde se van a obtener los documentos
files = dir(fullfile(input_folder, '*.gdf')); % lectura de todos los documentos con una especifica extension, en nuestro caso gdf
file_paths = fullfile({files.folder}, {files.name});% it will specify images names with full path and extension

for i = 1 : numel(file_paths)
         
    [SIG,HEAD] = sload(file_paths{i}, ''); %Cargar todos los documentos y guardarlos en una matriz 
    typ_l = find(HEAD.EVENT.TYP==769);%Encontrar el movimiento de mano derecha (769)en la matriz head.event.typ
    pos_l=(HEAD.EVENT.POS(typ_l)); %
    typ_r = find(HEAD.EVENT.TYP==770);
    pos_r = (HEAD.EVENT.POS(typ_r));
    typ_f = find(HEAD.EVENT.TYP==771);
    pos_f = (HEAD.EVENT.POS(typ_f));
    typ_t = find(HEAD.EVENT.TYP==772);
    pos_t = (HEAD.EVENT.POS(typ_t));
    
     c = unique(HEAD.EVENT.TYP,'stable'); %contar cuantas veces hay un numero repetido
    for h=1:length(c)
    nv(h)=length( find(HEAD.EVENT.TYP==c(h)));
    end
            
    y = movimiento(SIG,pos_l);
    d = movimiento(SIG,pos_r);
    e = movimiento(SIG,pos_f);
    f = movimiento(SIG,pos_l);
    
    figure; 
    subplot( 2, 2,1);
    imagesc(y); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    title('Subplot 1: Mano izquierda')
    subplot( 2, 2,2 )
    imagesc(d);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    title('Subplot 2: Mano derecha')
    subplot( 2, 2,3 )
    imagesc(e);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    title('Subplot 3: Pies')
    subplot( 2, 2,4 )
    imagesc(f);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    title('Subplot 3: Lengua')
    
    
end
    
%% MAPAS DE COLOR POR BANDA DE FRECUENCIA Y CANAL (PROMEDIO ENTRE CANAL C3,C4 Y CZ
clc
clear all
input_folder = 'D:\Documentos N\Tesis\BCICIV_2a_gdf'; %Nombre de la carpeta de donde se van a obtener los documentos
files = dir(fullfile(input_folder, '*.gdf')); % lectura de todos los documentos con una especifica extension, en nuestro caso gdf
file_paths = fullfile({files.folder}, {files.name});% it will specify images names with full path and extension

for i = 1 : 9
         
    [SIG,HEAD] = sload(file_paths{i}, ''); %Cargar todos los documentos y guardarlos en una matriz 
    typ_l = find(HEAD.EVENT.TYP==769);%Encontrar el movimiento de mano derecha (769)en la matriz head.event.typ
    pos_l=(HEAD.EVENT.POS(typ_l)); %
    typ_r = find(HEAD.EVENT.TYP==770);
    pos_r = (HEAD.EVENT.POS(typ_r));
    typ_f = find(HEAD.EVENT.TYP==771);
    pos_f = (HEAD.EVENT.POS(typ_f));
    typ_t = find(HEAD.EVENT.TYP==772);
    pos_t = (HEAD.EVENT.POS(typ_t));   
    
    %Bandas de frecuencia delta, alfa y beta del movimiento de mano derecha
   
    Delta=FiltroDelta(SIG,pos_l);
    Alfa=FiltroAlfa(SIG,pos_l);
    Beta=FiltroBeta(SIG,pos_l);
    
    Delta2=FiltroDelta(SIG,pos_r);
    Alfa2=FiltroAlfa(SIG,pos_r);
    Beta2=FiltroBeta(SIG,pos_r);
    
    Delta3=FiltroDelta(SIG,pos_f);
    Alfa3=FiltroAlfa(SIG,pos_f);
    Beta3=FiltroBeta(SIG,pos_f);
    
    Delta4=FiltroDelta(SIG,pos_t);
    Alfa4=FiltroAlfa(SIG,pos_t);
    Beta4=FiltroBeta(SIG,pos_t);   
 
    figure; 
    subplot(3,1,1)
    imagesc(Delta); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Mano izquierda) Sujeto ',num2str(i)];
    title(str)
    subplot(3,1,2)
    imagesc(Alfa);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str2 = ['Banda de frecuencia Alfa (Mano izquierda) Sujeto ',num2str(i)];
    title(str2)
    subplot(3,1,3)
    imagesc(Beta);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str3 = ['Banda de frecuencia Beta (Mano izquierda) Sujeto ',num2str(i)];
    title(str3)
    
    figure;     
    subplot(3,1,1)
    imagesc(Delta2); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str4 = ['Banda de frecuencia Delta (Mano derecha) Sujeto ',num2str(i)];
    title(str4)
    subplot(3,1,2)
    imagesc(Alfa2);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str5 = ['Banda de frecuencia Alfa (Mano derecha) Sujeto ',num2str(i)];
    title(str5)
    subplot(3,1,3)
    imagesc(Beta2);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str6 = ['Banda de frecuencia Beta (Mano derecha) Sujeto ',num2str(i)];
    title(str6)

    figure;
    subplot(3,1,1)
    imagesc(Delta3); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str7 = ['Banda de frecuencia Delta (Pies) Sujeto ',num2str(i)];
    title(str7)
    subplot(3,1,2)
    imagesc(Alfa3);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str8 = ['Banda de frecuencia Alfa (Pies) Sujeto ',num2str(i)];
    title(str8)
    subplot(3,1,3)
    imagesc(Beta3);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str9 = ['Banda de frecuencia Beta (Pies) Sujeto ',num2str(i)];
    title(str9)
    
    figure;
    subplot(3,1,1)
    imagesc(Delta4); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str10 = ['Banda de frecuencia Delta (Lengua) Sujeto ',num2str(i)];
    title(str10)
    subplot(3,1,2)
    imagesc(Alfa4);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str11 = ['Banda de frecuencia Alfa (Lengua) Sujeto ',num2str(i)];
    title(str11)
    subplot(3,1,3)
    imagesc(Beta4);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str12 = ['Banda de frecuencia Beta (Lengua) Sujeto ',num2str(i)];
    title(str12)
    
end

%% MAPAS DE COLOR POR BANDA DE FRECUENCIA Y CANAL C3,C4 Y CZ;
clc
clear all
input_folder = 'D:\Documentos N\Tesis\BCICIV_2a_gdf'; %Nombre de la carpeta de donde se van a obtener los documentos
files = dir(fullfile(input_folder, '*.gdf')); % lectura de todos los documentos con una especifica extension, en nuestro caso gdf
file_paths = fullfile({files.folder}, {files.name});% it will specify images names with full path and extension

for i = 1 : 9
         
    [SIG,HEAD] = sload(file_paths{i}, ''); %Cargar todos los documentos y guardarlos en una matriz 
    typ_l = find(HEAD.EVENT.TYP==769);%Encontrar el movimiento de mano derecha (769)en la matriz head.event.typ
    pos_l=(HEAD.EVENT.POS(typ_l)); %
    typ_r = find(HEAD.EVENT.TYP==770);
    pos_r = (HEAD.EVENT.POS(typ_r));
    typ_f = find(HEAD.EVENT.TYP==771);
    pos_f = (HEAD.EVENT.POS(typ_f));
    typ_t = find(HEAD.EVENT.TYP==772);
    pos_t = (HEAD.EVENT.POS(typ_t));   
    
    %Bandas de frecuencia delta, alfa y beta del movimiento de mano derecha
    [Delta_C3L,Delta_C4L,Delta_CzL]=FiltroDelta2(SIG,pos_l);
    [Beta_C3L,Beta_C4L,Beta_CzL]=FiltroBeta2(SIG,pos_l);
    [Alfa_C3L,Alfa_C4L,Alfa_CzL]=FiltroAlfa2(SIG,pos_l);
    
    [Delta_C3R,Delta_C4R,Delta_CzR]=FiltroDelta2(SIG,pos_r);
    [Beta_C3R,Beta_C4R,Beta_CzR]=FiltroBeta2(SIG,pos_r);
    [Alfa_C3R,Alfa_C4R,Alfa_CzR]=FiltroAlfa2(SIG,pos_r);
    
    [Delta_C3F,Delta_C4F,Delta_CzF]=FiltroDelta2(SIG,pos_f);
    [Beta_C3F,Beta_C4F,Beta_CzF]=FiltroBeta2(SIG,pos_f);
    [Alfa_C3F,Alfa_C4F,Alfa_CzF]=FiltroAlfa2(SIG,pos_f);
    
    [Delta_C3T,Delta_C4T,Delta_CzT]=FiltroDelta2(SIG,pos_t);
    [Beta_C3T,Beta_C4T,Beta_CzT]=FiltroBeta2(SIG,pos_t);
    [Alfa_C3T,Alfa_C4T,Alfa_CzT]=FiltroAlfa2(SIG,pos_t);
    
    %Figura Mano izquierda
    figure; 
    subplot(3,3,1)
    imagesc(Delta_C3L); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Mano izq-C3) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,2)
    imagesc(Delta_C4L); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Mano izq-C4) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,3)
    imagesc(Delta_CzL); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Mano izq-Cz) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,4)
    imagesc(Alfa_C3L); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Mano izq-C3) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,5)
    imagesc(Alfa_C4L); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Mano izq-C4) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,6)
    imagesc(Alfa_CzL); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Mano izq-Cz) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,7)
    imagesc(Beta_C3L); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Beta (Mano izq-C3) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,8)
    imagesc(Beta_C4L); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Beta (Mano izq-C4) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,9)
    imagesc(Beta_CzL); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Beta (Mano izq-Cz) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    
    %Figura Mano derecha
    figure; 
    subplot(3,3,1)
    imagesc(Delta_C3R); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Mano der-C3) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,2)
    imagesc(Delta_C4R); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Mano der-C4) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,3)
    imagesc(Delta_CzR); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Mano der-Cz) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,4)
    imagesc(Alfa_C3R); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Mano der-C3) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,5)
    imagesc(Alfa_C4R); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Mano der-C4) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,6)
    imagesc(Alfa_CzR); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Mano der-Cz) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,7)
    imagesc(Beta_C3R); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
%     caxis([-8*10^(-5) 8*10^(-5)]);
    str = ['Banda de frecuencia Beta (Mano der-C3) Sujeto ',num2str(i)];
    title(str)
    subplot(3,3,8)
    imagesc(Beta_C4R); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
%     caxis([-8*10^(-5) 8*10^(-5)]);
    str = ['Banda de frecuencia Beta (Mano der-C4) Sujeto ',num2str(i)];
    title(str)    
    subplot(3,3,9)
    imagesc(Beta_CzR); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Beta (Mano der-Cz) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    
    %Figura Pies
    figure; 
    subplot(3,3,1)
    imagesc(Delta_C3F); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Pies-C3) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,2)
    imagesc(Delta_C4F); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Pies-C4) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,3)
    imagesc(Delta_CzF); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
%     caxis([-8*10^(-5) 8*10^(-5)]);
    str = ['Banda de frecuencia Delta (Pies-Cz) Sujeto ',num2str(i)];
    title(str)
    subplot(3,3,4)
    imagesc(Alfa_C3F); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Pies-C3) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,5)
    imagesc(Alfa_C4F); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Pies-C4) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,6)
    imagesc(Alfa_CzF); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Pies-Cz) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,7)
    imagesc(Beta_C3F); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
%     caxis([-8*10^(-5) 8*10^(-5)]);
    str = ['Banda de frecuencia Beta (Pies-C3) Sujeto ',num2str(i)];
    title(str)
    subplot(3,3,8)
    imagesc(Beta_C4F); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
%     caxis([-8*10^(-5) 8*10^(-5)]);
    str = ['Banda de frecuencia Beta (Pies-C4) Sujeto ',num2str(i)];
    title(str)    
    subplot(3,3,9)
    imagesc(Beta_CzF); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Beta (Pies-Cz) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);    
    
    %Figura Lengua
    figure; 
    subplot(3,3,1)
    imagesc(Delta_C3T); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Lengua-C3) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,2)
    imagesc(Delta_C4T); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Delta (Lengua-C4) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,3)
    imagesc(Delta_CzT); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
%     caxis([-8*10^(-5) 8*10^(-5)]);
    str = ['Banda de frecuencia Delta (Lengua-Cz) Sujeto ',num2str(i)];
    title(str)
    subplot(3,3,4)
    imagesc(Alfa_C3T); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Lengua-C3) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,5)
    imagesc(Alfa_C4T); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Lengua-C4) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,6)
    imagesc(Alfa_CzT); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Alfa (Lengua-Cz) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);
    subplot(3,3,7)
    imagesc(Beta_C3T); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
%     caxis([-8*10^(-5) 8*10^(-5)]);
    str = ['Banda de frecuencia Beta (Lengua-C3) Sujeto ',num2str(i)];
    title(str)
    subplot(3,3,8)
    imagesc(Beta_C4T); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
%     caxis([-8*10^(-5) 8*10^(-5)]);
    str = ['Banda de frecuencia Beta (Lengua-C4) Sujeto ',num2str(i)];
    title(str)    
    subplot(3,3,9)
    imagesc(Beta_CzT); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str = ['Banda de frecuencia Beta (Lengua-Cz) Sujeto ',num2str(i)];
    title(str)
%     caxis([-8*10^(-5) 8*10^(-5)]);  

end

%% MAPAS DE COLOR POR BANDA DE FRECUENCIA Y CANAL C3-C4 visualización por movimiento 
clc
clear all
input_folder = 'D:\Documentos N\Tesis\BCICIV_2a_gdf'; %Nombre de la carpeta de donde se van a obtener los documentos
files = dir(fullfile(input_folder, '*.gdf')); % lectura de todos los documentos con una especifica extension, en nuestro caso gdf
file_paths = fullfile({files.folder}, {files.name});% it will specify images names with full path and extension

for i = 1 : 9
         
    [SIG,HEAD] = sload(file_paths{i}, ''); %Cargar todos los documentos y guardarlos en una matriz 
    typ_l = find(HEAD.EVENT.TYP==769);%Encontrar el movimiento de mano derecha (769)en la matriz head.event.typ
    pos_l=(HEAD.EVENT.POS(typ_l)); %
    typ_r = find(HEAD.EVENT.TYP==770);
    pos_r = (HEAD.EVENT.POS(typ_r));
    typ_f = find(HEAD.EVENT.TYP==771);
    pos_f = (HEAD.EVENT.POS(typ_f));
    typ_t = find(HEAD.EVENT.TYP==772);
    pos_t = (HEAD.EVENT.POS(typ_t));   
    
    %Bandas de frecuencia delta, alfa y beta del movimiento de mano derecha
    [Delta_L]=FiltroDelta3(SIG,pos_l);
    [Beta_L]=FiltroBeta3(SIG,pos_l);
    [Alfa_L]=FiltroAlfa3(SIG,pos_l);
    
    [Delta_R]=FiltroDelta3(SIG,pos_r);
    [Beta_R]=FiltroBeta3(SIG,pos_r);
    [Alfa_R]=FiltroAlfa3(SIG,pos_r);
    
    [Delta_F]=FiltroDelta3(SIG,pos_f);
    [Beta_F]=FiltroBeta3(SIG,pos_f);
    [Alfa_F]=FiltroAlfa3(SIG,pos_f);
    
    [Delta_T]=FiltroDelta3(SIG,pos_t);
    [Beta_T]=FiltroBeta3(SIG,pos_t);
    [Alfa_T]=FiltroAlfa3(SIG,pos_t);
    
%Figura Mano izquierda   
    figure;     
    subplot(3,1,1)
    imagesc(Delta_L); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str1 = ['Banda de frecuencia Delta (Mano izquierda) C3-C4 Sujeto ',num2str(i)];
    title(str1)
    subplot(3,1,2)
    imagesc(Alfa_L);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str2 = ['Banda de frecuencia Alfa (Mano izquierda) C3-C4 Sujeto ',num2str(i)];
    title(str2)
    subplot(3,1,3)
    imagesc(Beta_L);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str3 = ['Banda de frecuencia Beta (Mano izquierda) C3-C4 Sujeto ',num2str(i)];
    title(str3)

%Figura Mano derecha   
    figure;     
    subplot(3,1,1)
    imagesc(Delta_R); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str4 = ['Banda de frecuencia Delta (Mano derecha) C3-C4 Sujeto ',num2str(i)];
    title(str4)
    subplot(3,1,2)
    imagesc(Alfa_R);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str5 = ['Banda de frecuencia Alfa (Mano derecha) C3-C4 Sujeto ',num2str(i)];
    title(str5)
    subplot(3,1,3)
    imagesc(Beta_R);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str6 = ['Banda de frecuencia Beta (Mano derecha) C3-C4 Sujeto ',num2str(i)];
    title(str6)
    
    %Figura pies 
    figure;
    subplot(3,1,1)
    imagesc(Delta_F); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str7 = ['Banda de frecuencia Delta (Pies) C3-C4 Sujeto ',num2str(i)];
    title(str7)
    subplot(3,1,2)
    imagesc(Alfa_F);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str8 = ['Banda de frecuencia Alfa (Pies) C3-C4 Sujeto ',num2str(i)];
    title(str8)
    subplot(3,1,3)
    imagesc(Beta_F);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str9 = ['Banda de frecuencia Beta (Pies) C3-C4 Sujeto ',num2str(i)];
    title(str9)
    
   % Figura lengua 
    figure;
    subplot(3,1,1)
    imagesc(Delta_T); 
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str10 = ['Banda de frecuencia Delta (Lengua) C3-C4 Sujeto ',num2str(i)];
    title(str10)
    subplot(3,1,2)
    imagesc(Alfa_T);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str11 = ['Banda de frecuencia Alfa (Lengua) C3-C4 Sujeto ',num2str(i)];
    title(str11)
    subplot(3,1,3)
    imagesc(Beta_T);
    xlabel ('Tiempo (s)');
    ylabel ('Número de repeticiones');
    colorbar;
    str12 = ['Banda de frecuencia Beta (Lengua) C3-C4 Sujeto ',num2str(i)];
    title(str12)
    
end

%% MAPAS DE COLOR POR BANDA DE FRECUENCIA Y CANAL C3-C4 visualización por banda de frecuencia
clc
clear all
input_folder = 'D:\Documentos N\Tesis\BCICIV_2a_gdf'; %Nombre de la carpeta de donde se van a obtener los documentos
files = dir(fullfile(input_folder, '*.gdf')); % lectura de todos los documentos con una especifica extension, en nuestro caso gdf
file_paths = fullfile({files.folder}, {files.name});% it will specify images names with full path and extension

for i = 1 : 9
         
    [SIG,HEAD] = sload(file_paths{i}, ''); %Cargar todos los documentos y guardarlos en una matriz 
    typ_l = find(HEAD.EVENT.TYP==769);%Encontrar el movimiento de mano derecha (769)en la matriz head.event.typ
    pos_l=(HEAD.EVENT.POS(typ_l)); %
    typ_r = find(HEAD.EVENT.TYP==770);
    pos_r = (HEAD.EVENT.POS(typ_r));
    typ_f = find(HEAD.EVENT.TYP==771);
    pos_f = (HEAD.EVENT.POS(typ_f));
    typ_t = find(HEAD.EVENT.TYP==772);
    pos_t = (HEAD.EVENT.POS(typ_t));   
    
    %Bandas de frecuencia delta, alfa y beta del movimiento de mano derecha
    [Delta_L]=FiltroDelta3(SIG,pos_l);
    [Beta_L]=FiltroBeta3(SIG,pos_l);
    [Alfa_L]=FiltroAlfa3(SIG,pos_l);
    
    [Delta_R]=FiltroDelta3(SIG,pos_r);
    [Beta_R]=FiltroBeta3(SIG,pos_r);
    [Alfa_R]=FiltroAlfa3(SIG,pos_r);
    
    [Delta_F]=FiltroDelta3(SIG,pos_f);
    [Beta_F]=FiltroBeta3(SIG,pos_f);
    [Alfa_F]=FiltroAlfa3(SIG,pos_f);
    
    [Delta_T]=FiltroDelta3(SIG,pos_t);
    [Beta_T]=FiltroBeta3(SIG,pos_t);
    [Alfa_T]=FiltroAlfa3(SIG,pos_t);
    
%Figura de banda de frecuencia Delta 
    figure;
    subplot(4,1,1)
    imagesc(Delta_L);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str1 = ['Banda de frecuencia Delta (Mano izquierda) C3-C4 Sujeto ',num2str(i)];
    title(str1)
    subplot(4,1,2)
    imagesc(Delta_R);
    xlabel ('Tiempo');
    ylabel ('Núm de repeticiones');
    colorbar;
    str2 = ['Banda de frecuencia Delta (Mano derecha) C3-C4 Sujeto ',num2str(i)];
    title(str2)    
    subplot(4,1,3)
    imagesc(Delta_F);
    xlabel ('Tiempo');
    ylabel ('Núm de repeticiones');
    colorbar;
    str3 = ['Banda de frecuencia Delta (Pies) C3-C4 Sujeto ',num2str(i)];
    title(str3)
    subplot(4,1,4)
    imagesc(Delta_T);
    xlabel ('Tiempo');
    ylabel ('Núm de repeticiones');
    colorbar;
    str4 = ['Banda de frecuencia Delta (Lengua) C3-C4 Sujeto ',num2str(i)];
    title(str4)
    
%Figura de banda de frecuencia Beta 
    figure;
    subplot(4,1,1)
    imagesc(Beta_L);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str5 = ['Banda de frecuencia Beta (Mano izquierda) C3-C4 Sujeto ',num2str(i)];
    title(str5)
    subplot(4,1,2)
    imagesc(Beta_R);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str6 = ['Banda de frecuencia Beta (Mano derecha) C3-C4 Sujeto ',num2str(i)];
    title(str6)    
    subplot(4,1,3)
    imagesc(Beta_F);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str7 = ['Banda de frecuencia Beta (Pies) C3-C4 Sujeto ',num2str(i)];
    title(str7)
    subplot(4,1,4)
    imagesc(Beta_T);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str8 = ['Banda de frecuencia Beta (Lengua) C3-C4 Sujeto ',num2str(i)];
    title(str8)

%Figura de banda de frecuencia Alfa 
    figure;
    subplot(4,1,1)
    imagesc(Alfa_L);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str9 = ['Banda de frecuencia Alfa (Mano izquierda) C3-C4 Sujeto ',num2str(i)];
    title(str9)
    subplot(4,1,2)
    imagesc(Alfa_R);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str10 = ['Banda de frecuencia Alfa (Mano derecha) C3-C4 Sujeto ',num2str(i)];
    title(str10)    
    subplot(4,1,3)
    imagesc(Alfa_F);
    xlabel ('Tiempo');
    ylabel ('Núm de repeticiones');
    colorbar;
    str11 = ['Banda de frecuencia Alfa (Pies) C3-C4 Sujeto ',num2str(i)];
    title(str11)
    subplot(4,1,4)
    imagesc(Alfa_T);
    xlabel ('Tiempo');
    ylabel ('Número de repeticiones');
    colorbar;
    str12 = ['Banda de frecuencia Alfa (Lengua) C3-C4 Sujeto ',num2str(i)];
    title(str12)
    
end

%% Metodo 1 con filtro kaiser (aplicando la funcion filtro_pasabanda)

clc
clear all
input_folder = 'D:\Documentos N\Tesis\BCICIV_2a_gdf'; %Nombre de la carpeta de donde se van a obtener los documentos
files = dir(fullfile(input_folder, '*.gdf')); % lectura de todos los documentos con una especifica extension, en nuestro caso gdf
file_paths = fullfile({files.folder}, {files.name});% it will specify images names with full path and extension

[Band1,h1]=filtro_pasabanda(3,4,8,9); %Banda 4-8Hz
[Band2,h2]=filtro_pasabanda(7,8,12,13); %Banda 8-12Hz
[Band3,h3]=filtro_pasabanda(11,12,16,17); %Banda 12-16Hz
[Band4,h4]=filtro_pasabanda(15,16,20,21);%Banda 16-20Hz
[Band5,h5]=filtro_pasabanda(19,20,24,25);%Banda 20-24Hz
[Band6,h6]=filtro_pasabanda(23,24,28,29);%Banda 24-28Hz
[Band7,h7]=filtro_pasabanda(27,28,32,33);%Banda 28-32Hz
[Band8,h8]=filtro_pasabanda(31,32,36,37);%Banda 32-36Hz
[Band9,h9]=filtro_pasabanda(35,36,40,41);%Banda 36-40Hz
for i = 1 : 1
         
    [SIG,HEAD] = sload(file_paths{i}, ''); %Cargar todos los documentos y guardarlos en una matriz 
    A(:,:,1)=filter(h1,1,SIG);
    A(:,:,2)=filter(h2,1,SIG);
    A(:,:,3)=filter(h3,1,SIG);
    A(:,:,4)=filter(h4,1,SIG);
    A(:,:,5)=filter(h5,1,SIG);
    A(:,:,6)=filter(h6,1,SIG);
    A(:,:,7)=filter(h7,1,SIG);
    A(:,:,8)=filter(h8,1,SIG);
    A(:,:,9)=filter(h9,1,SIG);
    
end

%% Metodo 1 con filtro fir (No es posible hacer el filtro 9, frecuencias normalizadas)

clc
clear all
input_folder = 'D:\Documentos N\Tesis\BCICIV_2a_gdf'; %Nombre de la carpeta de donde se van a obtener los documentos
files = dir(fullfile(input_folder, '*.gdf')); % lectura de todos los documentos con una especifica extension, en nuestro caso gdf
file_paths = fullfile({files.folder}, {files.name});% it will specify images names with full path and extension

f=2*pi/250;
    
D_Filt1 = designfilt('bandpassfir','StopbandFrequency1',(3*f),...
    'PassbandFrequency1',(4*f),'PassbandFrequency2',(8*f),...
    'StopbandFrequency2',(9*f),'DesignMethod','kaiserwin','ScalePassband',true);
    
D_Filt2 = designfilt('bandpassfir','StopbandFrequency1',(7*f),...
    'PassbandFrequency1',(8*f),'PassbandFrequency2',(12*f),...
    'StopbandFrequency2',(13*f),'DesignMethod','kaiserwin','ScalePassband',true);
    
D_Filt3 = designfilt('bandpassfir','StopbandFrequency1',(11*f),...
    'PassbandFrequency1',(12*f),'PassbandFrequency2',(16*f),...
    'StopbandFrequency2',(17*f),'DesignMethod','kaiserwin','ScalePassband',true);
    
D_Filt4 = designfilt('bandpassfir','StopbandFrequency1',(15*f),...
    'PassbandFrequency1',(16*f),'PassbandFrequency2',(20*f),...
    'StopbandFrequency2',(21*f),'DesignMethod','kaiserwin','ScalePassband',true);
    
D_Filt5 = designfilt('bandpassfir','StopbandFrequency1',(19*f),...
    'PassbandFrequency1',(20*f),'PassbandFrequency2',(24*f),...
    'StopbandFrequency2',(25*f),'DesignMethod','kaiserwin','ScalePassband',true);
    
D_Filt6 = designfilt('bandpassfir','StopbandFrequency1',(23*f),...
    'PassbandFrequency1',(24*f),'PassbandFrequency2',(28*f),...
    'StopbandFrequency2',(29*f),'DesignMethod','kaiserwin','ScalePassband',true);
    
D_Filt7 = designfilt('bandpassfir','StopbandFrequency1',(27*f),...
    'PassbandFrequency1',(28*f),'PassbandFrequency2',(32*f),...
    'StopbandFrequency2',(33*f),'DesignMethod','kaiserwin','ScalePassband',true);
    
D_Filt8 = designfilt('bandpassfir','StopbandFrequency1',(31*f),...
    'PassbandFrequency1',(32*f),'PassbandFrequency2',(36*f),...
    'StopbandFrequency2',(37*f),'DesignMethod','kaiserwin','ScalePassband',true);
    
% D_Filt9 = designfilt('bandpassfir','StopbandFrequency1',(35*f),...
%     'PassbandFrequency1',(36*f),'PassbandFrequency2',(40*f),...
%     'StopbandFrequency2',(41*f),'DesignMethod','kaiserwin','ScalePassband',true);

for i = 1 : 9
         
    [SIG,HEAD] = sload(file_paths{i}, '');
     SIG(isnan(SIG))=0;
     A(:,:,1)= filtfilt(D_Filt1,SIG);
     A(:,:,2)= filtfilt(D_Filt2,SIG);
     A(:,:,3)= filtfilt(D_Filt3,SIG);
     A(:,:,4)= filtfilt(D_Filt4,SIG);
     A(:,:,5)= filtfilt(D_Filt5,SIG);
     A(:,:,6)= filtfilt(D_Filt6,SIG);
     A(:,:,7)= filtfilt(D_Filt7,SIG);
     A(:,:,8)= filtfilt(D_Filt8,SIG);
%      A(:,:,9)= filtfilt(D_Filt9,s);
end


% y1 = filter(Filt1,s);
% 
% subplot(2,1,1)
% plot([y y1])
% title('Filtered Waveforms')
% legend('Zero-phase Filtering','Conventional Filtering')
% 
% subplot(2,1,2)
% plot(s)
% title('Original Waveform')

%% Metodo 1 con filtro fir (frecuencias normales)
clc
clear all
input_folder = 'D:\Documentos N\Tesis\BCICIV_2a_gdf'; %Nombre de la carpeta de donde se van a obtener los documentos
files = dir(fullfile(input_folder, '*.gdf')); % lectura de todos los documentos con una especifica extension, en nuestro caso gdf
file_paths = fullfile({files.folder}, {files.name});% it will specify images names with full path and extension
  
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

for i = 1 : 9
         
    [SIG,HEAD] = sload(file_paths{i}, '');
     SIG(isnan(SIG))=0;
     A(:,:,1)= filtfilt(D_Filt1,SIG);
     A(:,:,2)= filtfilt(D_Filt2,SIG);
     A(:,:,3)= filtfilt(D_Filt3,SIG);
     A(:,:,4)= filtfilt(D_Filt4,SIG);
     A(:,:,5)= filtfilt(D_Filt5,SIG);
     A(:,:,6)= filtfilt(D_Filt6,SIG);
     A(:,:,7)= filtfilt(D_Filt7,SIG);
     A(:,:,8)= filtfilt(D_Filt8,SIG);
     A(:,:,9)= filtfilt(D_Filt9,SIG);
end

% fvtool(D_Filt1);
% fvtool(D_Filt2);
% fvtool(D_Filt3);
% fvtool(D_Filt4);
% fvtool(D_Filt5);
% fvtool(D_Filt6);
% fvtool(D_Filt7);
% fvtool(D_Filt8);
% fvtool(D_Filt9);
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
% plot(sig_o); hold on;
% plot(sig_1+50); hold on;
% plot(sig_2+100); hold on;
% plot(sig_3+150); hold on;
% plot(sig_4+200); hold on;
% plot(sig_5+250); hold on;
% plot(sig_6+300); hold on;
% plot(sig_7+350); hold on;
% plot(sig_8+400); hold on;
% plot(sig_9+450); hold on;
% str = ['Señal Original y Señales Filtradas Sujeto ',num2str(i)];
% title(str)
% xlabel ('Tiempo(ms)');
% ylabel ('Amplitud(mV)');

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

%     typ_l = find(HEAD.EVENT.TYP==769);%Encontrar el movimiento de mano derecha (769)en la matriz head.event.typ
%     pos_l=(HEAD.EVENT.POS(typ_l)); %
%     typ_r = find(HEAD.EVENT.TYP==770);
%     pos_r = (HEAD.EVENT.POS(typ_r));
%     typ_f = find(HEAD.EVENT.TYP==771);
%     pos_f = (HEAD.EVENT.POS(typ_f));
%     typ_t = find(HEAD.EVENT.TYP==772);
%     pos_t = (HEAD.EVENT.POS(typ_t));  
% 
%     %CSP
%     
% for j=1:72
%         for k=1:9
%         p1=(A(pos_l(j):pos_l(j)+1252,:,(k)));
%         p2=(A(pos_r(j):pos_r(j)+1252,:,(k)));
%         p3=(A(pos_f(j):pos_f(j)+1252,:,(k)));
%         p4=(A(pos_t(j):pos_t(j)+1252,:,(k)));
%         S1 = cov(p1);% S1~[C x C]
%         S2 = cov(p2);
%         S3 = cov(p3);
%         S4 = cov(p4);
%         [W1,L1] = eig(S1,S1+S2+S3+S4);   % Mixing matrix W (spatial filters are columns)
%         [W2,L2] = eig(S2,S1+S2+S3+S4);
%         [W3,L3] = eig(S3,S1+S2+S3+S4);
%         [W4,L4] = eig(S4,S1+S2+S3+S4);% Mixing matrix W (spatial filters are columns)
%         lambda1 = diag(L1); 
%         lambda2 = diag(L2); 
%         lambda3 = diag(L3); 
%         lambda4 = diag(L4); % Eigenvalues
%         s1=sort(lambda1,'descend');
%         s2=sort(lambda2,'descend');
%         s3=sort(lambda3,'descend');
%         s4=sort(lambda4,'descend');
%         max_1(:,k)=s1(1); max2_1(:,k)=s1(2); min_1(:,k)=s1(25); min2_1(:,k)=s1(24);
%         max_2(:,k)=s2(1); max2_2(:,k)=s2(2); min_2(:,k)=s2(25); min2_2(:,k)=s2(24);
%         max_3(:,k)=s3(1); max2_3(:,k)=s3(2); min_3(:,k)=s3(25); min2_3(:,k)=s3(24);
%         max_4(:,k)=s4(1); max2_4(:,k)=s4(2); min_4(:,k)=s4(25); min2_4(:,k)=s4(24);
%         end        
%      V1(j,:)=horzcat(max_1,max2_1,min_1,min2_1);
%      V2(j,:)=horzcat(max_2,max2_2,min_2,min2_2);
%      V3(j,:)=horzcat(max_3,max2_3,min_3,min2_3);
%      V4(j,:)=horzcat(max_4,max2_4,min_4,min2_4);
% end
% 
% VTotal=horzcat(V1,V2,V3,V4);
% VTransTotal=VTotal';
% VT1=V1';
% VT2=V2';
% VT3=V3';
% VT4=V4';
% 
% Nclass1=numel(pos_l);
% Nclass2=numel(pos_r);
% Nclass3=numel(pos_f);
% Nclass4=numel(pos_t);
% Ntotal=numel(HEAD.Classlabel);
% Px=[Nclass1 Nclass2 Nclass3 Nclass4]/Ntotal;
% 
% a(1:36,1)=1;
% a2(37:72,1)=2;
% a3(73:108,1)= 3;
% a4(109:144,1)= 4;
% Y=[a(1:36,1);a2(37:72,1);a3(73:108,1);a4(109:144,1)]';
% 
% target=VTransTotal(Y==1,:);
% nontarget=VTransTotal(Y~=1,:);
% 










