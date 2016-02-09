%% Programa de Intersecção Espacial Directa
%Realizado por Rui Jorge Abrunhosa Nunes nº32092
%FOTOGRAMETRIA ANALÍTICA 2011/2012
%MESTRADO EM ENGENHARIA GEOGRÁFICA
%--------------------------------------------------------------------------
clc;
close all;
format long g;

%% Abertura de Ficheiros
resultados = fopen('resultados.dat','w');

%% Leitura dos dados dos ficheiros e guardar dados em vectores-------------
%Leitura dos dados de ficheiros
DELIMITER = ' ';
HEADERLINES = 1;
newData1 = importdata('cfoto.dat', DELIMITER, HEADERLINES);
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
coordenadas = data;

newData2 = importdata('oriext.dat');
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData2.(vars{i}));
end
ori_externa = data;

%Criar vectores independentes para cada variável
    %Coordenadas de pontos nas fotos
    x1 = coordenadas(:,2)/1000;
    y1 = coordenadas(:,3)/1000;
    x2 = coordenadas(:,4)/1000;
    y2 = coordenadas(:,5)/1000;

    %Parâmetros externos
    fi = ori_externa(:,2);
    omega = ori_externa(:,3);
    kappa = ori_externa(:,4);
    X0 = ori_externa(:,5);
    Y0 = ori_externa(:,6);
    Z0 = ori_externa(:,7);
    
    %Identificações
    id_ponto = coordenadas(:,1);
    id_foto = ori_externa(:,1);

    %Parâmetros internos
    x0 = 0.002/1000;
    y0 = 0.004/1000;
    c = 153.070/1000;

%Escolha dos parâmetros internos das fotos que nos interessam
num_foto1 = input('Qual o numero da primeira fotografia?: ');
for i = 1:size(id_foto,1)
    if (id_foto(i,:) == num_foto1)
        fi_foto1 = fi(i,1) ;
        omega_foto1 = omega(i,1) ;
        kappa_foto1 = kappa(i,1) ;
        X0_foto1 = X0(i,1) ;
        Y0_foto1 = Y0(i,1) ;
        Z0_foto1 = Z0(i,1) ;
    end
end

num_foto2 = input('Qual o numero da segunda fotografia?:  ');
for j = 1:size(id_foto,1)
    if (id_foto(j,:) == num_foto2)
        fi_foto2 = fi(j,1) ;
        omega_foto2 = omega(j,1) ;
        kappa_foto2 = kappa(j,1) ;
        X0_foto2 = X0(j,1) ;
        Y0_foto2 = Y0(j,1) ;
        Z0_foto2 = Z0(j,1) ;
    end
end
%--------------------------------------------------------------------------
%% Calculo das Coordenadas Terreno segundo as Equações de Colineariedade---
    %Conversor de grados para garus
    conv = (pi/200);
    %Matriz de Rotação para a Fotografia 1
    R_omega_foto1 = [1         0                      0             ;
                     0 cos(omega_foto1*conv) -sin(omega_foto1*conv) ;
                     0 sin(omega_foto1*conv) cos(omega_foto1*conv)] ;
           
    R_fi_foto1 = [cos(fi_foto1*conv)   0  sin(fi_foto1*conv) ;
                          0            1          0          ;
                  -sin(fi_foto1*conv)  0  cos(fi_foto1*conv)];
        
    R_kappa_foto1 = [cos(kappa_foto1*conv) -sin(kappa_foto1*conv)  0 ;
                     sin(kappa_foto1*conv)  cos(kappa_foto1*conv)  0 ;
                             0                      0              1];
    
    Rot_foto1 = R_fi_foto1 * R_omega_foto1 * R_kappa_foto1
            
    %Matriz de Rotação A Fotografia 2
    R_omega_foto2 = [1         0                      0             ;
                     0 cos(omega_foto2*conv) -sin(omega_foto2*conv) ;
                     0 sin(omega_foto2*conv)  cos(omega_foto2*conv)];
           
    R_fi_foto2 = [cos(fi_foto2*conv)   0  sin(fi_foto2*conv) ;
                             0         1            0        ;
                  -sin(fi_foto2*conv)  0  cos(fi_foto2*conv)];
        
    R_kappa_foto2 = [cos(kappa_foto2*conv) -sin(kappa_foto2*conv)  0 ;
                     sin(kappa_foto2*conv)  cos(kappa_foto2*conv)  0 ;
                              0                      0             1];
    
    Rot_foto2 = R_fi_foto2 * R_omega_foto2 * R_kappa_foto2
    
    %Escrita do cabeçalho do ficheiro de resultados
    fprintf(resultados,'INTERSECCAO ESPACIAL DIRECTA\n');
    fprintf(resultados,'Realizado por: Rui Jorge Abrunhosa Nunes No32092\n');
    fprintf(resultados,'Mestrado em Engenharia Geografica: Fotogrametria Analitica 2011/2012\n');
    fprintf(resultados,'----------------------------------------------------\n');
    fprintf(resultados,'Fotografias usadas: %d e %d\n', num_foto1, num_foto2);
    fprintf(resultados,'Ponto \t\tX\t\t\t  Y\t\t\t   Z\n');
    fprintf(resultados,'----------------------------------------------------\n');
    
    for k = 1:size(id_foto,1)
        %CÁLCULO DAS APROXIMAÇÕES INICIAIS
        
        %Simplificação das equações
        %Fotografia 1
        Kx_foto1 = (Rot_foto1(1,1)*(x1(k)-x0) + Rot_foto1(1,2)*(y1(k)-y0) - Rot_foto1(1,3)*c) / (Rot_foto1(3,1)*(x1(k)-x0) + Rot_foto1(3,2)*(y1(k)-y0) - Rot_foto1(3,3)*c); 
        Ky_foto1 = (Rot_foto1(2,1)*(x1(k)-x0) + Rot_foto1(2,2)*(y1(k)-y0) - Rot_foto1(2,3)*c) / (Rot_foto1(3,1)*(x1(k)-x0) + Rot_foto1(3,2)*(y1(k)-y0) - Rot_foto1(3,3)*c);
        %Fotografia 2
        Kx_foto2 = (Rot_foto2(1,1)*(x2(k)-x0) + Rot_foto2(1,2)*(y2(k)-y0) - Rot_foto2(1,3)*c) / (Rot_foto2(3,1)*(x2(k)-x0) + Rot_foto2(3,2)*(y2(k)-y0) - Rot_foto2(3,3)*c);
        Ky_foto2 = (Rot_foto2(3,1)*(x2(k)-x0) + Rot_foto2(3,2)*(y2(k)-y0) - Rot_foto2(3,3)*c) / (Rot_foto2(3,1)*(x2(k)-x0) + Rot_foto2(3,2)*(y2(k)-y0) - Rot_foto2(3,3)*c);
     
        Z = 0;
        
        %Cálculo das coordenadas dos Pontos no terreno, da Fotografia 1
        X_foto1(k) = X0_foto1 + (Z-Z0_foto1) * Kx_foto1
        Y_foto1(k) = Y0_foto1 + (Z-Z0_foto1) * Ky_foto1
        
        %Cálculo das coordenadas dos Pontos no terreno, da Fotografia 1
        X_foto2(k) = X0_foto2 + (Z-Z0_foto2) * Kx_foto2;
        Y_foto2(k) = Y0_foto2 + (Z-Z0_foto2) * Ky_foto2;
        
        
        
%--------------------------------------------------------------------------
%% Ajustamento Paramétrico Não Linear
%Número total de observações
n = 4;

%Número de observações necessárias
n0 = 3;

%Graus de Liberdade
df = n-n0;

%Variância à priori
sigma0 = 1;

%Matriz de Variâncias/Covariâncias das observações
Cl = eye(n);

%Matriz dos pesos da observações
Pl = sigma0*inv(Cl);

%Vector de Observações
l = [x1(k);y1(k);x2(k);y2(k)];

%Vector de parâmetros desconhecidos (com aproximações inicias)
X_desc = (X_foto1(k) + X_foto2(k))/2;
Y_desc = (Y_foto1(k) + Y_foto2(k))/2;
Z_desc = (X0_foto2-Z0_foto2 * Kx_foto2+Z0_foto1*Kx_foto1-X0_foto1) / (Kx_foto1 - Kx_foto2);
x = [X_desc; Y_desc; Z_desc];

 
%----------------------------INICIO DO 1ºCICLO------------------------------
precisao = 1e-4;
delta = 1;
interaccao = 0;
while max(abs(delta)) > precisao
    
    %SIMPLIFICAÇÃO DO MODELO MATEMÁTICO
    %Fotografia 1
    Nx_foto1 = Rot_foto1(1,1)*(x(1)-X0_foto1) + Rot_foto1(2,1)*(x(2)-Y0_foto1) + Rot_foto1(3,1)*(x(3)-Z0_foto1);
    Ny_foto1 = Rot_foto1(1,2)*(x(1)-X0_foto1) + Rot_foto1(2,2)*(x(2)-Y0_foto1) + Rot_foto1(3,2)*(x(3)-Z0_foto1);
    D_foto1 =  Rot_foto1(1,3)*(x(1)-X0_foto1) + Rot_foto1(2,3)*(x(2)-Y0_foto1) + Rot_foto1(3,3)*(x(3)-Z0_foto1);
    
    %Fotografia 2
    Nx_foto2 = Rot_foto2(1,1)*(x(1)-X0_foto2) + Rot_foto2(2,1)*(x(2)-Y0_foto2) + Rot_foto2(3,1)*(x(3)-Z0_foto2);
    Ny_foto2 = Rot_foto2(1,2)*(x(1)-X0_foto2) + Rot_foto2(2,2)*(x(2)-Y0_foto2) + Rot_foto2(3,2)*(x(3)-Z0_foto2);
    D_foto2 =  Rot_foto2(1,3)*(x(1)-X0_foto2) + Rot_foto2(2,3)*(x(2)-Y0_foto2) + Rot_foto2(3,3)*(x(3)-Z0_foto2);
    
    %Matriz dos coeficientes das correcções aos parâmetros desconhecidos (Primeira
    %Matriz de Configuração)
    A(1,1) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(1,1) + Nx_foto1*Rot_foto1(1,3));
    A(1,2) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(2,1) + Nx_foto1*Rot_foto1(2,3));
    A(1,3) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(3,1) + Nx_foto1*Rot_foto1(3,3));
    A(2,1) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(1,2) + Ny_foto1*Rot_foto1(1,3));
    A(2,2) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(2,2) + Ny_foto1*Rot_foto1(2,3));
    A(2,3) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(3,2) + Ny_foto1*Rot_foto1(3,3));
    A(3,1) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(1,1) + Nx_foto2*Rot_foto2(1,3));
    A(3,2) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(2,1) + Nx_foto2*Rot_foto2(2,3));
    A(3,3) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(3,1) + Nx_foto2*Rot_foto2(3,3));
    A(4,1) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(1,2) + Ny_foto2*Rot_foto2(1,3));
    A(4,2) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(2,2) + Ny_foto2*Rot_foto2(2,3));
    A(4,3) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(3,2) + Ny_foto2*Rot_foto2(3,3));
   
    %Vector de fecho
    w(1,1) = x0 - c*(Nx_foto1/D_foto1) - l(1);   
    w(2,1) = y0 - c*(Ny_foto1/D_foto1) - l(2);
    w(3,1) = x0 - c*(Nx_foto2/D_foto2) - l(3);
    w(4,1) = y0 - c*(Ny_foto2/D_foto2) - l(4);
    
    %Vector de correcções aos parâmetros desconhecidos
    delta = -inv(A'*Pl*A)*A'*Pl*w
    
    %Vector de parâmetros desconhecidos ajustados
    x = x + delta;
    
    %Vector de Residuos
    v = A*delta + w;
    
    %Vector de observações ajustadas
    l = l + v;
    
    interaccao = interaccao + 1;
end
%----------------------------FIM DO 1º Ciclo-------------------------------
%%
%Matriz de variâncias/covariâncias dos parâmetros desconhecidos
C_delta = sigma0*inv(A'*Pl*A);
C_x = C_delta;

%Matriz de variâncias/covariâncias dos resíduos
C_v = sigma0*(inv(Pl) - A*inv(A'*Pl*A)*A');

%Matriz de variâncias/covariâncias das observações
C_l = Cl-C_v;

%Teste do Factor de Variância (Bilateral)
alpha = 0.95;
teste1 = chi2pdf(alpha/2,df);
teste2 = chi2pdf((1-(alpha/2)),df);
sigma = (v'*Pl*v)/df;

%Matriz de variâncias/covariâncias dos parâmetros desconhecidos
%ajustada
C_delta_ajustado = sigma*inv(A'*Pl*A);
C_x_ajustado = C_delta_ajustado;
    
%Matriz de variâncias/covariâncias dos resíduos ajustada
C_v_ajustado = sigma*(inv(Pl) - A*inv(A'*Pl*A)*A');

%Matriz de variâncias/covariâncias das observações
C_l_ajustado = Cl-C_v_ajustado;
 
%% Escrita dos valores em ficheiro
    fprintf(resultados,'%d\t',id_ponto(k));
    for l = 1:size(x,1) 
        fprintf(resultados,'|%2.3f\t',x(l));
    end
    fprintf(resultados,'\n');
end

fprintf(resultados,'----------------------------------------------------\n');
if (sigma<teste1 && sigma>teste2)
    fprintf(resultados,'\n');
    fprintf(resultados,'Teste do factor de Variancia: Rejeitado\n');
else
    fprintf(resultados,'Teste do factor de Variancia: Aceite\n');
end
%% Fecho de ficheiros
fclose(resultados);