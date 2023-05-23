% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%  Estruturas de Concreto I - 2023.1 - Prof. Augusto Albuquerque
%  Departamento de Engenharia Estrutural e Construção Civil - DEECC
%  Universidade Federal do Ceará - UFC
% ------------------------------------------------------------------------
%
%  LAJE
%  Programa para cálculo de momento em Lajes Maciças segundo as tabelas de
%  Marcus
%
% ------------------------------------------------------------------------
%
%  por:
%  Luiz Carlos Matias Teixeira
%
% ------------------------------------------------------------------------
% to-do list:
% - corrigir armadura mínima
%
% ---------------------------------------------------------------------
clc
%-------------------------------------------------------------------------
%Tabelas de Marcus
tipo_marcus = ["tipo1" "tipo2" "tipo3" "tipo4" "tipo5" "tipo6"];

global fsig;
fsig = fopen(strcat('Resultados.txt'),'w');

%faz a leitura da planilha
filename = 'marcus.xlsx';
arq_momento = 'momento.xlsx';
momento = readmatrix(arq_momento);
nm = size(momento(:,1),1);
M = zeros(nm,4);
T = zeros(nm,1);

for i=1:nm
    % ---------------------------
    % Calculos dos esforços
    tipo = momento(i+4*nm);
    h = momento(i+3*nm);
    marcus = readmatrix(filename,"Sheet",tipo_marcus(tipo));
    nx = size(marcus(:,1),1);
    lxx = momento(i+1*nm);
    lyy = momento(i+2*nm);
    lambm = momento(i+2*nm)/momento(i+1*nm);
    lambm = round(lambm*100)/100;
    if lambm > 2
        M(i,1) = 0.00;
        M(i,2) = 0.00;
        M(i,3) = 0.00;
        M(i,4) = 0.00;
    else
    %localiza a posição da linha de lambda
    indices = find(marcus(:,1) == round(lambm*100)/100);
    
    %localiza a posição da linha de kx, mx, ny, my, ny
    ikx = (indices + 1*nx);
    imx = (indices + 2*nx);
    inx = (indices + 3*nx);
    imy = (indices + 4*nx);
    iny = (indices + 5*nx);
    
    %e mx, ny, my, ny
    kxm = marcus(ikx);
    mxm = marcus(imx);
    mym = marcus(imy);
    
    if marcus(inx)==1e24
        nxm = inf;    
    else
        nxm = marcus(inx);
    end
 
    if marcus(iny)==1e24
        nym = inf;   
    else
        nym = marcus(iny);
    end


    M(i,1) = 1.4*(momento(i+10*nm))*(momento(i+1*nm))^2/marcus(imx);
    M(i,2) = 1.4*(momento(i+10*nm))*(momento(i+1*nm))^2/marcus(imy);
    M(i,3) = 1.4*(momento(i+10*nm))*(momento(i+1*nm))^2/marcus(inx);
    M(i,4) = 1.4*(momento(i+10*nm))*(momento(i+1*nm))^2/marcus(iny);
    
    fprintf(fsig, '--------------------------------------------------------\n');
    fprintf(fsig, 'CÁLCULO DOS ESFORÇOS - LAJE L%d',i);
    fprintf(fsig, '\n--------------------------------------------------------\n');
    fprintf(fsig, 'lx = %.2f m',lxx);
    fprintf(fsig, '\nly = %.2f m',lyy);
    fprintf(fsig, '\nλ = %f',lyy/lxx);
    fprintf(fsig, '\nλ = %.2f (Valor utilizado nas tabelas de Marcus)',round(lyy/lxx*100)/100);
    fprintf(fsig, '\n\nCondição de Apoio: Tipo %d',tipo);
    fprintf(fsig, '\nmx = %.4f',mxm);
    fprintf(fsig, '\nmy = %.4f',mym);
    fprintf(fsig, '\nnx = %.4f',nxm);
    fprintf(fsig, '\nny = %.4f',nym);
    %Cargas
    fprintf(fsig, '\n\n-----Cargas-----\nPeso Próprio + Pavimentação = %.2f kN/m²',momento(i+7*nm));
    fprintf(fsig, '\nCarga de Alvenaria = %.2f kN/m²',momento(i+8*nm));
    fprintf(fsig, '\nCarga Permanente (g) = %.2f kN/m²',momento(i+8*nm)+momento(i+7*nm));
    fprintf(fsig, '\nCarga Acidental (q)= %.2f kN/m²',momento(i+9*nm));
    fprintf(fsig, '\nCarga total (p) = g + q = %.2f kN/m²',momento(i+10*nm));
    %Momentos
    fprintf(fsig, '\n\n-----Momentos-----\nMx = p*lx^2/mx\nMy = p*lx^2/my\nXx = p*lx^2/nx\nXy = p*lx^2/ny');
    fprintf(fsig, '\n\n-----Momentos Positivos-----\nMx = %.4f kNm/m',M(i,1)/1.4);
    fprintf(fsig, '\nMy = %.4f kNm/m',M(i,2)/1.4);
    fprintf(fsig, '\n\n-----Momentos Negativos-----\nXx = %.4f kNm/m',M(i,3)/1.4);
    fprintf(fsig, '\nXy = %.4f kNm/m\n',M(i,4)/1.4);    
    

    %dimensionamento
    num_ferrox = 1:2:2*nm;
    num_ferroy = 2:2:2*nm+1;
    %h = 14.93;
    cobrimento = 3;
    d = h-cobrimento-1;
    fck = 30;
    fyk = 500;
    %bitx = 8;
   % bity = 8;
    

    msdx = M(i,1)*100;
    msdxt = msdx;
    msdy = M(i,2)*100; 
    msdyt = msdy;
    
    Mdmin = 0.8*(100/6)*(h/100)^2*1.3*0.3*(fck)^(2/3)*100;
    %fprintf('\nMdmin = %.2f kN.cm/m\n', Mdmin);
    if Mdmin > msdx
        msdx = Mdmin;
        %fprintf('\nMx(%d) = %.2f kN.cm/m\n',i, msdx);
    end
    if Mdmin > msdy
        msdy = Mdmin;
        %fprintf('\nMy(%d) = %.2f kN.cm/m\n',i, msdy);
    end


    %fprintf('\nMx(%d) = %.2f kN.cm/m\n',i, msdx);
    %fprintf('\nMy(%d) = %.2f kN.cm/m\n',i, msdy);
    %faz a leitura da planilha
    bitola = 'bitola.xlsx';
    bit50 = readmatrix(bitola,"Sheet","ca-50");
    nb = size(bit50(:,1),1);
    
   
    xx = (0.68*d-sqrt((0.68*d)^2-4*0.272*(msdx/(100*fck*0.1/1.4))))/0.544;
    %fprintf('\nxx = %.2f cm\n', xx);
    
    xy = (0.68*d-sqrt((0.68*d)^2-4*0.272*(msdy/(100*fck*0.1/1.4))))/0.544;
    %fprintf('\nxy = %.2f cm\n', xy);
    
    zx = d - 0.4*xx;
    Asx = msdx/(zx*(fyk*0.1/1.15));
    %fprintf('\nAsx Antes(%d) = %.2f cm²\n',i, Asx);
    Asxt = Asx;
    
    zy = d - 0.4*xy;
    Asy = msdy/(zy*(fyk*0.1/1.15));
    Asyt = Asy;
    %fprintf('\nAsy Antes(%d) = %.2f cm²\n',i, Asy);  
    Asmin = 0.0015*100*h;
    %fprintf('\nAsmin(%d) = %.2f cm²\n',i, Asmin);

    if Asmin > Asx
        Asx = Asmin;
        %fprintf('\nAsx(%d) = %.2f cm²\n',i, Asx);
    end
    %fprintf('\nAsx depois(%d) = %.2f cm²\n',i, Asx);
    if Asmin > Asy
        Asy = Asmin;
        %fprintf('\nAsy(%d) = %.2f cm²\n',i, Asy);
    end
    %fprintf('\nAsy depois(%d) = %.2f cm²\n',i, Asy);
    
    matS = zeros(nb,4);
    for ib=1:nb
        matS(ib,1)= bit50(ib,2)*100/Asx;
        matS(ib,2)= bit50(ib,2)*100/Asy;
        matS(ib,3)= bit50(ib,1);
        matS(ib,4)= bit50(ib,2);
    end 
    
    fprintf('\nib = %d e i = %d\n',ib, i);
    matS;
    sbx = max(matS(matS(:,2) <= 23,1));
    sby = max(matS(matS(:,2) <= 23,2));
    
    
    %localiza a posição da linha da bitola na tabela de área
    indx = find(matS(:,1) == sbx);
    indy = find(matS(:,2) == sby);

    Ax = matS(indx,4);
    Ay = matS(indy,4);
    
    sx = Ax*100/Asx;
    SX = floor(sx);
    
    sy = Ay*100/Asy;
    SY = floor(sy);
    
    qtdx = (lyy*100-14)/SX;
    QTDX = ceil(qtdx);
    lxf = ceil(lxx*100-2*cobrimento);
    
    qtdy = (lxx*100-14)/sy;
    QTDY = ceil(qtdy);
    lyf = ceil(lyy*100-2*cobrimento);

    fprintf(fsig, '\n--------------------------------------------------------\n');
    %fprintf(fsig, 'DIMENSIONAMENTO DO FERRO - LAJE L%d',i);
    fprintf(fsig,'DIMENSIONAMENTO DO FERRO\nLAJE L%d - (%.2f m x %.2f m)', i,lxx,lyy);
    fprintf(fsig, '\n--------------------------------------------------------');   
    fprintf(fsig,'\nAltura Útil (d) = %.2f cm', d);
    fprintf(fsig,'\ncobrimento (c) = %.2f cm', cobrimento);
    fprintf(fsig,'\nConcreto C%d (fck = %.2f MPa)', fck,fck);
    fprintf(fsig,'\nAço CA-%d (fyk = %.2f MPa)', fyk/10,fyk);
    
    fprintf(fsig, '\n\n--------------------------------------------------------');    
    fprintf(fsig,'\nCálculo da armadura em X');
    fprintf(fsig, '\n--------------------------------------------------------');
    fprintf(fsig,'\nCOLOCAR A FÓRMULA DO X(linha neutra)\n');
    fprintf(fsig,'\nlinha neutra x = %.2f cm', xx);
    fprintf(fsig,'\nMsdmin = %.4f kN.cm/m', Mdmin);
    fprintf(fsig,'\nMsd = 1.4 ⋅ Mx = %.4f kN.cm/m', msdx);
    fprintf(fsig,'\nBraço de alavanca (z) = d - 0.4x = %.2f cm', zx);
    fprintf(fsig,'\nÁrea de aço (Asx) = %.2f cm²', Asxt);
    fprintf(fsig,'\nFerro utilizado: Ø%d (área = %.2f cm²)',bitx, Ax);
    fprintf(fsig,'\nEspaçamento (Espx) = Øs ⋅ 100/Asx = %.2f ⋅ 100/%.2f = %.2f cm',Ax, Asx,sx);
    fprintf(fsig,'\nEspx = %.2f cm (projeto)', SX);
    fprintf(fsig, '\n--------------------------------------------------------');
    fprintf(fsig,'\nFERRO EM X: %d N%d Ø%.1f C/%d C = %d\n', QTDX,num_ferrox(1,i),bitx,SX,lxf);
    fprintf(fsig, '--------------------------------------------------------\n');
 
    fprintf(fsig, '\n\n--------------------------------------------------------');    
    fprintf(fsig,'\nCálculo da armadura em Y');
    fprintf(fsig, '\n--------------------------------------------------------');  
    fprintf(fsig,'\nCOLOCAR A FÓRMULA DO X(linha neutra)\n');
    fprintf(fsig,'\nlinha neutra x = %.2f cm', xy);
    fprintf(fsig,'\nMsdmin = %.4f kN.cm/m', Mdmin);
    fprintf(fsig,'\nMsd = 1.4 ⋅ My = %.4f kN.cm/m', msdy);
    fprintf(fsig,'\nBraço de alavanca (z) = d - 0.4x = %.2f cm', zy); 
    fprintf(fsig,'\nÁrea de aço (Asy) = %.2f cm²', Asyt); 
    fprintf(fsig,'\nFerro utilizado: Ø%d (área = %.2f cm²)',bity, Ay);
    fprintf(fsig,'\nEspaçamento (Espy) = Øs ⋅ 100/Asy = %.2f ⋅ 100/%.2f = %.2f cm',Ay, Asy,sy);
    fprintf(fsig,'\nEspy = %.2f cm (projeto)', SY); 
    fprintf(fsig,'\n--------------------------------------------------------');
    fprintf(fsig,'\nFERRO EM Y: %d N%d Ø%.1f C/%d C = %d\n', QTDY,num_ferroy(1,i),bity,SY,lyf);
    fprintf(fsig,'--------------------------------------------------------\n');
    end


    % --------------------------------------------------------------------
    
 





end

fprintf('\n--------------------------------------------------------\n')
fprintf(' OK');
fprintf('\n-------------------------------------------------------- \n');

M
