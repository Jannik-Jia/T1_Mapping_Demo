clear all
close all
clc
% Using S1 and S2 reconstructed by step d, T1 (filtered by T2*) is calculated.

%% Required Input data
DicomDatei_Lowflip           = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_lowflip\'];
DicomDatei_flipmap           = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\FlipMaps\'];
DicomDatei_Highflip          = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_highflip_orig\'];
load('ReferenzMatrix.mat')                  % See: b_T1_Hirnmaskierung.m
load('Highflip_Datei_nachdem_T2_Filterung.mat')
                                            % See: d_T2_Stern_Filterung.m
%% Output Informations  
TargetData_T1_Map_mit_T2_Fliterung     = ['T1_nachdem_T2_Filterung.mat']
                                            % The calculated results will be saved to the same folder with this file name.

%% Main Codes                                           
[greLow,iNoFiles,infoFirst]  = loadDicomFolder(DicomDatei_Lowflip,-1,0); 
FA1nom                       = infoFirst.FlipAngle 
TR1                          = infoFirst.RepetitionTime;

[flipmap,iNoFiles,infoFirst] = loadDicomFolder(DicomDatei_flipmap,-1,0);

[greHigh,iNoFiles,infoFirst] = loadDicomFolder(DicomDatei_Highflip,-1,0);
FA2nom                       = infoFirst.FlipAngle;
TR2                          = infoFirst.RepetitionTime;

group_lowflip                = greLow (:,:,1:15).*ReferenzMatrix;
group_highflip               = greHigh(:,:,1:15).*ReferenzMatrix;
group_flipmap                = flipmap(:,:,1:15).*ReferenzMatrix;
     

ratio                        = group_lowflip./group_highflip;
                                            % verhaltnis zwischen S1 und S2;
ratio_new                    = group_lowflip./group_highflip_1_new;   
                                            % verhaltnis zwischen S1 und S2;


%T1 herausfinden
syms t1;
[x,y,z] = size(ratio);
T1_S1   = zeros(x,y,z);

for c=1:1:z                                     % Schichten
    for a=2:1:x                                 % Zeilen
        for b=2:1:y                             % Spalten          
            for t1 = 50:1:3500
                if group_flipmap(a,b,c) ~= 0
                    if ratio(a,b,c)   == Inf    % there are some errors in matrix"ratio", a small number of values are displayed as"Inf".
                        ratio(a,b,c)  = 0;
                    end
                    Eq_lowflip        = ((1-exp(-TR1/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA1nom)*exp(-TR1/t1))))*sind((group_flipmap(a,b,c)/1000)*FA1nom);
                    Eq_Highflip       = ((1-exp(-TR2/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA2nom)*exp(-TR2/t1))))*sind((group_flipmap(a,b,c)/1000)*FA2nom);
                    Equ               = Eq_lowflip/Eq_Highflip;
                    Abweichung        = abs(Equ-ratio(a,b,c));
                    if (t1==50)
                        Abw(a,b,c)    = Abweichung;
                        T1_S1(a,b,c)  = t1;
                    elseif (Abweichung<Abw(a,b,c))
                        Abw(a,b,c)    = Abweichung;
                        T1_S1(a,b,c)  = t1;
                    end
                end
            end
        end
	end
end

%T1 herausfinden
syms t1;
[x,y,z] = size(ratio_new);
T1_S2   = zeros(x,y,z);
for c=1:1:z                                  % Schichten
    for a=2:1:x                              % Zeilen
        for b=2:1:y                          % Spalten          
            for t1 = 50:1:3500
                if group_flipmap(a,b,c) ~= 0
                    if ratio_new(a,b,c) == Inf % there are some errors in matrix "ratio", a small number of values are displayed as"Inf".
                       ratio_new(a,b,c) = 0;
                    end
                    Eq_lowflip          = ((1-exp(-TR1/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA1nom)*exp(-TR1/t1))))*sind((group_flipmap(a,b,c)/1000)*FA1nom);
                    Eq_Highflip         = ((1-exp(-TR2/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA2nom)*exp(-TR2/t1))))*sind((group_flipmap(a,b,c)/1000)*FA2nom);
                    Equ                 = Eq_lowflip/Eq_Highflip;
                    Abweichung          = abs(Equ-ratio_new(a,b,c));
                    if (t1==50)
                       Abw(a,b,c)       = Abweichung;
                       T1_S2(a,b,c)     = t1;
                    elseif (Abweichung<Abw(a,b,c))
                       Abw(a,b,c)       = Abweichung;
                       T1_S2(a,b,c)     = t1;
                    end
                end
            end
        end
	end
end

T1_M = (T1_S1+T1_S2)/2;                                         
                                        % The average value of T1 from both
                                        % steps

%% Rename the calculated data and export it.  The exported file name can be customized at the beginning of the file.

save(TargetData_T1_Map_mit_T2_Fliterung,'T1_M')               
                                        % Save the target T1 Map as mat form in the same folder, which can be imported in the next step.



figure (1)
subplot (351)
imagesc(T1_M(:,:,1), [50 3500]);  colorbar
subplot (352)
imagesc(T1_M(:,:,2), [50 3500]);  colorbar
subplot (353)
imagesc(T1_M(:,:,3), [50 3500]);  colorbar
subplot (354)
imagesc(T1_M(:,:,4), [50 3500]);  colorbar
subplot (355)
imagesc(T1_M(:,:,5), [50 3500]);  colorbar
subplot (356)
imagesc(T1_M(:,:,6), [50 3500]);  colorbar
subplot (357)
imagesc(T1_M(:,:,7), [50 3500]);  colorbar
subplot (358)
imagesc(T1_M(:,:,8), [50 3500]);  colorbar
subplot (359)
imagesc(T1_M(:,:,9), [50 3500]);  colorbar
subplot (3,5,10)
imagesc(T1_M(:,:,10), [50 3500]);  colorbar
subplot (3,5,11)
imagesc(T1_M(:,:,11), [50 3500]);  colorbar
subplot (3,5,12)
imagesc(T1_M(:,:,12), [50 3500]);  colorbar
subplot (3,5,13)
imagesc(T1_M(:,:,13), [50 3500]);  colorbar
subplot (3,5,14)
imagesc(T1_M(:,:,14), [50 3500]);  colorbar
subplot (3,5,15)
imagesc(T1_M(:,:,15), [50 3500]);  colorbar;
sgtitle('MSLOW0015: T1 Relaxationszeit, nachdem T2* Filterung');