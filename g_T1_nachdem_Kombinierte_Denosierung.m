clear all
close all
clc

%% Required Input data
DicomDatei_Lowflip           = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_lowflip\'];
DicomDatei_flipmap           = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\FlipMaps\'];
DicomDatei_Highflip          = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_highflip_orig\'];
load('ReferenzMatrix.mat')                  % See: b_T1_Hirnmaskierung.m
load('Highflip_Datei_nachdem_Raumlichen_Denoising.mat')
                                            % See: e_T1_Deniosing.m
%% Output Informations  
TargetData_T1_Map_mit_Kombinierte_Denosierung     = ['T1_nachdem_Kombinierte_Denosierung.mat']
                                            % The calculated results will be saved to the same folder with this file name.

%% Main Codes                                           
[greLow,iNoFiles,infoFirst]        = loadDicomFolder(DicomDatei_Lowflip,-1,0); 
FA1nom                             = infoFirst.FlipAngle 
TR1                                = infoFirst.RepetitionTime;
      
[flipmap,iNoFiles,infoFirst]       = loadDicomFolder(DicomDatei_flipmap,-1,0);
      
[greHigh,iNoFiles,infoFirst]       = loadDicomFolder(DicomDatei_Highflip,-1,0);
FA2nom                             = infoFirst.FlipAngle;
TR2                                = infoFirst.RepetitionTime;

group_lowflip                      = greLow (:,:,1:15).*ReferenzMatrix;
group_highflip_1_new_raum_8nachbar = group_highflip_1_new_raum_8nachbar.*ReferenzMatrix;
group_highflip_2_new_raum_8nachbar = group_highflip_2_new_raum_8nachbar.*ReferenzMatrix;
group_flipmap                      = flipmap(:,:,1:15).*ReferenzMatrix;



ratio_S1                           = group_lowflip./group_highflip_1_new_raum_8nachbar;
ratio_S2                           = group_lowflip./group_highflip_2_new_raum_8nachbar;



%T1 herausfinden
syms t1;
[x,y,z]=size(ratio_S1);
T1_S1=zeros(x,y,z);

for c=1:1:z         %Schichten
    for a=2:1:x     %Zeilen
        for b=2:1:y %Spalten          
            for t1 = 50:1:3500
                if group_flipmap(a,b,c)  ~= 0
                    if ratio_S1(a,b,c)   ==Inf          % there are some errors in matrix"ratio_S1", a small number of values are displayed as"Inf".
                        ratio_S1(a,b,c)  =0;
                    end 
                    Eq_lowflip           = ((1-exp(-TR1/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA1nom)*exp(-TR1/t1))))*sind((group_flipmap(a,b,c)/1000)*FA1nom);
                    Eq_Highflip          = ((1-exp(-TR2/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA2nom)*exp(-TR2/t1))))*sind((group_flipmap(a,b,c)/1000)*FA2nom);
                    Equ                  = Eq_lowflip/Eq_Highflip;
                    Abweichung           = abs(Equ-ratio_S1(a,b,c));
                    if (t1==50) 
                        Abw(a,b,c)       = Abweichung;
                        T1_S1(a,b,c)     = t1;
                    elseif (Abweichung<Abw(a,b,c))
                            Abw(a,b,c)   = Abweichung;
                            T1_S1(a,b,c) =t1;
                    end
                end
            end
        end
	end
end

%T1 herausfinden
syms t1;
[x,y,z]=size(ratio_S2);
T1_S2=zeros(x,y,z);
for c=1:1:z         %Schichten
    for a=2:1:x     %Zeilen
        for b=2:1:y %Spalten          
            for t1 = 50:1:3500
                if group_flipmap(a,b,c)  ~= 0
                    if ratio_S2(a,b,c)   == Inf     %there are some errors in matrix"ratio_S2", a small number of values are displayed as"Inf".
                        ratio_S2(a,b,c)  = 0;
                    end
                    Eq_lowflip           = ((1-exp(-TR1/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA1nom)*exp(-TR1/t1))))*sind((group_flipmap(a,b,c)/1000)*FA1nom);
                    Eq_Highflip          = ((1-exp(-TR2/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA2nom)*exp(-TR2/t1))))*sind((group_flipmap(a,b,c)/1000)*FA2nom);
                    Equ                  = Eq_lowflip/Eq_Highflip;
                    Abweichung           = abs(Equ-ratio_S2(a,b,c));
                    if (t1==50)
                        Abw(a,b,c)       = Abweichung;
                        T1_S2(a,b,c)     = t1;
                    elseif (Abweichung < Abw(a,b,c))
                            Abw(a,b,c)   = Abweichung;
                            T1_S2(a,b,c) = t1;
                    end
                end
            end
        end
	end
end


%% Rename the calculated data and export it.  The exported file name can be customized at the beginning of the file.
T1_After_Denoising_Raum_8nachbar         = (T1_S1+T1_S2)/2;

save(TargetData_T1_Map_mit_Kombinierte_Denosierung,'T1_After_Denoising_Raum_8nachbar')               
                                                            % Save the target T1 Map as mat form in the same folder, which can be imported in the next step.

figure (1)
subplot (351)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,1), [50 3500]);  colorbar
subplot (352)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,2), [50 3500]);  colorbar
subplot (353)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,3), [50 3500]);  colorbar
subplot (354)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,4), [50 3500]);  colorbar
subplot (355)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,5), [50 3500]);  colorbar
subplot (356)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,6), [50 3500]);  colorbar
subplot (357)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,7), [50 3500]);  colorbar
subplot (358)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,8), [50 3500]);  colorbar
subplot (359)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,9), [50 3500]);  colorbar
subplot (3,5,10)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,10), [50 3500]);  colorbar
subplot (3,5,11)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,11), [50 3500]);  colorbar
subplot (3,5,12)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,12), [50 3500]);  colorbar
subplot (3,5,13)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,13), [50 3500]);  colorbar
subplot (3,5,14)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,14), [50 3500]);  colorbar
subplot (3,5,15)
imagesc(T1_After_Denoising_Raum_8nachbar(:,:,15), [50 3500]);  colorbar;
sgtitle('MSLOW0015: T1 Map nachdem kombinierter Rauschunterdrückung, 8 Nachbarschaft');
