clear all
close all
clc
%% The unfiltered T1 relaxation time map filters the noise around the head and retains only the T1 data within the range of the head.

%% Required Input data
DicomDatei_Lowflip         = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_lowflip\'];
DicomDatei_flipmap         = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\FlipMaps\'];
DicomDatei_Highflip        = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_highflip_orig\'];

%% Output Informations  
TargetData_T1_HirnMaskierung = ['T1_Gehirn_Maskierung.mat'];  
                                                                % The calculated results will be saved to the same folder with this file name.
TargetData_ReferenzMatrix  = ['ReferenzMatrix.mat']

%% Main Codes        
% Two different T1 is obtained in two steps, and the average value is calculated, so that a more accurate T1 can be obtained.        
[M_flipmap]                = loadDicomFolder(DicomDatei_flipmap,-1,0);
[M_lowflip]                = loadDicomFolder(DicomDatei_Lowflip,-1,0); 
[M_highflip]               = loadDicomFolder(DicomDatei_Highflip,-1,0);

highflip_info_S1           = dicominfo([DicomDatei_Highflip,'File0001.ima']) % FlipAngle: 14; RepetitionTime: 200; EchoTime: 3.9500
lowflip_info_S1            = dicominfo([DicomDatei_Lowflip,'File0001.ima'])  % FlipAngle: 50; RepetitionTime: 300; EchoTime: 3.9500
highflip_info_S2           = dicominfo([DicomDatei_Highflip,'File0016.ima']) % FlipAngle: 14; RepetitionTime: 200; EchoTime: 11
lowflip_info_S2            = dicominfo([DicomDatei_Lowflip,'File0061.ima'])  % FlipAngle: 50; RepetitionTime: 300; EchoTime: 11

group_flipmap              = M_flipmap(:,:, 1:15);              % Extract the first 15 dicom files, that is, a complete brain measurement
group_lowflip_S1           = M_lowflip(:,:, 1:15);
group_lowflip_S2           = M_lowflip(:,:,61:75);
group_highflip_S1          = M_highflip(:,:, 1:15);
group_highflip_S2          = M_highflip(:,:,16:30);

%% DICOM_Info auslesen
% highflip info:
TR2_S1                     = highflip_info_S1.RepetitionTime;   % Schritt 1
TR2_S2                     = highflip_info_S2.RepetitionTime;   % Schritt 2
TE_Highflip_S1             = highflip_info_S1.EchoTime 
TE_Highflip_S2             = highflip_info_S2.EchoTime 
FA2nom_S1                  = highflip_info_S1.FlipAngle;        % Schritt 1
FA2nom_S2                  = highflip_info_S2.FlipAngle;        % Schritt 2
%lowflip info:      
TR1_S1                     = lowflip_info_S1.RepetitionTime;    % Schritt 1
TR1_S2                     = lowflip_info_S2.RepetitionTime;    % Schritt 2
TE_Lowflip_S1              = lowflip_info_S1.EchoTime 
TE_Lowflip_S2              = lowflip_info_S2.EchoTime 
FA1nom_S1                  = lowflip_info_S1.FlipAngle;         % Schritt 1
FA1nom_S2                  = lowflip_info_S2.FlipAngle;         % Schritt 2
            
%verhaltnis zwischen S1 und S2;       
%Schritt 1
ratio_S1                   = group_lowflip_S1./group_highflip_S1;
%Schritt 2               
ratio_S2                   = group_lowflip_S2./group_highflip_S2;
       
%T1 herausfinden       
[x,y,z]                    = size(ratio_S2);                    
T1_S1                      = zeros(x,y,z);                      % T1 calculated in the first step 
T1_S2                      = zeros(x,y,z);                      % T1 calculated in the second step  
T1_M                       = zeros(x,y,z);                      % The average value of T1 from both steps
Abw1                       = zeros(x,y,z);                      % Abweichung  = abs((Eq_lowflip/Eq_Highflip)-ratio)
Abw2                       = zeros(x,y,z);
ReferenzMatrix             = zeros(x,y,z);
noise_lowfilp_S1           = group_lowflip_S1 (7:11,7:11,1);    % Random background noise: 3*3 Matrix
noise_highflip_S1          = group_highflip_S1(10:14,4:8,1);    % Random background noise: 3*3 Matrix
NoiseLevel_lowflip_S1      = mean(noise_lowfilp_S1(:));
NoiseLevel_highflip_S1     = mean(noise_highflip_S1(:));
noise_Verstaerkung_S1      = 5;
denoise_parameter_low_S1   = NoiseLevel_lowflip_S1  * noise_Verstaerkung_S1;
denoise_parameter_high_S1  = NoiseLevel_highflip_S1 * noise_Verstaerkung_S1;

noise_lowfilp_S2           = group_lowflip_S2 (7:11,7:11,1);    % Random background noise: 3*3 Matrix
noise_highflip_S2          = group_highflip_S2(10:14,4:8,1);    % Random background noise: 3*3 Matrix
NoiseLevel_lowflip_S2      = mean(noise_lowfilp_S2(:));
NoiseLevel_highflip_S2     = mean(noise_highflip_S2(:));
noise_Verstaerkung_S2      = 5;
denoise_parameter_low_S2   = NoiseLevel_lowflip_S2   * noise_Verstaerkung_S2;
denoise_parameter_high_S2  = NoiseLevel_highflip_S2  * noise_Verstaerkung_S2;



for c=1:1:z                                                     % Schichten
% for c=1    
    for a=2:1:x                                                 % Zeilen
        for b=2:1:y                                             % Spalten          
            for t1 = 50:1:3500
                if group_flipmap(a,b,c) ~= 0
                                                                % T1 Mapping: Schritt1
                   if group_lowflip_S1(a,b,c) > denoise_parameter_low_S1 && group_highflip_S1(a,b,c) > denoise_parameter_high_S1
                        if  ratio_S1(a,b,c) ==Inf               % there are some errors in matrix"ratio", a small number of values are displayed as"Inf".
                            ratio_S1(a,b,c) = 0;
                        end
                            Eq_lowflip_S1   = ((1-exp(-TR1_S1/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA1nom_S1)*exp(-TR1_S1/t1))))*sind((group_flipmap(a,b,c)/1000)*FA1nom_S1);
                            Eq_Highflip_S1  = ((1-exp(-TR2_S1/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA2nom_S1)*exp(-TR2_S1/t1))))*sind((group_flipmap(a,b,c)/1000)*FA2nom_S1);
                            Equ_S1          = Eq_lowflip_S1/Eq_Highflip_S1;
                            Abweichung_S1   = abs(Equ_S1-ratio_S1(a,b,c));
                        if (t1==50)
                            Abw1(a,b,c)     = Abweichung_S1;
                            T1_S1(a,b,c)    = t1;
                        elseif (Abweichung_S1<Abw1(a,b,c))
                            Abw1(a,b,c)     = Abweichung_S1;
                            T1_S1(a,b,c)    = t1;
                        end
                    end
                    %Schtitt1 end
                    
                    %T1 Mapping: Schritt2
                    if group_lowflip_S2(a,b,c) > denoise_parameter_low_S2 && group_highflip_S2(a,b,c) > denoise_parameter_high_S2
                        if  ratio_S2(a,b,c) == Inf              % there are some errors in matrix"ratio", a small number of values are displayed as"Inf".
                            ratio_S2(a,b,c) =  0;
                        end
                            Eq_lowflip_S2   = ((1-exp(-TR1_S2/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA1nom_S2)*exp(-TR1_S2/t1))))*sind((group_flipmap(a,b,c)/1000)*FA1nom_S2);
                            Eq_Highflip_S2  = ((1-exp(-TR2_S2/t1))/(1-(cosd((group_flipmap(a,b,c)/1000)*FA2nom_S2)*exp(-TR2_S2/t1))))*sind((group_flipmap(a,b,c)/1000)*FA2nom_S2);
                            Equ_S2          = Eq_lowflip_S2/Eq_Highflip_S2;
                            Abweichung_S2   = abs(Equ_S2-ratio_S2(a,b,c));
                        if (t1==50)
                            Abw2(a,b,c)     = Abweichung_S2;
                            T1_S2(a,b,c)    = t1;
                        elseif (Abweichung_S2<Abw2(a,b,c))
                            Abw2(a,b,c)     = Abweichung_S2;
                            T1_S2(a,b,c)    = t1;
                        end
                    end
                    %Schtitt2 end
                end
            end
        end
	end
end

%% ReferenzMatrix_Hirnmaskierung
for c=1:1:z        %Schichten
% for c=1    
    for a=2:1:x     %Zeilen
        for b=2:1:y %Spalten          
            for t1 = 50:1:3500
                if group_flipmap(a,b,c) ~= 0
                   if group_lowflip_S1(a,b,c) > denoise_parameter_low_S1 && group_highflip_S1(a,b,c) > denoise_parameter_high_S1
                        if  ratio_S1(a,b,c)==Inf
                            ratio_S1(a,b,c) =0;
                        end
                       ReferenzMatrix(a,b,c) = 1;
                   else
                       ReferenzMatrix(a,b,c) = 0;
                   end
                end
            end
        end
	end
end


T1_S1;
T1_S2;
T1_M = (T1_S1+T1_S2)/2;                                         % The average value of T1 from both steps


%% Rename the calculated data and export it.  The exported file name can be customized at the beginning of the file.

save(TargetData_T1_HirnMaskierung,'T1_M')               % Save the target T1 Map as mat form in the same folder, which can be imported in the next step.
save(TargetData_ReferenzMatrix,'ReferenzMatrix') 

%% The results are shown below.
figure (1)
subplot (351)
imagesc(T1_S1(:,:,1), [50 3500]);  colorbar
subplot (352)
imagesc(T1_S1(:,:,2), [50 3500]);  colorbar
subplot (353)
imagesc(T1_S1(:,:,3), [50 3500]);  colorbar
subplot (354)
imagesc(T1_S1(:,:,4), [50 3500]);  colorbar
subplot (355)
imagesc(T1_S1(:,:,5), [50 3500]);  colorbar
subplot (356)
imagesc(T1_S1(:,:,6), [50 3500]);  colorbar
subplot (357)
imagesc(T1_S1(:,:,7), [50 3500]);  colorbar
subplot (358)
imagesc(T1_S1(:,:,8), [50 3500]);  colorbar
subplot (359)
imagesc(T1_S1(:,:,9), [50 3500]);  colorbar
subplot (3,5,10)
imagesc(T1_S1(:,:,10), [50 3500]);  colorbar
subplot (3,5,11)
imagesc(T1_S1(:,:,11), [50 3500]);  colorbar
subplot (3,5,12)
imagesc(T1_S1(:,:,12), [50 3500]);  colorbar
subplot (3,5,13)
imagesc(T1_S1(:,:,13), [50 3500]);  colorbar
subplot (3,5,14)
imagesc(T1_S1(:,:,14), [50 3500]);  colorbar
subplot (3,5,15)
imagesc(T1_S1(:,:,15), [50 3500]);  colorbar;
sgtitle('MSLOW0015: T1 Relaxationszeit, nachdem 1.Optimierung  und Hirnmaskierung, TE = 3,95ms');

figure (2)
subplot (351)
imagesc(T1_S2(:,:,1), [50 3500]);  colorbar
subplot (352)
imagesc(T1_S2(:,:,2), [50 3500]);  colorbar
subplot (353)
imagesc(T1_S2(:,:,3), [50 3500]);  colorbar
subplot (354)
imagesc(T1_S2(:,:,4), [50 3500]);  colorbar
subplot (355)
imagesc(T1_S2(:,:,5), [50 3500]);  colorbar
subplot (356)
imagesc(T1_S2(:,:,6), [50 3500]);  colorbar
subplot (357)
imagesc(T1_S2(:,:,7), [50 3500]);  colorbar
subplot (358)
imagesc(T1_S2(:,:,8), [50 3500]);  colorbar
subplot (359)
imagesc(T1_S2(:,:,9), [50 3500]);  colorbar
subplot (3,5,10)
imagesc(T1_S2(:,:,10), [50 3500]);  colorbar
subplot (3,5,11)
imagesc(T1_S2(:,:,11), [50 3500]);  colorbar
subplot (3,5,12)
imagesc(T1_S2(:,:,12), [50 3500]);  colorbar
subplot (3,5,13)
imagesc(T1_S2(:,:,13), [50 3500]);  colorbar
subplot (3,5,14)
imagesc(T1_S2(:,:,14), [50 3500]);  colorbar
subplot (3,5,15)
imagesc(T1_S2(:,:,15), [50 3500]);  colorbar;
sgtitle('MSLOW0015: T1 Relaxationszeit, nachdem 1.Optimierung  und Hirnmaskierung, TE = 11ms');

figure (3)
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
sgtitle('MSLOW0015: T1 Relaxationszeit, nachdem Optimierung  und Hirnmaskierung');


