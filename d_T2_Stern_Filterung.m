clear all
close all
clc
% Calculate the S0 filtered by T2* and reconstruct the new S1 and S2. 
% ** T1 is not calculated in this step!**

%% Required Input data
DicomDatei_Lowflip          = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_lowflip\'];
DicomDatei_flipmap          = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\FlipMaps\'];
DicomDatei_Highflip         = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_highflip_orig\'];
load('ReferenzMatrix.mat')                  % See: b_T1_Hirnmaskierung.m
load('T2Sternmap.mat')                      % See: c_T2Stern_Map_Creating.m

%% Output Informations  
TargetData_new_S1S2_Info    = ['Highflip_Datei_nachdem_T2_Filterung.mat']
                                            % The calculated results will be saved to the same folder with this file name.

%% Main Codes 
% DICOM_Matrix einladen
[M_highflip]                = loadDicomFolder(DicomDatei_Highflip,-1,0);

% Segment the data of the head and remove the noise around the head
group_highflip_1            = M_highflip(:,:, 1:15).*ReferenzMatrix;
                                            % Segment the data of the head and remove the noise around the head
group_highflip_2            = M_highflip(:,:,16:30).*ReferenzMatrix;
                                            % Segment the data of the head and remove the noise around the head
% DICOM_Info auslesen
highflip_info_S1            = dicominfo([DicomDatei_Highflip,'File0001.ima']) 
                                            % FlipAngle: 14; RepetitionTime: 200; EchoTime: 3.9500
highflip_info_S2            = dicominfo([DicomDatei_Highflip,'File0016.ima']) 
                                            % FlipAngle: 14; RepetitionTime: 200; EchoTime: 11
TE_Highflip_S1              = highflip_info_S1.EchoTime
TE_Highflip_S2              = highflip_info_S2.EchoTime
        

[x,y,z]                     = size(group_highflip_1);
group_highflip_1_new        = zeros(x,y,z);
group_highflip_2_new        = zeros(x,y,z);
S0_Map                      = zeros(x,y,z);
varianz_Map                 = zeros(x,y,z);

for c=1:1:z                                 % Schichten
    for a=2:1:x                             % Zeilen
        for b=2:1:y                         % Spalten
            if group_highflip_1(a,b,c) ~=0
                m1                          = exp((-TE_Highflip_S1)/T2_Stern_Map(a,b,c)); 
                                                % Variables used to simplify the formula
                n1                          = exp((-TE_Highflip_S2)/T2_Stern_Map(a,b,c));
                                                % Variables used to simplify the formula
                S1                          = group_highflip_1(a,b,c);
                S2                          = group_highflip_2(a,b,c);
                S0                          = (m1*S1+n1*S2)/(m1^2+n1^2);
                S0_Map(a,b,c)               = S0;
                group_highflip_1_new(a,b,c) = round(S0*exp((-TE_Highflip_S1)/T2_Stern_Map(a,b,c)));%S1 new;
                group_highflip_2_new(a,b,c) = round(S0*exp((-TE_Highflip_S2)/T2_Stern_Map(a,b,c)));%S2 new;
                varianz_Map(a,b,c)          = (group_highflip_1_new(a,b,c)-(S0*exp((-TE_Highflip_S1)/T2_Stern_Map(a,b,c))))^2 + (group_highflip_2_new(a,b,c)-(S0*exp((-TE_Highflip_S2)/T2_Stern_Map(a,b,c))))^2;
            end
        end
    end
end

%% Rename the calculated data and export it.  The exported file name can be customized at the beginning of the file.
save(TargetData_new_S1S2_Info,'group_highflip_1_new','group_highflip_2_new')              
                                            % Save the new created S1 and S2 as mat form in the same folder, which can be imported in the next step.


figure (1)
subplot (351)
imagesc(group_highflip_1_new(:,:,1), [0 800]); colorbar
subplot (352)
imagesc(group_highflip_1_new(:,:,2), [0 800]); colorbar
subplot (353)
imagesc(group_highflip_1_new(:,:,3), [0 800]); colorbar
subplot (354)
imagesc(group_highflip_1_new(:,:,4), [0 800]); colorbar
subplot (355)
imagesc(group_highflip_1_new(:,:,5), [0 800]); colorbar
subplot (356)
imagesc(group_highflip_1_new(:,:,6), [0 800]); colorbar
subplot (357)
imagesc(group_highflip_1_new(:,:,7), [0 800]); colorbar
subplot (358)
imagesc(group_highflip_1_new(:,:,8), [0 800]); colorbar
subplot (359)
imagesc(group_highflip_1_new(:,:,9), [0 800]); colorbar
subplot (3,5,10)
imagesc(group_highflip_1_new(:,:,10), [0 800]); colorbar
subplot (3,5,11)
imagesc(group_highflip_1_new(:,:,11), [0 800]); colorbar
subplot (3,5,12)
imagesc(group_highflip_1_new(:,:,12), [0 800]); colorbar
subplot (3,5,13)
imagesc(group_highflip_1_new(:,:,13), [0 800]); colorbar
subplot (3,5,14)
imagesc(group_highflip_1_new(:,:,14), [0 800]); colorbar
subplot (3,5,15)
imagesc(group_highflip_1_new(:,:,15), [0 800]); colorbar;
sgtitle('MSLOW0015:  Gre Highflip gefiltert und rekonstruiert nach T2*, TE = 3.95ms');


figure (2)
subplot (351)
imagesc(group_highflip_2_new(:,:,1), [0 600]); colorbar
subplot (352)
imagesc(group_highflip_2_new(:,:,2), [0 600]); colorbar
subplot (353)
imagesc(group_highflip_2_new(:,:,3), [0 600]); colorbar
subplot (354)
imagesc(group_highflip_2_new(:,:,4), [0 600]); colorbar
subplot (355)
imagesc(group_highflip_2_new(:,:,5), [0 600]); colorbar
subplot (356)
imagesc(group_highflip_2_new(:,:,6), [0 600]); colorbar
subplot (357)
imagesc(group_highflip_2_new(:,:,7), [0 600]); colorbar
subplot (358)
imagesc(group_highflip_2_new(:,:,8), [0 600]); colorbar
subplot (359)
imagesc(group_highflip_2_new(:,:,9), [0 600]); colorbar
subplot (3,5,10)
imagesc(group_highflip_2_new(:,:,10), [0 600]); colorbar
subplot (3,5,11)
imagesc(group_highflip_2_new(:,:,11), [0 600]); colorbar
subplot (3,5,12)
imagesc(group_highflip_2_new(:,:,12), [0 600]); colorbar
subplot (3,5,13)
imagesc(group_highflip_2_new(:,:,13), [0 600]); colorbar
subplot (3,5,14)
imagesc(group_highflip_2_new(:,:,14), [0 600]); colorbar
subplot (3,5,15)
imagesc(group_highflip_2_new(:,:,15), [0 600]); colorbar;
sgtitle('MSLOW0015:  Gre Highflip gefiltert und rekonstruiert nach T2*, TE = 11 ms');

figure (3)
subplot (351)
imagesc(varianz_Map(:,:,1) ); colorbar
subplot (352)
imagesc(varianz_Map(:,:,2) ); colorbar
subplot (353)
imagesc(varianz_Map(:,:,3) ); colorbar
subplot (354)
imagesc(varianz_Map(:,:,4) ); colorbar
subplot (355)
imagesc(varianz_Map(:,:,5) ); colorbar
subplot (356)
imagesc(varianz_Map(:,:,6) ); colorbar
subplot (357)
imagesc(varianz_Map(:,:,7) ); colorbar
subplot (358)
imagesc(varianz_Map(:,:,8) ); colorbar
subplot (359)
imagesc(varianz_Map(:,:,9) ); colorbar
subplot (3,5,10)
imagesc(varianz_Map(:,:,10) ); colorbar
subplot (3,5,11)
imagesc(varianz_Map(:,:,11) ); colorbar
subplot (3,5,12)
imagesc(varianz_Map(:,:,12) ); colorbar
subplot (3,5,13)
imagesc(varianz_Map(:,:,13) ); colorbar
subplot (3,5,14)
imagesc(varianz_Map(:,:,14) ); colorbar
subplot (3,5,15)
imagesc(varianz_Map(:,:,15) ); colorbar;
sgtitle('MSLOW0015:  Varianz zwischen neu erstellten MR-Bildern');