clear all
close all
clc



%% Required Input data
DicomDatei_Lowflip           = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_lowflip\'];
DicomDatei_flipmap           = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\FlipMaps\'];
DicomDatei_Highflip          = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_highflip_orig\'];

load('ReferenzMatrix.mat')                  % See: b_T1_Hirnmaskierung.m
load('T2Sternmap.mat')                      % See: c_T2Stern_Map_Creating.m                    

%% Output Informations  
TargetData_S1S2_Raumlichen_Denoising_Info = ['Highflip_Datei_nachdem_Raumlichen_Denoising.mat']
                                            % The calculated results will be saved to the same folder with this file name.
TargetData_S0_After_Raumlichen_Denoising  = ['S0_Map.mat']

%% Main Codes
% DICOM_Matrix einladen
 
[M_highflip]                = loadDicomFolder(DicomDatei_Highflip,-1,0);
Gre_highflip_group1         = M_highflip(:,:, 1:15).*ReferenzMatrix;
                                            % Hirnmaskierung
Gre_highflip_group2         = M_highflip(:,:,16:30).*ReferenzMatrix;

% DICOM_Info auslesen
highflip_info_S1            = dicominfo([DicomDatei_Highflip,'File0001.ima']) 
                                            % FlipAngle: 14; RepetitionTime: 200; EchoTime: 3.9500
highflip_info_S2            = dicominfo([DicomDatei_Highflip,'File0016.ima']) 
                                            % FlipAngle: 14; RepetitionTime: 200; EchoTime: 11
TE1                         = highflip_info_S1.EchoTime;
TE2                         = highflip_info_S2.EchoTime;
        
        
[nrow,ncol,nSlices]         = size(Gre_highflip_group1);
group_highflip_1_new_raum_8nachbar = zeros(nrow,ncol,nSlices);
group_highflip_2_new_raum_8nachbar = zeros(nrow,ncol,nSlices);
S0_Map                      = zeros(nrow,ncol,nSlices);
S0_new                      = zeros(nrow,ncol,nSlices);
var_Map                     = zeros(nrow,ncol,nSlices);
S0_Region                   = zeros(3,3,1);

for slice = 1:1:nSlices                      % Schichten
    for row = 2:1:nrow                       % Zeilen
        for col = 2:1:ncol                   % Spalten
            if Gre_highflip_group1(row,col,slice) ~=0
                x1          = Gre_highflip_group1(row-1,col-1,slice);  
                                             % (x-1, y-1), usw... 
                                             % 18 Nachbarschaft
                x2          = Gre_highflip_group1(row,col-1,slice);
                x3          = Gre_highflip_group1(row+1,col-1,slice);
                x4          = Gre_highflip_group1(row-1,col,slice); 
                x5          = Gre_highflip_group1(row,col,slice); 
                x6          = Gre_highflip_group1(row+1,col,slice); 
                x7          = Gre_highflip_group1(row-1,col+1,slice);
                x8          = Gre_highflip_group1(row,col+1,slice);
                x9          = Gre_highflip_group1(row+1,col+1,slice); 
                x10         = Gre_highflip_group2(row-1,col-1,slice);
                x11         = Gre_highflip_group2(row,col-1,slice);
                x12         = Gre_highflip_group2(row+1,col-1,slice);
                x13         = Gre_highflip_group2(row-1,col,slice); 
                x14         = Gre_highflip_group2(row,col,slice); 
                x15         = Gre_highflip_group2(row+1,col,slice); 
                x16         = Gre_highflip_group2(row-1,col+1,slice);
                x17         = Gre_highflip_group2(row,col+1,slice);
                x18         = Gre_highflip_group2(row+1,col+1,slice); 
                y1          = exp((-TE1)/T2_Stern_Map(row-1,col-1,slice));
                y2          = exp((-TE1)/T2_Stern_Map(row,col-1,slice));
                y3          = exp((-TE1)/T2_Stern_Map(row+1,col-1,slice));
                y4          = exp((-TE1)/T2_Stern_Map(row-1,col,slice));
                y5          = exp((-TE1)/T2_Stern_Map(row,col,slice));
                y6          = exp((-TE1)/T2_Stern_Map(row+1,col,slice));
                y7          = exp((-TE1)/T2_Stern_Map(row-1,col+1,slice));
                y8          = exp((-TE1)/T2_Stern_Map(row,col+1,slice));
                y9          = exp((-TE1)/T2_Stern_Map(row+1,col+1,slice));
                y10         = exp((-TE2)/T2_Stern_Map(row-1,col-1,slice));
                y11         = exp((-TE2)/T2_Stern_Map(row,col-1,slice));
                y12         = exp((-TE2)/T2_Stern_Map(row+1,col-1,slice));
                y13         = exp((-TE2)/T2_Stern_Map(row-1,col,slice));
                y14         = exp((-TE2)/T2_Stern_Map(row,col,slice));
                y15         = exp((-TE2)/T2_Stern_Map(row+1,col,slice));
                y16         = exp((-TE2)/T2_Stern_Map(row-1,col+1,slice));
                y17         = exp((-TE2)/T2_Stern_Map(row,col+1,slice));
                y18         = exp((-TE2)/T2_Stern_Map(row+1,col+1,slice));
                
                numerator   = (x1*y1)+(x2*y2)+(x3*y3)+(x4*y4)+(x5*y5)+(x6*y6)+(x7*y7)+(x8*y8)+(x9*y9)+(x10*y10)+(x11*y11)+(x12*y12)+(x13*y13)+(x14*y14)+(x15*y15)+(x16*y16)+(x17*y17)+(x18*y18);
                denominator = (y1*y1)+(y2*y2)+(y3*y3)+(y4*y4)+(y5*y5)+(y6*y6)+(y7*y7)+(y8*y8)+(y9*y9)+(y10*y10)+(y11*y11)+(y12*y12)+(y13*y13)+(y14*y14)+(y15*y15)+(y16*y16)+(y17*y17)+(y18*y18);
                S0_Equ      = numerator/denominator;
                S0_Map(row,col,slice)                             = S0_Equ;
                group_highflip_1_new_raum_8nachbar(row,col,slice) = round(S0_Equ*exp((-TE1)/T2_Stern_Map(row,col,slice)));%S1 new;
                group_highflip_2_new_raum_8nachbar(row,col,slice) = round(S0_Equ*exp((-TE2)/T2_Stern_Map(row,col,slice)));%S2 new;
                var_Map(row,col,slice)                            = (group_highflip_1_new_raum_8nachbar(row,col,slice)-(S0_Equ*exp((-TE1)/T2_Stern_Map(row,col,slice))))^2 + (group_highflip_2_new_raum_8nachbar(row,col,slice)-(S0_Equ*exp((-TE2)/T2_Stern_Map(row,col,slice))))^2;
             end
        end
    end
end

%% Rename the calculated data and export it.  The exported file name can be customized at the beginning of the file.

save(TargetData_S1S2_Raumlichen_Denoising_Info,'group_highflip_1_new_raum_8nachbar','group_highflip_2_new_raum_8nachbar')               
                                         % Save the new created S1 and S2 as mat form in the same folder, which can be imported in the next step.

save(TargetData_S0_After_Raumlichen_Denoising,'S0_Map')
                                         % Save the new created S0 Map as mat form in the same folder, which can be imported in the next step.
figure (1)
subplot (351)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,1), [0 800]); colorbar
subplot (352)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,2), [0 800]); colorbar
subplot (353)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,3), [0 800]); colorbar
subplot (354)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,4), [0 800]); colorbar
subplot (355)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,5), [0 800]); colorbar
subplot (356)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,6), [0 800]); colorbar
subplot (357)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,7), [0 800]); colorbar
subplot (358)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,8), [0 800]); colorbar
subplot (359)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,9), [0 800]); colorbar
subplot (3,5,10)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,10), [0 800]); colorbar
subplot (3,5,11)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,11), [0 800]); colorbar
subplot (3,5,12)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,12), [0 800]); colorbar
subplot (3,5,13)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,13), [0 800]); colorbar
subplot (3,5,14)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,14), [0 800]); colorbar
subplot (3,5,15)
imagesc(group_highflip_1_new_raum_8nachbar(:,:,15), [0 800]); colorbar;
sgtitle('MSLOW0015:  Gre Highflip gefiltert und rekonstruiert nach Raumlichen Filterung, TE = 3.95ms');


figure (2)
subplot (351)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,1), [0 600]); colorbar
subplot (352)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,2), [0 600]); colorbar
subplot (353)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,3), [0 600]); colorbar
subplot (354)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,4), [0 600]); colorbar
subplot (355)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,5), [0 600]); colorbar
subplot (356)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,6), [0 600]); colorbar
subplot (357)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,7), [0 600]); colorbar
subplot (358)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,8), [0 600]); colorbar
subplot (359)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,9), [0 600]); colorbar
subplot (3,5,10)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,10), [0 600]); colorbar
subplot (3,5,11)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,11), [0 600]); colorbar
subplot (3,5,12)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,12), [0 600]); colorbar
subplot (3,5,13)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,13), [0 600]); colorbar
subplot (3,5,14)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,14), [0 600]); colorbar
subplot (3,5,15)
imagesc(group_highflip_2_new_raum_8nachbar(:,:,15), [0 600]); colorbar;
sgtitle('MSLOW0015:  Gre Highflip gefiltert und rekonstruiert nach Raumlichen Filterung, TE = 11 ms');

figure (3)
subplot (351)
imagesc(var_Map(:,:,1) ); colorbar
subplot (352)
imagesc(var_Map(:,:,2) ); colorbar
subplot (353)
imagesc(var_Map(:,:,3) ); colorbar
subplot (354)
imagesc(var_Map(:,:,4) ); colorbar
subplot (355)
imagesc(var_Map(:,:,5) ); colorbar
subplot (356)
imagesc(var_Map(:,:,6) ); colorbar
subplot (357)
imagesc(var_Map(:,:,7) ); colorbar
subplot (358)
imagesc(var_Map(:,:,8) ); colorbar
subplot (359)
imagesc(var_Map(:,:,9) ); colorbar
subplot (3,5,10)
imagesc(var_Map(:,:,10) ); colorbar
subplot (3,5,11)
imagesc(var_Map(:,:,11) ); colorbar
subplot (3,5,12)
imagesc(var_Map(:,:,12) ); colorbar
subplot (3,5,13)
imagesc(var_Map(:,:,13) ); colorbar
subplot (3,5,14)
imagesc(var_Map(:,:,14) ); colorbar
subplot (3,5,15)
imagesc(var_Map(:,:,15) ); colorbar;
sgtitle('MSLOW0015:  Varianz zwischen neu erstellten MR-Bildern nach Raumlichen Filterung');