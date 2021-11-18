clear all
close all
clc
% This code is mainly used to calculate T2 * Map
%% Required Input data
DicomDatei_Lowflip         = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_lowflip\'];
DicomDatei_flipmap         = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\FlipMaps\'];
DicomDatei_Highflip        = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_highflip_orig\'];

%% Output Informations  
TargetData_T2_Stern_Info   = ['T2Sternmap.mat']
                                            % The calculated results will be saved to the same folder with this file name.


%% Main Codes                                           
% Two different T1 is obtained in two steps, and the average value is calculated, so that a more accurate T1 can be obtained.        
[M_flipmap]                = loadDicomFolder(DicomDatei_flipmap,-1,0);
[M_lowflip]                = loadDicomFolder(DicomDatei_Lowflip,-1,0); 
[M_highflip]               = loadDicomFolder(DicomDatei_Highflip,-1,0);

%--Automatically collect all TE data: BEGIN
[curDS,curNoFiles,infoFirst,DcmInfoList]  = loadDicomFolder(DicomDatei_Lowflip,-1,0);
i_TE = 1;
for group = 1:15:420
    TE_group(1,i_TE)       = DcmInfoList{group,1}.EchoTime;
    i_TE                   = i_TE+1;
end
%--Automatically collect all TE data: END

group_lowflip    = M_lowflip(:,:, 1:15);   % Extract the first 15 dicom files, that is, a complete brain measurement

%T2 herausfinden
[x,y,z_M]               = size(M_lowflip);
Signal_Voxel            = zeros();         % Collect each voxel in the same location and create a four-dimensional matrix as a temporary container. 
                                           % The next step is to make a linear regression of the voxels in the same position.
Signal_Voxel_polyfit    = zeros();         % temporary container to make linear regression
Ln_Signal_Voxel_polyfit = zeros();         % Find the natural logarithm of the value of each element of A
T2_Stern_Map            = zeros(size(group_lowflip));

% Collect each voxel in the same position to create a four-dimensional matrix
for d = 1:15                               % Das ganze Gehirn wird in 15 Schichten geschnitten
        for a = 2:1:x                      % Zeilen
            for b = 2:1:y                  % Spalten
                    e = 1;                          
                    for c = d:15:z_M       % Gleiches Voxel unter verschiedenen TE
                        Signal_Voxel(d,a,b,e) = M_lowflip(a,b,c) ; 
                        e = e+1;           
                    end
            end
        
        end
d = d+1;
end

% T2* Mapping
noise_lowfilp               = group_lowflip (7:11,7:11,1); 
                                           % Random background noise: 3*3 Matrix
NoiseLevel_lowflip          = mean(noise_lowfilp(:));
noise_Verstaerkung          = 5;
denoise_parameter_low       = NoiseLevel_lowflip * noise_Verstaerkung;

for d = 1:15   
    for a = 2:1:x
        for b = 2:1:y
            if  M_lowflip(a,b,d) > denoise_parameter_low
                Signal_Voxel_polyfit            = permute(Signal_Voxel(d,a,b,:),[1 4 2 3]);
                Ln_Signal_Voxel_polyfit         = log(Signal_Voxel_polyfit);
                                           % Find the natural logarithm of the value of each element of A
                    for g = 1:28
                        if  Ln_Signal_Voxel_polyfit(1,g) <0
                            Ln_Signal_Voxel_polyfit(1,g) = 0;
                        end
                            p                    = polyfit(TE_group,Ln_Signal_Voxel_polyfit,1);
                                           % The negative value will be calculated here, which is not consistent with the actual situation, and the absolute value will be taken.
                            R2_Stern             = abs(p(1,1));
                            T2_Stern             = 1/R2_Stern;
                            T2_Stern_Map(a,b,d)  = T2_Stern;
                            f                    = polyval(p,TE_group);
                    end
            end
        end
    end
d = d+1                                    % The calculation process is long. Use d to show the progress of the calculation
end

%% Rename the calculated data and export it.  The exported file name can be customized at the beginning of the file.
save(TargetData_T2_Stern_Info,'T2_Stern_Map')              
                                           % Save the target T2* Map as mat form in the same folder, which can be imported in the next step.

figure (2)
plot(TE_group,Ln_Signal_Voxel_polyfit,'o')
hold on
plot(TE_group,f,'r--')
axis([0  70  0  10])
xlabel('TE(ms)') 
ylabel('ln(Signal)') 
legend('ln(S(TE))','Angepasste Gerade')
title('Lineare Regression mit 28 Datenpunkten zur Anpassung an T2*');


figure (1)
subplot (351)
imagesc(T2_Stern_Map(:,:,1), [0 600]); colorbar
subplot (352)
imagesc(T2_Stern_Map(:,:,2), [0 600]); colorbar
subplot (353)
imagesc(T2_Stern_Map(:,:,3), [0 600]); colorbar
subplot (354)
imagesc(T2_Stern_Map(:,:,4), [0 600]); colorbar
subplot (355)
imagesc(T2_Stern_Map(:,:,5), [0 600]); colorbar
subplot (356)
imagesc(T2_Stern_Map(:,:,6), [0 600]); colorbar
subplot (357)
imagesc(T2_Stern_Map(:,:,7), [0 600]); colorbar
subplot (358)
imagesc(T2_Stern_Map(:,:,8), [0 600]); colorbar
subplot (359)
imagesc(T2_Stern_Map(:,:,9), [0 600]); colorbar
subplot (3,5,10)
imagesc(T2_Stern_Map(:,:,10), [0 600]); colorbar
subplot (3,5,11)
imagesc(T2_Stern_Map(:,:,11), [0 600]); colorbar
subplot (3,5,12)
imagesc(T2_Stern_Map(:,:,12), [0 600]); colorbar
subplot (3,5,13)
imagesc(T2_Stern_Map(:,:,13), [0 600]); colorbar
subplot (3,5,14)
imagesc(T2_Stern_Map(:,:,14), [0 600]); colorbar
subplot (3,5,15)
imagesc(T2_Stern_Map(:,:,15), [0 600]); colorbar;
sgtitle('MSLOW0015:  T2 Stern Relaxationszeit');