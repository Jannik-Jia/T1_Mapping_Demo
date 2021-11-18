clear all
close all
clc
%% Split the gray matter and white matter of the head

%% Required Input data
DicomDatei_Lowflip              = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_lowflip\'];
DicomDatei_flipmap              = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\FlipMaps\'];
DicomDatei_Highflip             = ['C:\Users\Jannik Jia\OneDrive - rheinahrcampus.de\Bachelorarbeit\T1Optimisation\MSLOW0015\gre_highflip_orig\'];
load('T1_nachdem_Kombinierte_Denosierung.mat')                  % See: g_T1_nachdem_Kombinierte_Denosierung
load('T2Sternmap.mat')                                          % See: d_T2_Stern_Filterung.m   
                                                                
%% Output Informations  
TargetData_GM_LogicMatrix       = ['GM.mat'];  
                                                                % The calculated results will be saved to the same folder with this file name.
TargetData_WM_LogicMatrix       = ['WM.mat']

%% Main Codes  
[greLowFlip,iNoFiles,infoFirst] = loadDicomFolder(DicomDatei_Lowflip,-1,0);
T1Maps = T1_After_Denoising_Raum_8nachbar;
T2Maps = T2_Stern_Map;

%--
[nrow,ncol,nSlices] = size(T1Maps);
SE = strel('disk',15);
graymatte = zeros(nrow,ncol,nSlices);
whitematte = zeros(nrow,ncol,nSlices);
greLow1   = greLowFlip(:,:,1:15);

%% Brain Mask
SE1=strel('disk',3);
for i=1:nSlices
    maske(:,:,i)=greLowFlip(:,:,i+nSlices*7)./greLowFlip(:,:,nSlices*3+i);
    maske(:,:,i)=maske(:,:,i)>1.5;
    maske(:,:,i)=imclose(maske(:,:,i),SE);
    a=~maske(:,:,i);
    L=bwlabel(a);
    stats = regionprops(L);
    Ar = cat(1, stats.Area);
    ind = find(Ar ==max(Ar));
    a(find(L~=ind))= 0;
    maske(:,:,i)=a;
    maske(:,:,i)=imclose(maske(:,:,i),SE);
    maske(:,:,i)=imerode(a,SE1);
    greLowFlip(:,:,i) = greLowFlip(:,:,i).*maske(:,:,i);
end

T1Maps = T1Maps.*maske;

%% Segmentation CSF, WM, GM
CSF = T1Maps>1600; %best with M0 filter
%CSF = T1Maps>2200; %default: 2300 best:1400-1500
WM = T1Maps <=400 & T1Maps>=100 & T2Maps>60;
GM = T1Maps <900 & T1Maps>500 & T2Maps>60;

whitematte = greLow1.*WM;
graymatte  = greLow1.*GM;

%% Rename the calculated data and export it.  The exported file name can be customized at the beginning of the file.

save(TargetData_WM_LogicMatrix,'WM')               % Save the target WM/GM Logic Matrix as mat form in the same folder, which can be imported in the next step.
save(TargetData_GM_LogicMatrix,'GM') 


%% Bild Vergleichung
figure (1)
subplot (351)
imagesc(whitematte(:,:,1) ); colorbar
subplot (352)
imagesc(whitematte(:,:,2) ); colorbar
subplot (353)
imagesc(whitematte(:,:,3) ); colorbar
subplot (354)
imagesc(whitematte(:,:,4) ); colorbar
subplot (355)
imagesc(whitematte(:,:,5) ); colorbar
subplot (356)
imagesc(whitematte(:,:,6) ); colorbar
subplot (357)
imagesc(whitematte(:,:,7) ); colorbar
subplot (358)
imagesc(whitematte(:,:,8) ); colorbar
subplot (359)
imagesc(whitematte(:,:,9) ); colorbar
subplot (3,5,10)
imagesc(whitematte(:,:,10) ); colorbar
subplot (3,5,11)
imagesc(whitematte(:,:,11) ); colorbar
subplot (3,5,12)
imagesc(whitematte(:,:,12) ); colorbar
subplot (3,5,13)
imagesc(whitematte(:,:,13) ); colorbar
subplot (3,5,14)
imagesc(whitematte(:,:,14) ); colorbar
subplot (3,5,15)
imagesc(whitematte(:,:,15) ); colorbar;
sgtitle('MSLOW0015:  Whitematter');

figure (2)
subplot (351)
imagesc(graymatte(:,:,1) ); colorbar
subplot (352)
imagesc(graymatte(:,:,2) ); colorbar
subplot (353)
imagesc(graymatte(:,:,3) ); colorbar
subplot (354)
imagesc(graymatte(:,:,4) ); colorbar
subplot (355)
imagesc(graymatte(:,:,5) ); colorbar
subplot (356)
imagesc(graymatte(:,:,6) ); colorbar
subplot (357)
imagesc(graymatte(:,:,7) ); colorbar
subplot (358)
imagesc(graymatte(:,:,8) ); colorbar
subplot (359)
imagesc(graymatte(:,:,9) ); colorbar
subplot (3,5,10)
imagesc(graymatte(:,:,10) ); colorbar
subplot (3,5,11)
imagesc(graymatte(:,:,11) ); colorbar
subplot (3,5,12)
imagesc(graymatte(:,:,12) ); colorbar
subplot (3,5,13)
imagesc(graymatte(:,:,13) ); colorbar
subplot (3,5,14)
imagesc(graymatte(:,:,14) ); colorbar
subplot (3,5,15)
imagesc(graymatte(:,:,15) ); colorbar;
sgtitle('MSLOW0015:  Graymatter');


% figure (1)
% subplot 321
% imagesc(GM(:,:,1)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,1)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,1), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,1), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,1), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,1), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% figure (2)
% subplot 321
% imagesc(GM(:,:,2)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,2)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,2), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,2), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,2), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,2), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% 
% 
% 
% figure (3)
% subplot 321
% imagesc(GM(:,:,3)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,3)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,3), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,3), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,3), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,3), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% figure (4)
% subplot 321
% imagesc(GM(:,:,4)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,4)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,4), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,4), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,4), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,4), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% figure (5)
% subplot 321
% imagesc(GM(:,:,5)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,5)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,5), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,5), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,5), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,5), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% 
% figure (6)
% subplot 321
% imagesc(GM(:,:,6)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,6)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,6), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,6), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,6), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,6), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% 
% figure (7)
% subplot 321
% imagesc(GM(:,:,7)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,7)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,7), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,7), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,7), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,7), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% figure (8)
% subplot 321
% imagesc(GM(:,:,8)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,8)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,8), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,8), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,8), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,8), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% figure (9)
% subplot 321
% imagesc(GM(:,:,9)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,9)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,9), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,9), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,9), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,9), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% figure (10)
% subplot 321
% imagesc(GM(:,:,10)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,10)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,10), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,10), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,10), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,10), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% figure (11)
% subplot 321
% imagesc(GM(:,:,11)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,11)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,11), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,11), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,11), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,11), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% figure (12)
% subplot 321
% imagesc(GM(:,:,12)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,12)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,12), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,12), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,12), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,12), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% figure (13)
% subplot 321
% imagesc(GM(:,:,13)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,13)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,13), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,13), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,13), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,13), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
% 
% figure (14)
% subplot 321
% imagesc(GM(:,:,14)),colorbar
% title('Logische Matrix der graue Substanz, 8 Nachbarschaft')
% subplot 322
% imagesc(WM(:,:,14)),colorbar
% title('Logische Matrix der weißen Substanz, 8 Nachbarschaft')
% subplot 323
% graymatte = greLow1.*GM;
% imagesc(graymatte(:,:,14), [0 2000]),colorbar 
% title('Graue Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 324
% whitematte = greLow1.*WM;
% imagesc(whitematte(:,:,14), [0 2000]),colorbar 
% title('Weiße Substanz nach Segmentierung, 8 Nachbarschaft')
% subplot 325
% imagesc(greLow1(:,:,14), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% subplot 326
% imagesc(greLow1(:,:,14), [0 2000]),colorbar 
% title('Bild des Gehirns als Vergleich')
% 
