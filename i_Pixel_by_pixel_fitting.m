clear all
close all
clc

%% Required Input data
Figure_TitleName1 = ['R1 und S0, Slope Martrix, WM'];
Figure_TitleName2 = ['R1 und S0, Offset Martrix, WM'];
Figure_TitleName3 = ['R2 und S0, Slope Martrix, WM'];
Figure_TitleName4 = ['R2 und S0, Offset Martrix, WM'];

GM_struct = load('GM.mat');                  % See: h_Segmentierung.m 
WM_struct = load('WM.mat');                  % See: h_Segmentierung.m        
S0_struct = load('S0_Map.mat');              % See: f_Kombinierte_Denosierung_auf_S0.m
T1_struct = load('T1_nachdem_Kombinierte_Denosierung.mat');
                                             % See: g_T1_nachdem_Kombinierte_Denosierung.m
T2_struct = load('T2Sternmap.mat'); 
                                             % See: c_T2Stern_Map_Creating  
%% Output Informations  
savename1 = ['R1-SO,Slope,WM.jpg'];          % Save the pixel-by-pixel fitted matrix as picture.
savename2 = ['R1-SO,Offset,WM.jpg'];
savename3 = ['R2-SO,Slope,WM.jpg'];
savename4 = ['R2-SO,Offset,WM.jpg'];

%% Main Code
%The split indication matrix is in Struck format, and the desired logical matrix should be obtained in the following ways:                                            
Pointer0 = fieldnames(GM_struct);
Pointer1 = fieldnames(WM_struct);
Pointer2 = fieldnames(S0_struct);
Pointer3 = fieldnames(T1_struct);
Pointer4 = fieldnames(T2_struct);
GM = getfield(GM_struct,char(Pointer0));
WM = getfield(WM_struct,char(Pointer1));
S0 = getfield(S0_struct,char(Pointer2));
T1 = getfield(T1_struct,char(Pointer3));
T2 = getfield(T2_struct,char(Pointer4)); 

T1_Map       = T1 .* WM;                     %  Split the required WM or GM data on the target data
S0_Map       = S0 .* WM;
T2_Stern_Map = T2 .* WM;

%% Exclude elements with a value of "Inf" in R1 and R2
R1 = 1/T1_Map;
[m1,n1,p1] = size(T1_Map);
    for i11 = 1:m1
        for j11 = 1:n1
            for k11 = 1:p1
                if isinf(R1(i11,j11,k11))
                    R1(i11,j11,k11) = 0;  
                end
            end
        end
    end

S0 = S0_Map;

R2 = 1/T2_Stern_Map;
    for i11 = 1:m1
        for j11 = 1:n1
            for k11 = 1:p1
                if isinf(R2(i11,j11,k11))
                    R2(i11,j11,k11) = 0;  
                end
            end
        end
    end


%% Carry out element-by-element linear regression
n = 5   %5*5像素点,为奇数
[a,b,c] = size(R1);
para1 = zeros(a,b,c);                            % New 0 matrix, the slope obtained  by R1 and S0 regression
para2 = zeros(a,b,c);                            % New 0 matrix, the offset obtained by R1 and S0 regression
para3 = zeros(a,b,c);                            % New 0 matrix, the slope obtained  by R2 and S0 regression
para4 = zeros(a,b,c);                            % New 0 matrix, the offset obtained by R2 and S0 regression
nn = (n-1)/2;
yeah = 0;
for k = 1:15
    for inn = nn+1:a-nn
        for j = nn+1:b-nn
                x1 = reshape(R1(inn-nn:inn+nn,j-nn:j+nn,k),1,[]);
                x2 = reshape(R2(inn-nn:inn+nn,j-nn:j+nn,k),1,[]);
                y  = reshape(S0(inn-nn:inn+nn,j-nn:j+nn,k),1,[]);
            if sum(sum(x1~=0))>5
                 A = polyfit(x1,y,1);             % slope and offset, fitted by R1 - S0
                 B = polyfit(x2,y,1);             % slope and offset, fitted by R2 - S0
                 para1(inn,j,k) = A(1,1);         % slope,  obtained by R1 and S0 regression
                 para2(inn,j,k) = A(1,2);         % offset, obtained by R1 and S0 regression
                 para3(inn,j,k) = B(1,1);         % slope,  obtained by R2 and S0 regression
                 para4(inn,j,k) = B(1,2);         % offset, obtained by R2 and S0 regression
                 yeah = yeah+1;                   % It is used to count the number of efficient points.
            end
         end
    end
end




figure ()
subplot (351)
imagesc(para1(:,:,1), [0 200000]); colormap(gray); colorbar
subplot (352)
imagesc(para1(:,:,2), [0 200000]); colormap(gray); colorbar
subplot (353)
imagesc(para1(:,:,3), [0 200000]); colormap(gray); colorbar
subplot (354)
imagesc(para1(:,:,4), [0 200000]); colormap(gray); colorbar
subplot (355)
imagesc(para1(:,:,5), [0 200000]); colormap(gray); colorbar
subplot (356)
imagesc(para1(:,:,6), [0 200000]); colormap(gray); colorbar
subplot (357)
imagesc(para1(:,:,7), [0 200000]); colormap(gray); colorbar
subplot (358)
imagesc(para1(:,:,8), [0 200000]); colormap(gray); colorbar
subplot (359)
imagesc(para1(:,:,9), [0 200000]); colormap(gray); colorbar
subplot (3,5,10)
imagesc(para1(:,:,10), [0 200000]); colormap(gray); colorbar
subplot (3,5,11)
imagesc(para1(:,:,11), [0 200000]); colormap(gray); colorbar
subplot (3,5,12)
imagesc(para1(:,:,12), [0 200000]); colormap(gray); colorbar
subplot (3,5,13)
imagesc(para1(:,:,13), [0 200000]); colormap(gray); colorbar
subplot (3,5,14)
imagesc(para1(:,:,14), [0 200000]); colormap(gray); colorbar
subplot (3,5,15)
imagesc(para1(:,:,15), [0 200000]); colormap(gray); colorbar;
sgtitle(Figure_TitleName1);
saveas(gcf,savename1)

figure ()
subplot (351)
imagesc(para2(:,:,1), [0 500]); colormap(gray); colorbar
subplot (352)
imagesc(para2(:,:,2), [0 500]); colormap(gray); colorbar
subplot (353)
imagesc(para2(:,:,3), [0 500]); colormap(gray); colorbar
subplot (354)
imagesc(para2(:,:,4), [0 500]); colormap(gray); colorbar
subplot (355)
imagesc(para2(:,:,5), [0 500]); colormap(gray); colorbar
subplot (356)
imagesc(para2(:,:,6), [0 500]); colormap(gray); colorbar
subplot (357)
imagesc(para2(:,:,7), [0 500]); colormap(gray); colorbar
subplot (358)
imagesc(para2(:,:,8), [0 500]); colormap(gray); colorbar
subplot (359)
imagesc(para2(:,:,9), [0 500]); colormap(gray); colorbar
subplot (3,5,10)
imagesc(para2(:,:,10), [0 500]); colormap(gray); colorbar
subplot (3,5,11)
imagesc(para2(:,:,11), [0 500]); colormap(gray); colorbar
subplot (3,5,12)
imagesc(para2(:,:,12), [0 500]); colormap(gray); colorbar
subplot (3,5,13)
imagesc(para2(:,:,13), [0 500]); colormap(gray); colorbar
subplot (3,5,14)
imagesc(para2(:,:,14), [0 500]); colormap(gray); colorbar
subplot (3,5,15)
imagesc(para2(:,:,15), [0 500]); colormap(gray); colorbar;
sgtitle( Figure_TitleName2);
saveas(gcf,savename2)




figure()
subplot (351)
imagesc(para3(:,:,1), [0 50000]); colormap(gray); colorbar
subplot (352)
imagesc(para3(:,:,2), [0 50000]); colormap(gray); colorbar
subplot (353)
imagesc(para3(:,:,3), [0 50000]); colormap(gray); colorbar
subplot (354)
imagesc(para3(:,:,4), [0 50000]); colormap(gray); colorbar
subplot (355)
imagesc(para3(:,:,5), [0 50000]); colormap(gray); colorbar
subplot (356)
imagesc(para3(:,:,6), [0 50000]); colormap(gray); colorbar
subplot (357)
imagesc(para3(:,:,7), [0 50000]); colormap(gray); colorbar
subplot (358)
imagesc(para3(:,:,8), [0 50000]); colormap(gray); colorbar
subplot (359)
imagesc(para3(:,:,9), [0 50000]); colormap(gray); colorbar
subplot (3,5,10)
imagesc(para3(:,:,10), [0 50000]); colormap(gray); colorbar
subplot (3,5,11)
imagesc(para3(:,:,11), [0 50000]); colormap(gray); colorbar
subplot (3,5,12)
imagesc(para3(:,:,12), [0 50000]); colormap(gray); colorbar
subplot (3,5,13)
imagesc(para3(:,:,13), [0 50000]); colormap(gray); colorbar
subplot (3,5,14)
imagesc(para3(:,:,14), [0 50000]); colormap(gray); colorbar
subplot (3,5,15)
imagesc(para3(:,:,15), [0 50000]); colormap(gray); colorbar;
sgtitle(Figure_TitleName3);
saveas(gcf,savename3)

figure ()
subplot (351)
imagesc(para4(:,:,1), [0 400]); colormap(gray); colorbar
subplot (352)
imagesc(para4(:,:,2), [0 400]); colormap(gray); colorbar
subplot (353)
imagesc(para4(:,:,3), [0 400]); colormap(gray); colorbar
subplot (354)
imagesc(para4(:,:,4), [0 400]); colormap(gray); colorbar
subplot (355)
imagesc(para4(:,:,5), [0 400]); colormap(gray); colorbar
subplot (356)
imagesc(para4(:,:,6), [0 400]); colormap(gray); colorbar
subplot (357)
imagesc(para4(:,:,7), [0 400]); colormap(gray); colorbar
subplot (358)
imagesc(para4(:,:,8), [0 400]); colormap(gray); colorbar
subplot (359)
imagesc(para4(:,:,9), [0 400]); colormap(gray); colorbar
subplot (3,5,10)
imagesc(para4(:,:,10), [0 400]); colormap(gray); colorbar
subplot (3,5,11)
imagesc(para4(:,:,11), [0 400]); colormap(gray); colorbar
subplot (3,5,12)
imagesc(para4(:,:,12), [0 400]); colormap(gray); colorbar
subplot (3,5,13)
imagesc(para4(:,:,13), [0 400]); colormap(gray); colorbar
subplot (3,5,14)
imagesc(para4(:,:,14), [0 400]); colormap(gray); colorbar
subplot (3,5,15)
imagesc(para4(:,:,15), [0 400]); colormap(gray); colorbar;
sgtitle(Figure_TitleName4);
saveas(gcf,savename4)


