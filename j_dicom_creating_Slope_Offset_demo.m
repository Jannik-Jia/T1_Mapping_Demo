        
% for i = 1:12
        i = 1   %COLOW0001
        %i = 2   %COLOW0002
        %...
        %i = 12  %COLOW0012
        
        %files automaticlly reading
        GM_file = load(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\Datei_JJia\GM_CO_000',num2str(i),'.mat']);
        S0_file = load(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\Datei_JJia\S0_Map_CO_000',num2str(i),'.mat']);
        T1_file = load(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\Datei_JJia\T1_After_Denoising_Raum_CO_000',num2str(i),'.mat']);
        T2_file = load(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\Datei_JJia\T2_Stern_Map_CO_000',num2str(i),'.mat']);

        name1 = fieldnames(GM_file);
        name2 = fieldnames(S0_file);
        name3 = fieldnames(T1_file);
        name4 = fieldnames(T2_file);
        GM = getfield(GM_file,char(name1)); %Logical Matrix from GM
        S0 = getfield(S0_file,char(name2));
        T1 = getfield(T1_file,char(name3));
        T2 = getfield(T2_file,char(name4));


        T1_Map       = T1 .* GM;
        S0_Map       = S0 .* GM;
        T2_Stern_Map = T2 .* GM;
%         T1_Map       = T1 .* WM;
%         S0_Map       = S0 .* WM;
%         T2_Stern_Map = T2 .* WM;


        %% Replace Inf as 0
        R1=1/T1_Map;
        [m1,n1,p1]=size(T1_Map);
            for i11=1:m1
                for j11=1:n1
                    for k11=1:p1
                        if isinf(R1(i11,j11,k11))
                            R1(i11,j11,k11)=0;  
                        end
                    end
                end
            end
        
        S0=S0_Map;

        R2=1/T2_Stern_Map;
            for i11=1:m1
                for j11=1:n1
                    for k11=1:p1
                        if isinf(R2(i11,j11,k11))
                            R2(i11,j11,k11)=0;  
                        end
                    end
                end
            end
        
        
        %% ployfit. Create a 5-to-5 matrix. Linear regression of the values of this elements
        n=5                 %5*5 Matrix
        [a,b,c]=size(R1);
        para1=zeros(a,b,c); % Slope between R1 and S0
        para2=zeros(a,b,c); % Offset between R1 and S0
        para3=zeros(a,b,c); % Slope between R2 and S0
        para4=zeros(a,b,c); % Offset between R2 and S0
        nn=(n-1)/2;
        yeah=0;
        for k=1:15
            for inn=nn+1:a-nn
                for j=nn+1:b-nn
                    x1=reshape(R1(inn-nn:inn+nn,j-nn:j+nn,k),1,[]);
                    x2=reshape(R2(inn-nn:inn+nn,j-nn:j+nn,k),1,[]);
                    y=reshape(S0(inn-nn:inn+nn,j-nn:j+nn,k),1,[]);
                    if sum(sum(x1~=0))>5        % It is not necessary to have a value in every position in the matrix of 5 to 5. 
                                                % When the number greater than 0 exceeds 5, a fitting can be carried out.
                        A=polyfit(x1,y,1);      % Slope and Offset: R1-S0
                        B=polyfit(x2,y,1);      % Slope and Offset: R2-S0
                        para1(inn,j,k)=A(1,1);  % Slope between R1 and S0    
                        para2(inn,j,k)=A(1,2);  % Offset between R1 and S0    
                        para3(inn,j,k)=B(1,1);  % Slope between R2 and S0    
                        para4(inn,j,k)=B(1,2);  % Offset between R2 and S0    
                        yeah=yeah+1;            % Count the number of efficient points
                    end
                 end
            end
        end



            MAX1 = mean(max(max(para1(:,:,:)))); % Calculate a series of maximum values and calculate their average values to avoid the interference of extreme values.
            MIN1 = min(para1(para1~=0));         % Because only the part of GM is calculated, and the other part is 0, it will cause unnecessary interference. 
                                                    %So leave out the parts that are zero.

            MAX2 = mean(max(max(para2(:,:,:))));
            MIN2 = min(para2(para2~=0));
            
            
            MAX3 = mean(max(max(para3(:,:,:))));
            MIN3 = min(para3(para3~=0));
             
            MAX4 = mean(max(max(para4(:,:,:))));
            MIN4 = min(para4(para4~=0));
           
                       
            para1_min(:,:,:) = para1(:,:,:) - MIN1;  % [A-MIN(A)]/MAX(A)   Narrow the scope of the data to [0,1]
            para1_min(para1(:,:,:)==0) = 0;          % The zero parts of the data should not be included in the calculation, so it should be restored.  
            para2_min(:,:,:) = para2(:,:,:) - MIN2;
            para2_min(para2(:,:,:)==0) = 0;
            para3_min(:,:,:) = para3(:,:,:) - MIN3;
            para3_min(para3(:,:,:)==0) = 0;
            para4_min(:,:,:) = para4(:,:,:) - MIN4;
            para4_min(para4(:,:,:)==0) = 0;
            
            para1_standard(:,:,:) = para1_min(:,:,:)./MAX1;  % [A-MIN(A)]/MAX(A)
            para2_standard(:,:,:) = para2_min(:,:,:)./MAX2;
            para3_standard(:,:,:) = para3_min(:,:,:)./MAX3;
            para4_standard(:,:,:) = para4_min(:,:,:)./MAX4;

%             para11(:,:,:) = para1_standard(:,:,:).*((2^10)-1);  %WM   
%             para22(:,:,:) = para2_standard(:,:,:).*((2^11)-1);
%             para33(:,:,:) = para3_standard(:,:,:).*((2^10)-1);
%             para44(:,:,:) = para4_standard(:,:,:).*((2^11)-1);
            
            para11(:,:,:) = para1_standard(:,:,:).*((2^10)-1);   %GM  % According to the difference of GM and WM, as well as the difference of Slope and Offset, they are enlarged to [0, 65535] respectively. 
                                                                      % Magnification standard: can clearly see the change of color in DICOM reader.
            para22(:,:,:) = para2_standard(:,:,:).*((2^12)-1);
            para33(:,:,:) = para3_standard(:,:,:).*((2^10)-1);
            para44(:,:,:) = para4_standard(:,:,:).*((2^12)-1);


            
            para_new_1(:,:,:) = uint16(para11(:,:,:));
            para_new_2(:,:,:) = uint16(para22(:,:,:));
            para_new_3(:,:,:) = uint16(para33(:,:,:));
            para_new_4(:,:,:) = uint16(para44(:,:,:));
        





            highflip_info_S1_1 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_1 = para_new_1(:,:,1);

         
            dicomwrite(para1_1, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_1.ima'], highflip_info_S1_1 );
            para2_1 = para_new_2(:,:,1);
            dicomwrite(para2_1, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_1.ima'] , highflip_info_S1_1 );
            para3_1 = para_new_3(:,:,1);
            dicomwrite(para3_1, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_1.ima'] , highflip_info_S1_1 );
            para4_1 = para_new_4(:,:,1);
            dicomwrite(para4_1, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_1.ima'] , highflip_info_S1_1 );
            
            
            
            
            highflip_info_S1_2 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_2 = para_new_1(:,:,2);
            dicomwrite(para1_2, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_2.ima' ], highflip_info_S1_2 );
            para2_2 = para_new_2(:,:,2);
            dicomwrite(para2_2, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_2.ima' ], highflip_info_S1_2 );
            para3_2 = para_new_3(:,:,2);
            dicomwrite(para3_2, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_2.ima' ], highflip_info_S1_2);
            para4_2 = para_new_4(:,:,2);
            dicomwrite(para4_2, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_2.ima' ], highflip_info_S1_2 );
            
            
            highflip_info_S1_3 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_3 = para_new_1(:,:,3);
            dicomwrite(para1_3, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_3.ima' ], highflip_info_S1_3 );
            para2_3 = para_new_2(:,:,3);
            dicomwrite(para2_3, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_3.ima' ], highflip_info_S1_3 );
            para3_3 = para_new_3(:,:,3);
            dicomwrite(para3_3, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_3.ima' ], highflip_info_S1_3);
            para4_3 = para_new_4(:,:,3);
            dicomwrite(para4_3, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_3.ima' ], highflip_info_S1_3 );
            
            
            
            highflip_info_S1_4 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_4 = para_new_1(:,:,4);
            dicomwrite(para1_4, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_4.ima' ], highflip_info_S1_4 );
            para2_4 = para_new_2(:,:,4);
            dicomwrite(para2_4, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_4.ima' ], highflip_info_S1_4 );
            para3_4 = para_new_3(:,:,4);
            dicomwrite(para3_4, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_4.ima' ], highflip_info_S1_4);
            para4_4 = para_new_4(:,:,4);
            dicomwrite(para4_4, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_4.ima' ], highflip_info_S1_4 );
            
            
            highflip_info_S1_5 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_5 = para_new_1(:,:,5);
            dicomwrite(para1_5, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_5.ima' ], highflip_info_S1_5 );
            para2_5 = para_new_2(:,:,5);
            dicomwrite(para2_5, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_5.ima' ], highflip_info_S1_5 );
            para3_5 = para_new_3(:,:,5);
            dicomwrite(para3_5, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_5.ima' ], highflip_info_S1_5);
            para4_5 = para_new_4(:,:,5);
            dicomwrite(para4_5, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_5.ima' ], highflip_info_S1_5 );
            

            highflip_info_S1_6 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_6 = para_new_1(:,:,6);
            dicomwrite(para1_6, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_6.ima' ], highflip_info_S1_6 );
            para2_6 = para_new_2(:,:,6);
            dicomwrite(para2_6, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_6.ima' ], highflip_info_S1_6 );
            para3_6 = para_new_3(:,:,6);
            dicomwrite(para3_6, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_6.ima' ], highflip_info_S1_6);
            para4_6 = para_new_4(:,:,6);
            dicomwrite(para4_6, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_6.ima' ], highflip_info_S1_6 );
            
            
            highflip_info_S1_7 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_7 = para_new_1(:,:,7);
            dicomwrite(para1_7, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_7.ima' ], highflip_info_S1_7 );
            para2_7 = para_new_2(:,:,7);
            dicomwrite(para2_7, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_7.ima' ], highflip_info_S1_7 );
            para3_7 = para_new_3(:,:,7);
            dicomwrite(para3_7, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_7.ima' ], highflip_info_S1_7);
            para4_7 = para_new_4(:,:,7);
            dicomwrite(para4_7, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_7.ima' ], highflip_info_S1_7 );
            
            highflip_info_S1_8 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_8 = para_new_1(:,:,8);
            dicomwrite(para1_8, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_8.ima' ], highflip_info_S1_8 );
            para2_8 = para_new_2(:,:,8);
            dicomwrite(para2_8, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_8.ima' ], highflip_info_S1_8 );
            para3_8 = para_new_3(:,:,8);
            dicomwrite(para3_8, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_8.ima' ], highflip_info_S1_8);
            para4_8 = para_new_4(:,:,8);
            dicomwrite(para4_8, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_8.ima' ], highflip_info_S1_8 );
            
            highflip_info_S1_9 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_9 = para_new_1(:,:,9);
            dicomwrite(para1_9, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_9.ima' ], highflip_info_S1_9 );
            para2_9 = para_new_2(:,:,9);
            dicomwrite(para2_9, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_9.ima' ], highflip_info_S1_9 );
            para3_9 = para_new_3(:,:,9);
            dicomwrite(para3_9, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_9.ima' ], highflip_info_S1_9);
            para4_9 = para_new_4(:,:,9);
            dicomwrite(para4_9, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_9.ima' ], highflip_info_S1_9 );

            
            highflip_info_S1_10 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_10 = para_new_1(:,:,10);
            dicomwrite(para1_10, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_10.ima' ], highflip_info_S1_10 );
            para2_10 = para_new_2(:,:,10);
            dicomwrite(para2_10, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_10.ima' ], highflip_info_S1_10 );
            para3_10 = para_new_3(:,:,10);
            dicomwrite(para3_10, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_10.ima' ], highflip_info_S1_10);
            para4_10 = para_new_4(:,:,10);
            dicomwrite(para4_10, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_10.ima' ], highflip_info_S1_10 );
            
            highflip_info_S1_11 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_11 = para_new_1(:,:,11);
            dicomwrite(para1_11, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_11.ima' ], highflip_info_S1_11 );
            para2_11 = para_new_2(:,:,11);
            dicomwrite(para2_11, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_11.ima' ], highflip_info_S1_11 );
            para3_11 = para_new_3(:,:,11);
            dicomwrite(para3_11, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_11.ima' ], highflip_info_S1_11);
            para4_11 = para_new_4(:,:,11);
            dicomwrite(para4_11, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_11.ima' ], highflip_info_S1_11 );

            
            highflip_info_S1_12 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
             para1_12 = para_new_1(:,:,12);
            dicomwrite(para1_12, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_12.ima' ], highflip_info_S1_12 );
            para2_12 = para_new_2(:,:,12);
            dicomwrite(para2_12, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_12.ima' ], highflip_info_S1_12 );
            para3_12 = para_new_3(:,:,12);
            dicomwrite(para3_12, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_12.ima' ], highflip_info_S1_12);
            para4_12 = para_new_4(:,:,12);
            dicomwrite(para4_12, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_12.ima' ], highflip_info_S1_12 );
            
            highflip_info_S1_13 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_13 = para_new_1(:,:,13);
            dicomwrite(para1_13, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_13.ima' ], highflip_info_S1_13 );
            para2_13 = para_new_2(:,:,13);
            dicomwrite(para2_13, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_13.ima' ], highflip_info_S1_13 );
            para3_13 = para_new_3(:,:,13);
            dicomwrite(para3_13, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_13.ima' ], highflip_info_S1_13);
            para4_13 = para_new_4(:,:,13);
            dicomwrite(para4_13, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_13.ima' ], highflip_info_S1_13 );
            
            highflip_info_S1_14 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_14 = para_new_1(:,:,14);
            dicomwrite(para1_14, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_14.ima' ], highflip_info_S1_14 );
            para2_14 = para_new_2(:,:,14);
            dicomwrite(para2_14, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_14.ima' ], highflip_info_S1_14 );
            para3_14 = para_new_3(:,:,14);
            dicomwrite(para3_14, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_14.ima' ], highflip_info_S1_14);
            para4_14 = para_new_4(:,:,14);
            dicomwrite(para4_14, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_14.ima' ], highflip_info_S1_14 );
            
            highflip_info_S1_15 = dicominfo(['C:\Users\...\T1Optimisation\COLOW000',num2str(i),'\gre_highflip_orig\File000',num2str(i),'.ima']);
            para1_15 = para_new_1(:,:,15);
            dicomwrite(para1_15, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R1_S0_Slope_15.ima' ], highflip_info_S1_15 );
            para2_15 = para_new_2(:,:,15);
            dicomwrite(para2_15, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R1_S0_Offset_15.ima' ], highflip_info_S1_15 );
            para3_15 = para_new_3(:,:,15);
            dicomwrite(para3_15, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Slope\COLOW000',num2str(i),'_GM_R2_S0_Slope_15.ima' ], highflip_info_S1_15);
            para4_15 = para_new_4(:,:,15);
            dicomwrite(para4_15, ['C:\Users\...\Schritt9 Dicom Datei\CO\COLOW000',num2str(i),'\GM\Offset\COLOW000',num2str(i),'_GM_R2_S0_Offset_15.ima' ], highflip_info_S1_15 );

%end