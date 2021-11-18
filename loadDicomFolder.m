function [mFrame,iNoFiles,infoFirst,dcmInfoList] = loadDicomFolder(Folder, doSort, maxFiles)

%***************************** start loading images from files ****************************** 
echo off
bVar=0;
tmpCnt=1;
cPref=Folder;                                             % full path to folder content
files=dir(cPref);                                        % list of files *.*
iNoFiles=max(size(files));                                % number of loaded images
for i=1:iNoFiles                                          % loop for checking of file format
    tmpName=files(i).name;
    if(tmpName(1)=='.') 
           bVar=bVar+1;    
     end
end                                                       % loop for checking of file format
cRoot=files(bVar+1).name;         
cFull=strcat(cPref,cRoot);
infoFirst = dicominfo(cFull);
%infoFirst.RepetitionTime;
%infoFirst.InversionTime
%infoFirst.EchoTime
%infoFirst.FlipAngle
iImgResolR=infoFirst.Rows;                                     % img resolution of the first file in list
iImgResolC=infoFirst.Columns; 
iNoFiles=iNoFiles-bVar;
if(maxFiles==0)
    maxFiles=iNoFiles;
end
mFrame=zeros(iImgResolR,iImgResolC,maxFiles);              % multi frame declaration 
iNoInstance=zeros(maxFiles,1);                            % offset number of sorted image files
dcmInfoList = cell(maxFiles,1);
for iCountFiles=1:iNoFiles                                % loop for loading images, size is iNo of files
     cRoot=files(iCountFiles+bVar).name;         
     cFull=strcat(cPref,cRoot);   
     info = dicominfo(cFull);
     if (info.Rows ~= iImgResolR || info.Columns ~= iImgResolC)  
          iNoFiles=iCountFiles-1;
          break;                                          % case of changed resolution breaks loop
     end
     if(info.InstanceNumber>maxFiles)
         continue;
     end;
     if( doSort == 1 )
         iNoInstance(tmpCnt)=info.InstanceNumber;
     end
     [tmp1, map] = dicomread(cFull);
     tmp1=double(tmp1);
     mFrame(:,:,tmpCnt)=tmp1;   
     dcmInfoList{iCountFiles}=info;
     tmpCnt=tmpCnt+1;
end                                                       % loop for loading images
if( doSort == 1 )
    [tmp,ISort]=sort(iNoInstance);                      % sort image list
    mTmpFrame=mFrame;
    for( i=1:maxFiles )
        mFrame(:,:,i)=mTmpFrame(:,:,ISort(i));
    end
end
return
