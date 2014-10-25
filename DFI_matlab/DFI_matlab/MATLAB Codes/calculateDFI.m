close all; clear
INVBKBT=importdata('1YMZ_INVBKBT.txt');
if size(INVBKBT,2)==1
    INVBKBT=reshape(INVBKBT,sqrt(size(INVBKBT,1)),[]);
end
resnum=size(INVBKBT,1)/3;
resnum2=resnum;
npb=1;
resperturb=[1];
directions=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1];
ndirect=size(directions,1);
delRperbDistMat=repmat(0, [resnum resnum ndirect]);
ratio=zeros(resnum2,ndirect);
bindcolumnsum=zeros(resnum2,ndirect);
allcolumnsum=zeros(resnum2,ndirect);
for k=1:ndirect       %k is the dimension of direction
    direct=directions(k,:)/norm(directions(k,:))        
    for j=1:resnum          % perturb j
        delforce=zeros(resnum*3,1);
        delforce(3*j-2:3*j)=direct;
        delXperbVec=INVBKBT*delforce;
        delXperbMat=reshape(delXperbVec,3,[])';
        delRperbVec=sqrt(sum(delXperbMat.^2,2));
        delRperbVecnormed=delRperbVec/sum(delRperbVec);
%        delRperbVecnormed=delRperbVec;
        for i=1:resnum           % response i
            if abs(j-i)>3
                delRperbDistMat(i,j,k)=delRperbVecnormed(i);
            end
        end
    end
    allcolumnsum(:,k)=sum(delRperbDistMat(1:resnum2,1:resnum2,k),2);
    bindcolumnsum(:,k)=sum(delRperbDistMat(1:resnum2,resperturb(1:npb),k),2);
    ratio(:,k)=(bindcolumnsum(:,k)/npb)./(allcolumnsum(:,k)/(resnum2-1));
end
delRperbDistMat
dRmat=zeros(resnum2, resnum2);
for i=1:resnum
    for j=resnum
        dRmat(i,j)=sum(delRperbDistMat(i,j,:));
    end
end
dRmat
dfi=sum(dRmat,1)/sum(dRmat(:))


%% addtional code of ARR calculation added for Banu, need double check.
% meanratio=mean(ratio,2);
% %arr=(tiedrank(meanratio)-1)/(length(meanratio)-1);
% arr=meanratio;
