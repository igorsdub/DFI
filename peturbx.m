close all; clear
INVBKBT=importdata('pinv_svd.debug');
if size(INVBKBT,2)==1
    INVBKBT=reshape(INVBKBT,sqrt(size(INVBKBT,1)),[]);
end
resnum=size(INVBKBT,1)/3;
resnum2=resnum;
directions = [1 0 0]; 
perbMat = zeros(resnum,resnum);
directions = [1 0  0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1]; 
numdirect = size(directions,1);
for k=1:numdirect 
direct = directions(k,:)/norm(directions(k,:))

for j=1:resnum
        delforce=zeros(resnum*3,1); %Force column 
        delforce(3*j-2:3*j)=direct; %Foce in the force colmn 
        delXperbVec=INVBKBT*delforce; %displacement 
        delXperbMat=reshape(delXperbVec,3,[])'; %put in a matrix 
        delRperbVec=sqrt(sum(delXperbMat.^2,2)); %get the magnitude 
        delRperbVecnormed=delRperbVec; %normalize, may be unecessay because giong to normalize in the end
        for i=1:resnum
           perbMat(i,j) = perbMat(i,j) + delRperbVec(i);
        end
end

end
nrml=perbMat.*(1/numdirect); %average by the number of directions 
nrml = nrml./(sum(nrml(:)));
dfi = sum(nrml,2)
mdfi = sum(nrml,1)

outfile = fopen('S1-Avg.dat','w');
for i=1:size(dfi,1)
   fprintf(outfile,'%f\n',dfi(i)); 
end
fclose(outfile)

outfile = fopen('S2-Avg.dat','w');
for i=1:size(mdfi,2)
    fprintf(outfile,'%f\n',mdfi(i));
end
fclose(outfile)