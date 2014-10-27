function [] = hingemdfiperturb(list)
    %close all; clear
    hlist = list 
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
            j
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
    fdfi = sum(nrml(:,hlist),2);
    %mdfi = sum(nrml(hlist,:),1) %hlist pulls out the hinges and does the sum over those 

    
    outfile = fopen('fdfi-Avg.dat','w');
    for i=1:size(fmdfi,2);
        fprintf(outfile,'%f\n',mdfi(i));
    end
    fclose(outfile);
    
    outfile = fopen('fdfiindex.debug','w');
    for i=hlist
       fprintf(outfile,'%d\n',i); 
    end
    fclose(outfile);
end
