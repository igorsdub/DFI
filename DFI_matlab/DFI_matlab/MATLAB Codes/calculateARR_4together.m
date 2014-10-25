clear; 
figure
INVBKBT=importdata('/Users/Taisong/Research/Rebekka/allostery/1atmCR8/tetramer_struct2_MD_LEA/4DXNmwcovarmat.dat');
%INVBKBT=importdata('/Users/Taisong/Research/Rebekka/allostery/1atmCR8/tetramer_struct2/4DXM_INVBKBT.txt');
if size(INVBKBT,2)==1
    INVBKBT=reshape(INVBKBT,sqrt(size(INVBKBT,1)),[]);
end
resnum=size(INVBKBT,1)/3;
resnum0=resnum/4;
npb=40;
npb0=npb/4;
%resperturb=[53,60,65,95,96,107,145,148,185,207,260,267,272,302,303,314,352,355,392,414,467,474,479,509,510,521,559,562,599,621,674,681,686,716,717,728,766,769,806,828];
%resperturb0=[53,60,65,95,96,107,145,148,185,207];
ndirect=10;
ncombine=100;
directions=zeros(ndirect,npb0*3);      % monomer
combinations=zeros(ncombine,npb0);    % monomer
for j=1:ncombine     % j is the dimension of residue combination    
    if j==1
        combine=[53,60,65,95,96,107,145,148,185,207];
        sampleset=1:207;
        sampleset(combine)=[];
    else
%        combine=randsample(resnum0,npb0)';  %include mutation sites
        z=randperm(numel(sampleset));
        combine=z(1:npb0)';
    end
    combinations(j,:)=combine;
end
dlmwrite('mycombinations.txt',combinations,'delimiter','\t');

delRperbDistMat=repmat(0, [resnum ncombine ndirect]);
ratio=zeros(resnum,ndirect);
bindcolumnsum=zeros(resnum,ndirect);
allcolumnsum=zeros(resnum,ndirect);
for k=1:ndirect       %k is the dimension of direction
    direct_tmp=rand(1,npb0*3);
%    direct_tmp=horzcat([1,0,0], repmat([0,0,0],1,9));
    direct=direct_tmp/norm(direct_tmp);
    directions(k,:)=direct;
    for j=1:ncombine     % j is the dimension of residue combination    
        combine=combinations(j,:);
        combine3X_tmp=vertcat(combine*3-2,combine*3-1,combine*3);
        combine3X=reshape(combine3X_tmp,1,[]);
        delforce0=zeros(resnum0*3,1);
        delforce0(combine3X)=direct;
        delforce=repmat(delforce0,4,1);
        delXperbVec=INVBKBT*delforce;
        delXperbMat=reshape(delXperbVec,3,[])';
        delRperbVec=sqrt(sum(delXperbMat.^2,2));
        delRperbVecnormed=delRperbVec/sum(delRperbVec);
    %    for i=1:resnum           % response i
    %        if abs(j-i)>3
    %            delRperbDistMat(i,j,k)=delRperbVecnormed(i);
    %        end
    %    end
        delRperbDistMat(:,j,k)=delRperbVecnormed;
    end
    allcolumnsum(:,k)=sum(delRperbDistMat(:,:,k),2);
    bindcolumnsum(:,k)=delRperbDistMat(:,1,k);
%    ratio(:,k)=bindcolumnsum(:,k)./(allcolumnsum(:,k)/ncombine);
    ratio(:,k)=bindcolumnsum(:,k);
end
dlmwrite('mydirections.txt',directions,'delimiter','\t');


dfi_4dxm=importdata('/Users/Taisong/Dropbox/thesis/mythesis/GFP/colordfi/dfiRank_4DXM');
dfi_4dxn=importdata('/Users/Taisong/Dropbox/thesis/mythesis/GFP/colordfi/dfiRank_4DXN');
diff=dfi_4dxn(:,2)-dfi_4dxm(:,2);
%diff=importdata('/Users/Taisong/Dropbox/thesis/mythesis/GFP/colordfi/diffdfi');

meanratio=mean(ratio,2);
%arr=(tiedrank(meanratio)-1)/(length(meanratio)-1);
arr=meanratio;
% average over 4 units
arr_tmp=reshape(arr,[],4);
%arr0=mean(arr_tmp,2);
arr0=max(arr_tmp,[],2);
resid=dfi_4dxm(:,1);

redres=[67,68,69,70,71,142,187,188,189,190,191,192,193,213,214,215];
blueres=[19,20,21,22,50,97,98,125,127,128,129,165,167];
redblueres=horzcat(redres,blueres);

redidx=[];
blueidx=[];
for i=1:length(redres)
    redidx(end+1)=find(resid==redres(i));
end
for i=1:length(blueres)
    blueidx(end+1)=find(resid==blueres(i));
end
redblueidx=horzcat(redidx,blueidx);

phandle = plot(arr0, diff,'ok');
set(phandle, 'MarkerSize', 5)
hold on

phandle = plot(arr0(redidx),diff(redidx),'or');
set(phandle, 'MarkerSize', 5, 'MarkerFaceColor', 'r')
hold on

phandle = plot(arr0(blueidx),diff(blueidx),'ok');
set(phandle, 'MarkerSize', 5, 'MarkerFaceColor', 'b')
hold on

temp = polyfit(arr0,diff,1);
fit = temp(1)*arr0+temp(2);
[R,p]=corrcoef(arr0, diff);
coor = sprintf('%.2f', R(1,2));
pvalue = sprintf('%.2e', p(1,2));

plot(arr0, fit, '-k');
hold on

t1=text(0.05, 0.7, strcat('R= ', coor), 'fontsize',16);
t2=text(0.05, 0.6, strcat('$p$= ', pvalue), 'fontsize',16);
set(t1, 'interpreter','latex','FontSize', 16)
set(t2, 'interpreter','latex','FontSize', 16)

%axis([0 40 -5 15])
set(gca,'FontSize',14)
xlabel('$ARR$','interpreter','latex', 'fontsize',20)
ylabel('$\Delta\%dfi$', 'interpreter','latex', 'fontsize',20)

%saveas(gcf,'/Users/Taisong/Dropbox/thesis/mythesis/GFP/dfiVSarr.eps', 'epsc');

edges=-0.7:0.3:0.9;
[n, bins]=histc(diff, edges);
arr_aver=zeros(size(n));
for i=1:length(n)
    idx=find(bins==i);
    arr_aver(i)=sum(arr0(idx))/length(idx);
end
figure
boxplot(arr0, edges(bins))
set(gca,'FontSize',14)
ylabel('$ARR$','interpreter','latex', 'fontsize',20)
xlabel('$\Delta\%dfi$', 'interpreter','latex', 'fontsize',20)
