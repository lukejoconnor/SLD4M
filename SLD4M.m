function [ maest,maerr,h2est,h2err,maenrichest,maenricherr,h2enrichest,h2enricherr,...
    ma_jk,h2_jk,kurtexpest,kurtexperr,propkurtest,propkurterr ] ...
    = SLD4M( chisq,ell2,ell4,idxldsc,annot,weights_ldsc,weights_ld4m,no_blocks,ref_annot,...
    report_annot_mat)
%SLD4M (stratified LD 4th moments regression).
%   INPUT ARGUMENTS:chisq,  Mx1 matrix of chi^2 statistics; ell2, MxP matrix of
%   LD scores (LD second moments), where P is #annotations; ell4, MxP
%   vector of LD 4th moments; idxldsc, Mx1 vector of indices corresponding
%   to the reference SNPs that will be used in the regression; annot,
%   MtotxP matrix of annotation values, where Mtot is the total number of
%   reference SNPs (may be greater than M); weights_ldsc, Mx1 vector of
%   weights for LDSC; weights_ld4m, Mx1 vector of weights for LD4M;
%   no_blocks, number of jackknife blocks; ref_annot <= P', index of the
%   annotation that will be used in denominator to compute enrichments;
%   report_annot_mat, a PxP' matrix whose columns correspond to linear
%   combinations of input annotations to output estimates for. For example,
%   if P=2 with annotations for common and LF SNPs, then set
%   report_annot_mat to [[1; 1], [1; 0], [0; 1]] and SLD4M will output Ma
%   estimates for 3 annotations: all SNPs, common SNPs, LF SNPs.
%   OUTPUT ARGUMENTS: maest, Estimated Ma for each category;
%   maerr, Ma standard error; h2est, heritability estimate; h2err,
%   heritability standard error; maenrichest, estimated Ma enrichment;
%   maenricherr, Ma enrichment standard error; h2enrichest, estimated
%   heritability enrichment; h2enricherr, heritability enrichment standard
%   error; ma_jk, leave-one-out jackknife estimates of Ma; h2_jk,
%   leave-one-out jackknife estimates of heritability; kurtexpest, kurtosis
%   explained by functional annotations; kurtexperr, standard error;
%   propkurtest, proportion of kurtosis explained by functional annotations
%   (log scale); propkurterr, standard error.


if ~exist('ref_annot')
    ref_annot=1;
end
if ~exist('report_annot_mat')
    report_annot=annot;
    report_annot_mat=speye(size(annot,2));
else
    report_annot=annot*report_annot_mat;
    report_l2=ell2*report_annot_mat;
end

%info_mat=diag(sparse(info));% Info score matrix
WW2_mat=diag(sparse(weights_ldsc));% 2nd-moment weights matrix
WW4_mat=diag(sparse(weights_ld4m));% 4th-moment weights matrix

mm=length(idxldsc);
[mm2,pp]=size(ell2);
pp=pp+1;
pp2=size(report_annot,2);
ell2=[ones(mm2,1) ell2];
ell4=[ones(mm2,1) ell4];

blocksize=floor(mm/no_blocks);
l2var=zeros(pp,pp,no_blocks);
l4var=l2var;
l2l4a2cov=l2var;
l2l4cov=l2var;
annotvar=zeros(pp-1,pp2,no_blocks);
l2annot=annotvar;
l2a2cov=zeros(pp,no_blocks);
l4a4cov=l2a2cov;
l4a2cov=l2a2cov;
l2mean=l2a2cov;
l4mean=l2a2cov;
a2annot_ref=zeros(no_blocks,pp-1);
annotmean_ref=zeros(no_blocks,pp-1);

% Compute various moments for each jackknife block
for jk=1:no_blocks
    ind=(jk-1)*blocksize+1:jk*blocksize;
    ind2=idxldsc(ind);
    l2var(:,:,jk)=ell2(ind2,:)'*WW2_mat(ind,ind)*ell2(ind2,:)/sum(weights_ldsc(ind));
    l4var(:,:,jk)=ell4(ind2,:)'*WW4_mat(ind,ind)*ell4(ind2,:)/sum(weights_ld4m(ind));
    l2l4a2cov(:,:,jk)=(ell2(ind2,:).*chisq(ind))'*WW4_mat(ind,ind)*ell4(ind2,:)/sum(weights_ld4m(ind));
    l2l4cov(:,:,jk)=(ell2(ind2,:))'*WW4_mat(ind,ind)*ell4(ind2,:)/sum(weights_ld4m(ind));
    l4a4cov(:,jk)=ell4(ind2,:)'*WW4_mat(ind,ind)*(chisq(ind).^2)/sum(weights_ld4m(ind));
    l2a2cov(:,jk)=ell2(ind2,:)'*WW2_mat(ind,ind)*chisq(ind)/sum(weights_ldsc(ind));
    l4a2cov(:,jk)=ell4(ind2,:)'*WW4_mat(ind,ind)*chisq(ind)/sum(weights_ld4m(ind));
    l4mean(:,jk)=mean(WW4_mat(ind,ind)*ell4(ind2,:))/mean(weights_ld4m(ind));
    l2mean(:,jk)=mean(WW2_mat(ind,ind)*ell2(ind2,:))/mean(weights_ldsc(ind));
    annotvar(:,:,jk)=annot(ind2,:)'*report_annot(ind2,:)./sum(report_annot(ind2,:));
    l2annot(:,:,jk)=ell2(ind2,2:end)'*report_annot(ind2,:)./sum(report_annot(ind2,:));
    a2annot_ref(jk,:)=mean(annot(ind2,:).*chisq(ind).*report_annot(ind2,ref_annot))/mean(report_annot(ind2,ref_annot));
    annotmean_ref(jk,:)=mean(annot(ind2,:).*report_annot(ind2,ref_annot))/mean(report_annot(ind2,ref_annot));
    
end

% Combine jackknife estimates to obtain leave-one-out regression
% coefficients
for jk=1:no_blocks
    ind=[1:jk-1,jk+1:no_blocks];
    temp=mean(l2var(:,:,ind),3)\mean(l2a2cov(:,ind),2);% l2var^-1 * l2a2
    tau(:,jk)=temp(2:end);
    intercept(jk)=temp(1);
    
    total_cov(jk,:)=mean(l4a4cov(2:end,ind)-...
        6*intercept(jk)*l4a2cov(2:end,ind)+...
        3*intercept(jk)^2*l4mean(2:end,ind),2)'-...
        3*tau(:,jk)'*mean(l2l4a2cov(2:end,2:end,ind),3)+...
        3*intercept(jk)*tau(:,jk)'*mean(l2l4cov(2:end,2:end,ind),3);
    
    total_var=mean(l4var(2:end,2:end,ind),3);
    
    kurtcoef(:,jk)=total_var\total_cov(jk,:)';
    
    kurttot(:,jk)=kurtcoef(:,jk)'*mean(annotvar(:,:,ind),3);
    beta2tot(:,jk)=tau(:,jk)'*mean(annotvar(:,:,ind),3);
    alpha2tot(:,jk)=tau(:,jk)'*mean(l2annot(:,:,ind),3);
    
end

% Point estimates and standard errors
h2_jk=beta2tot.*sum(report_annot)';
h2est=mean(h2_jk,2);
h2err=std(h2_jk')'*sqrt(no_blocks+1);

kurttot=(kurttot+3*beta2tot.*alpha2tot)./beta2tot.^2;

kurttot_normalized=kurttot.*beta2tot./alpha2tot;

kurtest=mean(kurttot_normalized,2);
kurterr=std(kurttot_normalized')'*sqrt(no_blocks+1);

ma_jk=3*sum(report_annot)'./kurttot;
maest=mean(ma_jk,2);
maerr=std(ma_jk')'*sqrt(no_blocks+1);
malogerr=std(log(ma_jk'))*sqrt(no_blocks+1);

h2genrich=h2_jk./h2_jk(ref_annot,:).*mean(report_annot(:,ref_annot))./mean(report_annot)';
maenrich=ma_jk./ma_jk(ref_annot,:)*mean(report_annot(:,ref_annot))./mean(report_annot)';

h2enrichest=mean(h2genrich,2);
h2enricherr=std(h2genrich')'*sqrt(no_blocks+1);

maenrichest=mean(maenrich,2);
maenricherr=std(maenrich')'*sqrt(no_blocks+1);

piest=mean(3./kurttot_normalized');
pierr=std(3./kurttot_normalized')*sqrt(no_blocks+1);

% Estimate kurtosis explained by functional annotations
if nargout>10
    kurt_exp=zeros(no_blocks,1);
    for jk=1:no_blocks
        ind=[1:jk-1,jk+1:no_blocks];
        kurt_exp(jk)=3*sum(mean(a2annot_ref(ind,:) - intercept(jk)*annotmean_ref(ind,:)).*tau(:,jk)')...
            ./beta2tot(ref_annot,jk)./alpha2tot(ref_annot,jk);
        ind=(jk-1)*blocksize+1:jk*blocksize;        
    end
    
    kurtexpest=mean(kurt_exp);
    kurtexperr=std(kurt_exp)*sqrt(no_blocks+1);
    
    prop_kurtexp_jk=log(kurt_exp/3)./...
        log(kurttot_normalized(ref_annot,:)/3);
    propkurtest=mean(prop_kurtexp_jk);
    propkurterr=std(prop_kurtexp_jk)*sqrt(no_blocks+1);
    
end




