%load ED
%ED is the expression data: raws are genes, columns are conditions/time
%points
R=abs(corrcoef(transpose(ED)));
R=R-diag(diag(R));
i=0;
y=0;
rk=0;
h=0;
% the first for loop is to determine hard-threshold and
% corresponding R2 and slope values
for ht=0.70:0.01:0.99;
    h=h+1;
    HT=logical(R>=ht)*1;
    % HT is a binary adjacency matrix obtained by applying threshold 'ht'
    D=sum(HT,2);
    % D is degrees of HT
    ps=histogram(D,7);
    for i=1:7
        z(i)=ps.BinEdges(i)+ps.BinWidth/2;
    end
    linefit=fitlm(log10(z),log10(ps.BinCounts));
    rkht(h)=linefit.Rsquared.Ordinary;
    % rkht is the R2 value corresponding to 'ht'
    C=table2array(linefit.Coefficients);
    slopeht(h)=C(2,1);    
    % slopeht is the slope of line-fit corresponding to 'ht'
end
% The best tau will be selected
tau=0;
for i=1:length(rkht)
    if rkht(i)>= 0.85
        if slopeht(i)<0
            tau=0.69+(i*0.01);
            break
        end
    end
end
% writing co-expression matrix created by tau on .csv file
if tau==0
        display('no tau satisfying scale-freeness could be detected')
else
    HT=logical(R>=tau)*1;
    save HT HT
    csvwrite('COVID1CSV-tau.csv', HT);
    % -----------------------------------------------------------------------
    % the second part is to determine soft-thresholding parameters for ASF
    % and corresponding R2 and slope values
    x=0;
    m=tau;
    for A=7:50;
        x=x+1;
        ET=(1+exp(-A*(R-m)));
        ET=ET.^-1;
        % ET is transformed correlations matrix by asymmetric sigmoid function
        D=sum(ET,2);
        % D is vector of degrees of ET
        ps=histogram(D,7);
        for i=1:7
            z(i)=ps.BinEdges(i)+ps.BinWidth/2;
        end
        linefit=fitlm(log10(z),log10(ps.BinCounts));
        rk(x,1)=linefit.Rsquared.Ordinary;
        % rk is the R2 value corresponding to 'A' & 'm'
        C=table2array(linefit.Coefficients);
        slope(x,1)=C(2,1);
        % slope is the slope of line-fit corresponding to 'A' & 'm'
    end
    % The best A (alpha) will be selected
    A=0;
    for i=1:length(rk)
        if rk(i,1)>= 0.85
            if slope(i,1)<0
                A=6+i;
                ET=(1+exp(-A*(R-m)));
                ET=ET.^-1;
                break
            end
        end
    end
    % writing co-expression matrix created by ASF on .csv file
    if A==0
        A=50;
        ET=(1+exp(-A*(R-m)));
        ET=ET.^-1;
        save ET ET
        csvwrite('COVID1CSV-ASF.csv', ET);    
    else
        save ET ET
        csvwrite('COVID1CSV-ASF.csv', ET);
    end
end
