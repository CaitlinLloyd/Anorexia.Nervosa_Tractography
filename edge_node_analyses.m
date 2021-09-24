cd('/DataLocation')
%NBS
%UI structure corresponding to the example data provided:
UI.method.ui='Run NBS'; 
UI.test.ui='t-test';
UI.size.ui='Extent';
UI.thresh.ui='3.1';
UI.perms.ui='10000';
UI.alpha.ui='0.05';
UI.contrast.ui='[-1,1,0,0,0]'; 
UI.design.ui='design_matrix.txt';
UI.exchange.ui=''; 
UI.matrices.ui='connectivity_matrices_combined.mat';
UI.node_coor.ui='atlasCORD.txt';                         
UI.node_label.ui='node_names.txt';
%
%   Remarks:
%       This function can be used as a command line version of NBS: 
%           1. Specify inputs in structure UI 
%           2. Run NBSrun(UI)
%           3. At completion, results stored in structure nbs
%              Type 'global nbs' before attempting to access the structure.
%


NBSrun(UI)


%node metrics
load('connectivity_matrices_combined.mat')
sublist=load('sublist.txt')
fname1=('node metrics.txt')
nodes=readtable('node_names.txt', 'ReadVariableNames',false)
covs=readtable('covariates.txt','ReadVariableNames',false)
group=covs(:,1)

for sub = (1:length(sublist))
    sub
    clear EF
    clear clust
    clear strength
    clear adj
    adj = x(:,:,sub);
    EF = efficiency_wei(adj,2);
    w=weight_conversion(adj, "normalize");
    clust=clustering_coef_wu(w); 
    trans=transitivity_wu(w);
    gl_ef=efficiency_wei(adj);
    strength = strengths_und(adj); 
    subname = sublist(sub);
    type=group(sub,1);
    type=table2array(type);
        if type == 1;
            diag="AN";
            else diag = "HC";
        end
    fid=fopen(fname1,'a');
for l = 1:84;
    if l==1 && sub==1;
    fprintf(fid,'\nsubname\tdiag\tEF\tstrength\tBET');   
    end
node=nodes(l,1);
node=table2array(node);
fprintf(fid,'\n%f\t%s\t%s\t%f\t%f\t%f', subname,diag, node{1}, EF(l), strength(l), clust(l));
end
glob="global_trans"
fprintf(fid,'\n%f\t%s\t%s\t%f', subname,diag, glob, trans);
glob="global_ef"
fprintf(fid,'\n%f\t%s\t%s\t%f', subname,diag, glob, gl_ef);
end

