% -- prior knowledge from interactome --
load('./PriorKnowledge/iRefIndex_adj_matrix.mat');
len_gene_net = length(GeneSymbol_net);
% -- normalization --
D_half_inv = sparse(diag(sum(network_adj_matrix,1).^-(0.5)));
Laplacian_mat = speye(len_gene_net) - D_half_inv*network_adj_matrix*D_half_inv;
clear D_half_inv Gene2Path_normalized len_gene_net network_adj_matrix

% -- prior knowledge from pathways --
load('./PriorKnowledge/Pathway_map_bipartite.mat');
Ind_vec = GetInd_in_Dict(Gene_list,GeneSymbol_net);
Gene2Path_cur = sparse(length(GeneSymbol_net),size(Gene2Path_map,2));
Gene2Path_cur(Ind_vec(Ind_vec~=0),:) = Gene2Path_map(Ind_vec~=0,:);
% -- normalization --
vec_sum_col = sum(Gene2Path_cur,1);
vec_sum_col(vec_sum_col==0) = 1;
Gene2Path_normalized = Gene2Path_cur*diag(vec_sum_col.^-1);
clear Ind_vec Gene_list Gene2Path_cur Gene2Path_map vec_sum_col

mkdir('./Output_data');
% --- configure ---
LLP = 0.15;
MC = 5;
Gene2Path.lambda_C = 0.01;
Gene2Path.lambda_V = 0.01;
Net_Conf.lambda_L = 0.1;

Gene2Path.Gene2Path = Gene2Path_normalized;
Net_Conf.Lap_mat = Laplacian_mat;
clear Laplacian_mat Gene2Path_normalized

CancerName = {'BRCA','GBM','THCA'};
for i_cancer = 1
    % -- load data --
    load(['./Input_data/' CancerName{i_cancer} '.mat']);

    [S_sample_indicator,G_gene_score,V_pathway_score] = ...
        DGPathinter(mutation_mat,Net_Conf,Gene2Path,LLP,MC);

    [~,ind_gene] = sort(max(G_gene_score,[],2),'descend');
    PotentialDriverGenes = GeneSymbol_net(ind_gene(1:100));
    [~,ind_path] = sort(max(V_pathway_score,[],2),'descend');
    TopScoredPathways = Pathway_list(ind_path(1:30));

    save(['./Output_data/Output_'  CancerName{i_cancer} '.mat'],...
        'S_sample_indicator','G_gene_score','V_pathway_score',...
        'PotentialDriverGenes','TopScoredPathways');
end