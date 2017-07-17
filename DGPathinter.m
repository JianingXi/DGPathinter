function [S_vec,G_vec,V_vec,rank_k] = DGPathinter(X,NetConf,PathConf,LLP,MC)

[TotalNum, TotalGene] = size(X);
TotalPath = size(PathConf.Gene2Path,2);

DefaultLLP = 0.15;
DefaultMC = 5;

if ~exist('CompLeastProportion','var')
    LLP = DefaultLLP;
elseif LLP <= 0 || LLP >= 1
    disp('CompLeastProportion: out of domain. Use default value instead.');
    LLP = DefaultLLP;
end

if ~exist('maxComponent','var')
    MC = DefaultMC;
elseif MC <= 0
    disp('Max Component number: out of domain. Use default value instead.');
    MC = DefaultMC;
end

leastGroupNum = max(ceil(LLP*TotalNum),10);

S_vec = zeros(TotalNum,DefaultMC);
G_vec = zeros(TotalGene,DefaultMC);
V_vec = zeros(TotalPath,DefaultMC);

selected_idx_sample = zeros(TotalNum,1);

Remain_sample_num = TotalNum;

rank_k = 0;
while sum(selected_idx_sample) < TotalNum && ...
        rank_k < (MC-1) && Remain_sample_num >= leastGroupNum
    
    [s_vec_1,g_vec_1,v_vec_1] = RankOneLayerFactoring(X,NetConf,PathConf,leastGroupNum);
    rank_k = rank_k + 1;
    
    cur_idx_s = (s_vec_1~=0);
    S_vec_cur = zeros(TotalNum,1);
    S_vec_cur(selected_idx_sample==0)=s_vec_1;
    selected_idx_sample(selected_idx_sample==0) = cur_idx_s;
    
    S_vec(:,rank_k) = S_vec_cur;
    G_vec(:,rank_k)  = g_vec_1;
    V_vec(:,rank_k)  = v_vec_1;
    X = X(~cur_idx_s,:);
    
    Remain_sample_num = sum(~selected_idx_sample);
    if Remain_sample_num == 0
        return;
    end
end

if Remain_sample_num > 0
    [s_vec_1,g_vec_1,v_vec_1] = RemainDecomp(X,NetConf,PathConf);
    
    rank_k = rank_k + 1;
    
    S_vec_cur = zeros(TotalNum,1);
    S_vec_cur(selected_idx_sample==0) = s_vec_1;

    S_vec(:,rank_k) = S_vec_cur;
    G_vec(:,rank_k)  = g_vec_1;
    V_vec(:,rank_k) = v_vec_1;
end

ind_sample = (sum(abs(S_vec),1)>0);
ind_gene = (sum(abs(G_vec),1)>0);
S_vec = S_vec(:,ind_sample);
G_vec = G_vec(:,ind_gene);

end


function [s_vec,g_vec,v_vec,iter] = RankOneLayerFactoring(X,NetConf,PathConf,leastGroupNum)

eps_denom = 10^-5;
[n_samp,p_gene] = size(X);
TotalPath = size(PathConf.Gene2Path,2);

if sum(sum(abs(X))) == 0
    s_vec = ones(n_samp,1);
    g_vec = zeros(p_gene,1);
    iter = 0;
    return;
end

try
    [s_vec_0,~,~] = svds(X,1);
    if sum(s_vec_0<-(10^-4)) > sum(s_vec_0>(10^-4))
        % ensure most element of u0 is positive
        s_vec_0 = -s_vec_0;
    end
    s_vec_0(s_vec_0<=0)=0;
    s_vec_0 = s_vec_0/max(s_vec_0);
    g_vec_0 = X'*s_vec_0;
catch
    s_vec_0 = ones(n_samp,1);
    g_vec_0 = sum(X,1)';
end
v_vec_0 = zeros(TotalPath,1);


Sample_Number_in = 0;

relative_merr = 0.005;
niter = 50;

relative_s_vec_d = 2;
relative_g_vec_d = 2;
relative_v_vec_d = 2;
iter = 0;
while (relative_s_vec_d > relative_merr && relative_g_vec_d > relative_merr &&...
        relative_v_vec_d > relative_merr || Sample_Number_in < leastGroupNum) ...
        && iter < niter
        
    iter = iter+1;
    % ---- V_vec updating ----
    v_vec_1 = v_vec_Est(g_vec_0,PathConf);
    
    % ---- G_vec updating ----
    g_vec_1 = g_vec_Est(X,s_vec_0,v_vec_1,NetConf,PathConf);
    
    % ---- S_vec updating ----
    s_vec_1 = S_vec_binary(X,g_vec_1,leastGroupNum);
    Sample_Number_in = sum(s_vec_1~=0);
    
    % ---- Residual ----
    relative_s_vec_d = sum((s_vec_0-s_vec_1).^2)/(sum(s_vec_0.^2) + eps_denom);
    relative_g_vec_d = sum((g_vec_0-g_vec_1).^2)/(sum(g_vec_0.^2) + eps_denom);
    relative_v_vec_d = sum((v_vec_0-v_vec_1).^2)/(sum(v_vec_0.^2) + eps_denom);
    
    s_vec_0 = s_vec_1;
    g_vec_0 = g_vec_1;
    v_vec_0 = v_vec_1;
end

s_vec = s_vec_1;
g_vec = g_vec_1;
v_vec = v_vec_1;

end



function [s_vec,g_vec,v_vec,iter] = RemainDecomp(X,NetConf,PathConf)

eps_denom = 10^-5;
[n_samp,p_gene] = size(X);
s_vec = ones(n_samp,1);

Num_m = size(PathConf.Gene2Path,2);

if sum(sum(abs(X))) == 0
    g_vec = zeros(p_gene,1);
    iter = 0;
    return;
end
g_vec_0 = sum(X,1)';
v_vec_0 = zeros(Num_m,1);

relative_merr = 0.005;
niter = 50;


relative_g_vec_d = 2;
relative_v_vec_d = 2;

iter = 0;
while (relative_g_vec_d > relative_merr && relative_v_vec_d > relative_merr) ...
        && iter < niter
    iter = iter+1;
    
    % ---- V_vec updating ----
    v_vec_1 = v_vec_Est(g_vec_0,PathConf);
    
    % ---- G_vec updating ----
    g_vec_1 = g_vec_Est(X,s_vec,v_vec_1,NetConf,PathConf);
    
    % ---- Residual ----
    relative_g_vec_d = sum((g_vec_0-g_vec_1).^2)/(sum(g_vec_0.^2) + eps_denom);
    relative_v_vec_d = sum((v_vec_0-v_vec_1).^2)/(sum(v_vec_0.^2) + eps_denom);
    
	g_vec_0 = g_vec_1;
    v_vec_0 = v_vec_1;
end

g_vec = g_vec_1;
v_vec = v_vec_1;

end

function g_vec_1_hat = g_vec_Est(X,s_vec_fix,v_vec_fix,NetConf,PathConf)

GeneLen = size(X,2);
vec_1 = X'*s_vec_fix;
vec_2 = PathConf.Gene2Path*v_vec_fix;


L_Path = PathConf.lambda_C;

A_Mat_all = sum(s_vec_fix.^2)*speye(GeneLen) + (NetConf.lambda_L)*(NetConf.Lap_mat);
b_vec_all = vec_1 + L_Path * vec_2;
g_vec_1_hat = A_Mat_all\b_vec_all;

end

function v_vec_hat = v_vec_Est(g_vec_fix,PathConf)

v_vec_hat = (PathConf.lambda_C/PathConf.lambda_V)*(PathConf.Gene2Path')*g_vec_fix;
v_vec_hat = max(v_vec_hat,0);

end


function s1_hat = S_vec_binary(X,g_vec_fix,leastGroupNum)
a = sum(g_vec_fix.^2);
b_vec = X*g_vec_fix;

vec_x = b_vec/a;

th_candi = sort(vec_x,'descend');
g_vec_th = min(th_candi(leastGroupNum),0.5);

s1_hat = (vec_x >= g_vec_th);
end

