function Ind_vec = GetInd_in_Dict(List_a,dict_b)
% List_a = {'b','c','d','b','a','e'};
% dict_b = {'c','b','e'};


if size(List_a,2) ~= 1
    List_a = List_a';
end
if size(dict_b,2) ~= 1
    dict_b = dict_b';
end
Len_Dict = size(dict_b,1);

[~,ind_remain] = setdiff(List_a,dict_b);

mark_remain = 0;
if ~isempty(ind_remain)
    mark_remain = 1;
    if isnumeric(dict_b)
        new_elem = 0;
    else
        new_elem = {'&*D&*SEF#'};
    end
    dict_b = [dict_b; new_elem];
    
    while ~isempty(ind_remain)
        List_a(ind_remain) = new_elem;
        [~,ind_remain] = setdiff(List_a,dict_b);
    end
end

[Sort_a,ind_a] = sort(List_a);
[~,ind_a_inv] = sort(ind_a);
[Sort_b,ind_b] = sort(dict_b);
[~,~,ib] = intersect(Sort_a,Sort_b);
dict_num = ind_b(ib);

[~,~,ind_123] = unique(Sort_a);
Vec_sort = dict_num(ind_123);
Ind_vec = Vec_sort(ind_a_inv);

if mark_remain == 1
    Ind_vec(Ind_vec == (Len_Dict+1)) = 0;
end