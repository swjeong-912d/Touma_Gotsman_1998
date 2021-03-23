function [success, print_list, vs_r, fs_r] = connectivity_encoding(fs,vs,print_mode)
addpath(genpath('C:\My_Research\4D_Mesh_Compression\Time_varying_mesh_comp\source_codes\utils'))
if size(vs,1) ~= 3
    vs = vs';
end
if size(fs,1) ~= 3
    fs = fs';
end
 t_ind = 0;
fs2 = fill_hole(fs);
vn = max(max(fs));
vn2 = max(max(fs2));
assert(length(unique(fs2)) == vn2);
v_info = cell(vn2,3);
[~,valence_n_list] = compute_nrings(unique(fs2),fs2,1);
v_info(:,1) = num2cell(1:vn2);
v_info(:,2) = num2cell(valence_n_list);% active_valance
v_info(:,3) = num2cell(valence_n_list);
visited_f = zeros(length(fs2),1);
f_order = [];
encoded_v = zeros(vn2,1);
al_stack = {cell(1),0};
print_list = cell(length(fs2),1);
nadd = 1;
nsplit = 2;
nmerge = 3;
ndummy = 4;
% print_mode = true;
iter = 0;

fs3 = [];
vs2 = vs;
a = compute_nrings(unique(fs2),fs2,1);
for i = (vn+1):vn2
vs2(:,end+1) = mean(vs(:,a{i}),2);
end

[HE_list, sepIdx] = make_HE(fs2,true);
v_info(:,4) = num2cell(sepIdx(1:vn2));
v_order = [];
HE_list(:,end+1) = 0; % add encoded_edge info
while prod(visited_f) == 0
    t_ind = find(visited_f == 0,1);
    v1 = fs2(1,t_ind);
    v2 = fs2(2,t_ind);
    v3 = fs2(3,t_ind);
    fs3 = [fs3,[v1 v2 v3]'];
    encoded_v([v1,v2,v3]) = 1;
    visited_f(t_ind) = 1;
    f_order(end+1) = t_ind;
    v_order(end+1:end+3) = [v1 v2 v3];
    print_list{iter+1} = [nadd, v_info{v1,3}];
    print_encode(print_list{iter+1},print_mode);
    print_list{iter+2} = [nadd, v_info{v2,3}];
    print_encode(print_list{iter+2},print_mode);
    print_list{iter+3} = [nadd, v_info{v3,3}];
    print_encode(print_list{iter+3},print_mode);
    iter = iter+3;

    [Al, v_info, HE_list] = Add(v1,v2,v3,v_info,HE_list,sepIdx);
    
    al_stack = push(al_stack,Al);
    while ~empty(al_stack)
        Al = top(al_stack);
        al_stack = pop(al_stack);
        
        al = Al{1};
        al_end = Al{2};
        al_focus = Al{3};
        al_fsub = Al{4};
        v_info{Al{1}(al_end),2} = v_info{Al{1}(Al{2}),2} + Al{4};
        Al{4} = 0;
        while Al{2} > 0 % al_end > 0
            al = Al{1};
            al_end = Al{2};
            al_focus = Al{3};
            al_fsub = Al{4};
            e_start = v_info{al_focus,4};
            e = next_free_edge(HE_list,e_start,visited_f);
            if e ~= -1
                u = HE_list(e,2);
                if ~encoded_v(u)
                    if length(al) < al_end+1
                        al_reserve(al,30);
                    end
                    al(al_end+1) = u;
                    al_end = al_end+1;
                    if u <= vn
                        print_list{iter+1} = [nadd , v_info{u,3}];
                    else
                        print_list{iter+1} = [ndummy , v_info{u,3}];
                    end
                    iter = iter+1;
                    encoded_v(u) = 1;
                    v_order(end+1) = u;
                    assert(length(v_order) == nnz(encoded_v));
                    print_encode(print_list{iter},print_mode)
                    
                    e1 = HE_list(e,5);
                    HE_list(e1,6) = 1;
                    e2 = HE_list(e1,4);
                    HE_list(e2,6) = 1;
                    e3 = HE_list(e2,4);
                    HE_list(e3,6) = 1;
                    
                    
                    al_prev = HE_list(e2,2);
                    v_info = decrease_active_valance(v_info,al_prev,1);
                    v_info = decrease_active_valance(v_info,al_focus,1);
                    v_info = decrease_active_valance(v_info,u,1);
                    v_info = set_start_edge(v_info,al_focus,e);
                    assert(visited_f(HE_list(e1,3)) == 0);
                    f_order(end+1) = HE_list(e1,3);
                    visited_f(HE_list(e1,3)) = 1;
                    fs3 = [fs3,[al_focus, al_prev,u]'];
                    e_last = find_edge(HE_list,sepIdx,[al_focus, al(al_end)]);
                    assert(e_last ~= -1);
                    Al = {al,al_end,al_focus,al_fsub};
                elseif ismember(u,al(1:al_end)) % split
                    Al = {al,al_end,al_focus,al_fsub};
                    [Al, Al1, offset, inner_val, v_info] = split(Al,e,HE_list,v_info);
                    al_stack = push(al_stack,Al1);
                    print_list{iter+1} = [nsplit, offset, inner_val];
                    iter = iter+1;
                    print_encode(print_list{iter},print_mode)
%                     disp(['In encoder, before add triangle: ',num2str(length(f_order))]);
                    f_order(end+1) = HE_list(HE_list(e,5),3);
                    e_l = e;
                    e_r = HE_list(HE_list(e_l,5),4);
                     [visited_f,v_info,HE_list,fs3] = check_triangle(visited_f,v_info,HE_list,e_r,e_l,fs3);
                    
                else
                    [Al1, al_stack, k] = find_and_remove_al1(al_stack,u);
                    Al = {al,al_end,al_focus,al_fsub};
                    [Al, offset, v_info] = merge(Al,Al1,u, v_info);
                    print_list{iter+1} = [nmerge, k, offset];
                    iter = iter+1;
                    print_encode(print_list{iter},print_mode)
                    f_order(end+1) = HE_list(HE_list(e,5),3);
                    e_l = e;
                    e_r = HE_list(HE_list(e_l,5),4);
                    [visited_f,v_info,HE_list,fs3] = check_triangle(visited_f,v_info,HE_list,e_r,e_l,fs3);

                end
            end
            [Al] = remove_full_vertices(Al,v_info);
            al = Al{1};
            al_end = Al{2};
            al_focus = Al{3};
            al_fsub = Al{4};
            if al_end ~= 0 && is_full_vert(al_focus,v_info) % focus is full. 
                prev_focus = al_focus;
                al_focus = al(1);
                e_start = v_info{prev_focus,4};
                e_r = find_edge(HE_list,sepIdx,[prev_focus al(al_end)]);
                if e_r > 0
                    e_l = find_edge(HE_list,sepIdx,[prev_focus al_focus]);
                    if e_r ~= e_l
                        f_order(end+1) = HE_list(e_r,3);
                        [visited_f,v_info,HE_list,fs3] = check_triangle(visited_f,v_info,HE_list,e_r,e_l,fs3);
                    end
                end
                e_last = find_edge(HE_list,sepIdx,[al_focus, al(al_end)]);
                if e_last ~= -1
                    v_info = set_start_edge(v_info,al_focus,e_last);
                end
            end
            Al = {al,al_end,al_focus,al_fsub};
        end
    end
end
disp('Done!');
success = true;
print_list = print_list(1:iter);
inv_v_order(v_order) = 1:vn2;
fs_r = inv_v_order(fs3);
vs_r = vs2(:,v_order);
end
%%%
%%
function [visited_f,v_info,HE_list,fs3] = check_triangle(visited_f,v_info,HE_list,e_r,e_l,fs3)
    assert(visited_f(HE_list(e_r,3)) == 0);
    visited_f(HE_list(e_r,3)) = 1;

    HE_list(e_r,6) = 1;
    v_1 = HE_list(e_r,1);
    v_info = decrease_active_valance(v_info,v_1,1);

    e_n = HE_list(e_r,4);
    HE_list(e_n,6) = 1;
    v_2 = HE_list(e_r,2);
    v_info = decrease_active_valance(v_info,v_2,1);

    HE_list(HE_list(e_l,5),6) = 1;
    v_3 = HE_list(e_l,2);
    v_info = decrease_active_valance(v_info,v_3,1);
    fs3 = [fs3,[v_1 v_2 v_3]'];
end

%%%

%%
function [Al, v_info, HE_list] = Add(v1,v2,v3,v_info,HE_list,sepIdx)
    al = zeros(30,1);
    al(1) = v1;
    al(2) = v2;
    al(3) = v3;
    e1 = find_edge(HE_list,sepIdx,[v1,v2]);
    e2 = find_edge(HE_list,sepIdx,[v2,v3]);
    e3 = find_edge(HE_list,sepIdx,[v3,v1]);
    HE_list(e1,6) = 1;
    HE_list(e2,6) = 1;
    HE_list(e3,6) = 1;
    v_info = decrease_active_valance(v_info,v1,1);
    v_info = decrease_active_valance(v_info,v2,1);
    v_info = decrease_active_valance(v_info,v3,1);
    v_info = set_start_edge(v_info,v1,HE_list(e3,5));
    al_focus = v1;
    al_end = 3;
    al_fsub = 0;

    Al = {al,al_end,al_focus,al_fsub};
end
%%%

%%
function [Al1, al_stack, k] = find_and_remove_al1(al_stack,u)
top = al_stack{2};
for st_i = 1:top
    al1 = al_stack{1}{st_i}{1};
    if ~isempty(find(al1 == u,1))
        Al1 = al_stack{1}{st_i};
        al_stack{1}(st_i) = [];
        al_stack{2} = top-1;
        k = st_i;
        return;
    end
end
end
%%%

%%
function [Al, Al1, offset, inner_val, v_info] = split(Al, e, HE_list, v_info)
al = Al{1};
al_end = Al{2};
al_focus = Al{3};
al_fsub = Al{4};

%%%%%%%%%%% split start %%%%%%%%%%%
u = HE_list(e,2);
al1_focus = al_focus;
al_focus = u;
sep_ind = find(al == al_focus,1);
al1 = [al(1:sep_ind);zeros(30,1)];
al = al((sep_ind):end);
al1_end = sep_ind;
al_end = al_end-al1_end+1;
offset = al_end;

v1 = al1(1);
v2 = al1(al1_end);
inner_val = 1;
next_e = HE_list(e,4);
while HE_list(next_e,5) ~= e
    if HE_list(next_e,6) == 0 && HE_list(HE_list(next_e,5),6) == 1
        break;
    end
    inner_val = inner_val + 1;
    next_e = HE_list(HE_list(next_e,5),4);
end
al1_fsub = inner_val;
v_info = decrease_active_valance(v_info,u,al1_fsub);
v_info = set_start_edge(v_info,u,HE_list(e,5));
v_info = set_start_edge(v_info,al1_focus,e);
%%%%%%%%%%% split end %%%%%%%%%%%

Al = {al,al_end,al_focus,al_fsub};
Al1 = {al1,al1_end,al1_focus,al1_fsub};
end
%%%
%%
function e = find_edge(HE_list,sepIdx,e_info)
v1 = e_info(1);
v2 = e_info(2);
e = sepIdx(v1)-1+find(HE_list(sepIdx(v1):sepIdx(v1+1)-1,2) == v2);
if isempty(e)
    e = -1;
end
end
%%%


%%
function [Al, offset, v_info] = merge(Al,Al1,u,v_info)
al = Al{1};
al_end = Al{2};
al_focus = Al{3};
al_fsub = Al{4};
al1 = Al1{1};
al1_end = Al1{2};
al1_focus = Al1{3};
al1_fsub = Al1{4};

%%%%%%%%%%% merge start %%%%%%%%%%%

v_info{al1(al1_end),2} = v_info{al1(al1_end),2} + al1_fsub;
offset = find(al1 == u);
al = [al(1:al_end);al1(offset:al1_end);al1(1:offset);zeros(30,1)];
al_end = al_end+al1_end+1;
al_focus = al_focus;
%%%%%%%%%%% merge end %%%%%%%%%%%

Al = {al,al_end,al_focus,al_fsub};
end
%%%

%%
function [Al] = remove_full_vertices(Al,vert_info)
al = Al{1};
al_end = Al{2};
al_focus = Al{3};
al_fsub = Al{4};

%%%%%%%%%%% remove start %%%%%%%%%%%
    
fullvert_list = zeros(al_end,1);
if is_full_vert(al_focus,vert_info)
    fullvert_list(1)= 1;
else
    return;
end
for al_iter = 2:al_end
    al_v = al(al_iter);
    if vert_info{al_v,2} == 0
        fullvert_list(al_iter) = 1;
    else
        break;
    end
end
al(logical(fullvert_list)) = [];
al_end = al_end-nnz(fullvert_list);
%%%%%%%%%%% remove end %%%%%%%%%%%

Al = {al,al_end,al_focus,al_fsub};
end
%%%

%%
function answer = is_full_vert(v,vert_info)
answer = vert_info{v,2} < 2;
end
%%%
function v_info = decrease_active_valance(v_info,v,num)
v_info{v,2} = v_info{v,2}-num;
end
function v_info = set_start_edge(v_info,v,e)
v_info{v,4} = e;
end
%%
function e_next = next_free_edge(HE_list,e_start,visited_f)
e_next = e_start;
    if HE_list(e_next,6) == 0 && HE_list(HE_list(e_next,5),6) == 0
        return;
    end
    if visited_f(HE_list(e_next,3)) == 0
        e_next = HE_list(HE_list(HE_list(e_next,4),4),5);
        if HE_list(e_next,6) == 0 && HE_list(HE_list(e_next,5),6) == 0
            return;
        else
            e_next = -1;
            return;
            % cannot be further traveled by decoder but there's no free
            % edge.
        end
    else
        e_next = HE_list(HE_list(HE_list(e_next,4),4),5);
    end
while e_next ~= e_start 
    if HE_list(e_next,6) == 0 && HE_list(HE_list(e_next,5),6) == 0
        return;
    end
    if visited_f(HE_list(e_next,3)) == 0
        e_next = HE_list(HE_list(HE_list(e_next,4),4),5);
        if HE_list(e_next,6) == 0 && HE_list(HE_list(e_next,5),6) == 0
            return;
        else
            % cannot be further traveled by decoder but there's no free
            % edge.
            break;
        end
    end
    e_next = HE_list(HE_list(HE_list(e_next,4),4),5);
end
e_next = -1;
end
%%%

%%
function al = al_reserve(al,n)
al = [al;zeros(n,1)];
end
%%%

%%
function stack = push(stack,elem)
top = stack{2};
if length(stack{1}) < top+5
    stack{1}(end+1:end+5) = cell(1,1);
end
stack{1}(top+1) = {elem};
stack{2} = top+1;
end
%%%

%%
function stack = pop(stack)
top = stack{2};
stack{2} = top-1;
end
%%%

%%
function elem = top(stack)
top = stack{2};
elem = stack{1}{top};
end
%%%

%%
function answer = empty(stack)
answer = stack{2} == 0;
end
%%%

%%
function print_encode(print_list,print_mode)
word_list = {'add','split','merge','dummy'};
if print_mode
    fprintf('%s',word_list{print_list(1)});
    for num = print_list(2:end)
        fprintf(' %d',num);
    end
    fprintf('\n');
end
end
%%%