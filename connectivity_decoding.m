%%

function [vs,fs] = connectivity_decoding(print_list, print_mode,vs2,fs_e)
addpath(genpath('C:\My_Research\4D_Mesh_Compression\Time_varying_mesh_comp\source_codes\utils'))
subnode_list = [];
al_stack = {cell(1),0};
fs2 = [];
iter = 0;
HE_list = [];
sepIdx = [1];
v_info = {};
al_stack = {cell(1),0};
vn = 0;
nadd = 1;
nsplit = 2;
nmerge = 3;
ndummy = 4;
if nargin < 2
    print_mode = false;
end
dummy_list = [];
iter = 0;
while iter < length(print_list)
    v1 = vn+1;
    v2 = vn+2;
    v3 = vn+3;
    vn = vn+3;
    print_encode(print_list{iter+1},print_mode);
    print_encode(print_list{iter+2},print_mode);
    print_encode(print_list{iter+3},print_mode);
    val1 = print_list{iter+1}(2);
    val2 = print_list{iter+2}(2);
    val3 = print_list{iter+3}(2);
    iter = iter+3;
    [fs2, Al, v_info, HE_list, sepIdx] = Addt(fs2,v1,v2,v3,val1,val2,val3,v_info,HE_list,sepIdx);
    
    al_stack = push(al_stack,Al);
    while ~empty(al_stack)
        Al = top(al_stack);
        al_stack = pop(al_stack);
        v_info{Al{1}(Al{2}),2} = v_info{Al{1}(Al{2}),2} + Al{4};
        Al{4} = 0;
        while Al{2} > 0 % al_end > 0
            al = Al{1};
            al_end = Al{2};
            al_focus = Al{3};
            al_fsub = Al{4};
            e_start = v_info{al_focus,4};
            [e,HE_list] = next_free_edge(HE_list,e_start,sepIdx,v_info);
            if e ~= -1
                iter = iter+1;
                print_encode(print_list{iter},print_mode);
                switch print_list{iter}(1)
                    case {nadd, ndummy}
                        if length(al) < al_end+1
                            al = reserve(al,30);
                        end
                        val = print_list{iter}(2);
                        u = vn+1;
                        if print_list{iter}(1) == ndummy
                            dummy_list(end+1) = u;
                        end
                        vn = vn+1;
                        [fs2, Al, v_info, HE_list, sepIdx] = Adde(fs2,Al,u,val,e,e_start,v_info,HE_list,sepIdx);
                        v_info = set_start_edge(v_info,al_focus,e);
                    case nsplit
                        %                         disp(['In decoder, before add triangle: ',num2str(length(fs2))]);
                        offset = print_list{iter}(2);
                        fsub = print_list{iter}(3);
                        [Al, Al1, v_info] = split(Al,offset,fsub,v_info);
                        al_stack = push(al_stack,Al1);
                        HE_list(e,2) = Al{1}(1);
                        ep = find_edge(HE_list,sepIdx,HE_list(e,[2 1]));
                        if ep == -1
                            [ep,HE_list] = add_edge(HE_list,sepIdx,HE_list(e,[2 1]));
                        end
                        HE_list([ep,e],5) = [e,ep];
                        e_r = e_start;
                        e_l = e;
                        [fs2,v_info,HE_list] = check_triangle(fs2,e_r,e_l,v_info,HE_list,sepIdx);
                        e_last = find_edge(HE_list,sepIdx,[Al{3}, Al{1}(Al{2})]);
                        if e_last ~= -1
                            v_info = set_start_edge(v_info,Al{3},e_last);
                        end
                        v_info = set_start_edge(v_info,HE_list(e,1),e);
                    case nmerge
                        k = print_list{iter}(2);
                        offset = print_list{iter}(3);
                        Al1 = al_stack{1}{k};
                        al_stack{1}(k) = [];
                        al_stack{2} = al_stack{2}-1;
                        Al = {al,al_end,al_focus,al_fsub};
                        [Al, v_info] = merge(Al,Al1,offset, v_info);
                        
                        
                        HE_list(e,2) = Al{1}(Al{2});
                        ep = find_edge(HE_list,sepIdx,HE_list(e,[2 1]));
                        if ep == -1
                            [ep,HE_list] = add_edge(HE_list,sepIdx,HE_list(e,[2 1]));
                        end
                        HE_list([ep,e],5) = [e,ep];
                        v_info = set_start_edge(v_info,HE_list(e,1),e);

                        e_r = e_start;
                        e_l = e;
                        [fs2,v_info,HE_list] = check_triangle(fs2,e_r,e_l,v_info,HE_list,sepIdx);
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
                        [fs2,v_info,HE_list] = check_triangle(fs2,e_r,e_l,v_info,HE_list,sepIdx);
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
subnode_list = 1:vn;
subnode_list(dummy_list) = [];
[vs,fs] = create_submesh(vs2,fs2,subnode_list);
disp('Done!');

end
%%%

%%
function [fs2,v_info,HE_list] = check_triangle(fs2,e_r,e_l,v_info,HE_list,sepIdx)
v_f = HE_list(e_r,1);
v_1 = HE_list(e_r,2);
v_2 = HE_list(e_l,2);

v_info{v_f,2} = v_info{v_f,2} - 1;
v_info{v_1,2} = v_info{v_1,2} - 1;
v_info{v_2,2} = v_info{v_2,2} - 1;

ef1 = e_r;
ef2 = e_l;
e2f = HE_list(ef2,5);
e12 = find_edge(HE_list,sepIdx,[v_1 v_2]);
e21 = find_edge(HE_list,sepIdx,[v_2 v_1]);
if e12 == -1
    [e12,HE_list] = add_edge(HE_list,sepIdx,[v_1 v_2]);
end
if e21 == -1
    [e21,HE_list] = add_edge(HE_list,sepIdx,[v_2 v_1]);
end

fs2 = [fs2,[v_f; v_1; v_2]];
fn = length(fs2);
HE_list([ef1,e12,e2f],3) = fn;
HE_list([ef1,e12,e2f],4) = [e12,e2f,ef1];
HE_list([e12,e21],5) = [e21,e12];
HE_list([ef1,e12,e2f],6) = 1;
end
%%%

%%
function [fs2,Al,v_info,HE_list,sepIdx] = Adde(fs2,Al,u,val,e,e_start,v_info,HE_list,sepIdx)
Al{1}(Al{2}+1) = u;
Al{2} = Al{2}+1;


HE_list(e,2) = u;
sepIdx(end+1) = sepIdx(end)+val;
if size(HE_list,1) < max(sepIdx)
    HE_list = [HE_list;zeros(val*10000,6)];
end
vn = length(sepIdx)-1;
if size(v_info,1) < vn
    v_info = [v_info;cell(10000,4)];
end
v_info(vn,:) = {u,val,val,sepIdx(u)};

v_f = Al{3};
v_u = u;
v_p = HE_list(e_start,2);
v_info = decrease_active_valance(v_info,v_f,1);
v_info = decrease_active_valance(v_info,v_u,1);
v_info = decrease_active_valance(v_info,v_p,1);
efu = e;
efp = e_start;
euf = find_edge(HE_list,sepIdx,[v_u v_f]);
[euf,HE_list] = add_edge(HE_list,sepIdx,[v_u v_f]);
[eup,HE_list] = add_edge(HE_list,sepIdx,[v_u v_p]);
[epu,HE_list] = add_edge(HE_list,sepIdx,[v_p v_u]);

fs2 = [fs2,[v_f;v_p;v_u]];
fn = length(fs2);
HE_list([efp,epu,euf],3) = fn;
HE_list([efp,epu,euf],4) = [epu,euf,efp];
HE_list([epu,eup, euf,efu],5) = [eup,epu, efu,euf];
HE_list([efp,epu,euf],6) = 1;
end

%%
function [fs2, Al, v_info, HE_list, sepIdx] = Addt(fs2,v1,v2,v3,val1,val2,val3,v_info,HE_list,sepIdx)
vn = length(sepIdx)-1;
fs2 = [fs2,[v1;v2;v3]];
fn = length(fs2);


if size(v_info,1) < vn
    v_info = [v_info;cell(10000,4)];
end
sepIdx(end+1) = sepIdx(end)+val1;
sepIdx(end+1) = sepIdx(end)+val2;
sepIdx(end+1) = sepIdx(end)+val3;
if size(HE_list,1) < max(sepIdx)
    HE_list = [HE_list;zeros(val1*10000,6)];
end
al = zeros(30,1);
al(1:3) = [v1,v2,v3];
[e1,HE_list] = add_edge(HE_list,sepIdx,[v1,v2]);
[e2,HE_list] = add_edge(HE_list,sepIdx,[v2,v3]);
[e3,HE_list] = add_edge(HE_list,sepIdx,[v3,v1]);
[ep1,HE_list] = add_edge(HE_list,sepIdx,[v2,v1]);
[ep2,HE_list] = add_edge(HE_list,sepIdx,[v3,v2]);
[ep3,HE_list] = add_edge(HE_list,sepIdx,[v1,v3]);

HE_list([e1,e2,e3],4) = [e2,e3,e1];
HE_list([e1,ep1, e2,ep2, e3,ep3],5) = [ep1,e1, ep2,e2, ep3,e3];
HE_list([e1 e2 e3],3) = fn;
HE_list([e1,e2,e3],6) = 1;

v_info(vn+1,:) = {v1,val1-1,val1,sepIdx(v1)};
v_info(vn+2,:) = {v2,val2-1,val2,sepIdx(v2)};
v_info(vn+3,:) = {v3,val3-1,val3,sepIdx(v3)};
al_focus = v1;
al_end = 3;
al_fsub = 0;

v_info = set_start_edge(v_info,al_focus,ep3);

Al = {al,al_end,al_focus,al_fsub};
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
function e = empty_edge(HE_list,sepIdx,v)
e = sepIdx(v)-1+find(HE_list(sepIdx(v):sepIdx(v+1)-1,1) == 0,1);
if isempty(e)
    e = -1;
end
end
%%%
%%
function [e,HE_list] = add_edge(HE_list,sepIdx,e_info)
v1 = e_info(1);
v2 = e_info(2);
e = empty_edge(HE_list,sepIdx,v1);
if e ~= -1
    HE_list(e,[1,2]) = [v1,v2];
end
end
%%%



%%
function [Al, Al1, v_info] = split(Al,offset,fsub,v_info)
al = Al{1};
al_end = Al{2};
al_focus = Al{3};
al_fsub = Al{4};

%%%%%%%%%%% split start %%%%%%%%%%%
sep_ind = al_end-offset+1;
al1 = [al(1:sep_ind);zeros(30,1)];
al = al((sep_ind):end);
al1_end = sep_ind;

al_end = al_end-al1_end+1;
al_focus = al(1);
al1_focus = al1(1);
v_info{al_focus,2} = v_info{al_focus,2} - fsub;
al1_fsub = fsub;

%%%%%%%%%%% split end %%%%%%%%%%%

Al = {al,al_end,al_focus,al_fsub};
Al1 = {al1,al1_end,al1_focus,al1_fsub};
end
%%%

%%
function [Al,v_info] = merge(Al,Al1,offset,v_info)
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
al = [al(1:al_end);al1(offset:al1_end);al1(1:offset);zeros(30,1)];
al_end = al_end+al1_end+1;
al_focus = al_focus;
%%%%%%%%%%% merge end %%%%%%%%%%%

Al = {al, al_end, al_focus,al_fsub};
end
%%%

%%
function [Al] = remove_full_vertices(Al,v_info)
al = Al{1};
al_end = Al{2};
al_focus = Al{3};
al_fsub = Al{4};

%%%%%%%%%%% remove start %%%%%%%%%%%

fullvert_list = zeros(al_end,1);
if is_full_vert(al_focus,v_info)
    fullvert_list(1)= 1;
end
for al_iter = 2:al_end
    al_v = al(al_iter);
    if v_info{al_v,2} == 0
        fullvert_list(al_iter) = 1;
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
function [e_next, HE_list] = next_free_edge(HE_list,e_start,sepIdx,v_info)
v = HE_list(e_start,1);
e_next = e_start;
if HE_list(e_next,6) == 0 && HE_list(HE_list(e_next,5),6) == 0
    return;
end
if HE_list(e_next,4) > 0
    e_next = HE_list(HE_list(HE_list(e_next,4),4),5);
elseif ~is_full_vert(v,v_info)
    [e_next,HE_list] = add_edge(HE_list,sepIdx,[v 0]);
    return;
end
while e_next ~= e_start
    if HE_list(e_next,6) == 0 && HE_list(HE_list(e_next,5),6) == 0
        return;
    end
    if HE_list(e_next,4) > 0
        e_next = HE_list(HE_list(HE_list(e_next,4),4),5);
    elseif ~is_full_vert(v,v_info)
        [e_next,HE_list] = add_edge(HE_list,sepIdx,[v 0]);
        return;
    elseif is_full_vert(v,v_info)
        break;
    end
end
e_next = -1;
end
%%%
%%
function al = reserve(al,n)
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