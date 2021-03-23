function fs2 = fill_hole(fs)
if size(fs,1) ~= 3
    fs = fs';
end
number_of_original_vertices = max(max(fs));
number_of_holefilled_vertices = number_of_original_vertices;
[HE_list sepIdx] = make_HE(fs,true);
boundary_edge_index = find(HE_list(:,5) == 0);
number_of_boundary_edges = length(boundary_edge_index);
vertices_in_boundary = unique(HE_list(boundary_edge_index,1:2));
boundary_HE_sublist = HE_list(boundary_edge_index,:);

% boundary edges separation
sepbIdx = zeros(length(vertices_in_boundary)+1,1);
sepbIdx(1) = 1;
vi = 2;
for i = 2:number_of_boundary_edges
    if boundary_HE_sublist(i-1,1) ~= boundary_HE_sublist(i,1)
        v = boundary_HE_sublist(i,1);
        sepbIdx(vi) = i;
        vi = vi+1;
    end
end
sepbIdx(vi) = number_of_boundary_edges+1;

original_v_index_to_boundary_v_index(vertices_in_boundary) = 1:length(vertices_in_boundary);



hole_filled = zeros(number_of_boundary_edges,1);
fs_add = [];
while sum(hole_filled) < number_of_boundary_edges
    boundary_edge_with_hole = find(hole_filled==0,1);
    boundary_edges_with_hole = [];
    edges_with_unfilled_hole = [boundary_edge_with_hole];
    number_of_edges_with_hole = 1;
    while number_of_edges_with_hole > 0
        boundary_edges_list = edges_with_unfilled_hole;
        edges_with_unfilled_hole = []; % reset
        for b_i = 1:number_of_edges_with_hole
            boundary_e_ind = boundary_edges_list(b_i);
            boundary_v_ind = original_v_index_to_boundary_v_index(boundary_HE_sublist(boundary_e_ind,1:2));
            edges_with_unfilled_hole = [edges_with_unfilled_hole,...
                            unique([sepbIdx(boundary_v_ind(1)):(sepbIdx(boundary_v_ind(1)+1)-1),sepbIdx(boundary_v_ind(2)):(sepbIdx(boundary_v_ind(2)+1)-1)])];
        end
        edges_with_unfilled_hole = unique(edges_with_unfilled_hole);
        edges_with_unfilled_hole = edges_with_unfilled_hole(~ismember(edges_with_unfilled_hole,boundary_edges_with_hole));
        number_of_edges_with_hole = length(edges_with_unfilled_hole);
        boundary_edges_with_hole = [boundary_edges_with_hole,edges_with_unfilled_hole];
    end
    fs_add = [fs_add,[boundary_HE_sublist(boundary_edges_with_hole,2:-1:1),repmat(number_of_holefilled_vertices+1,length(boundary_edges_with_hole),1)]'];
    hole_filled(boundary_edges_with_hole) = 1;
    number_of_holefilled_vertices = number_of_holefilled_vertices+1;
end
fs2 = [fs,fs_add];