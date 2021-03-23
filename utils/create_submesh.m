function [vs_sub fs_sub] = create_submesh(vs,fs,subnodes,tight)
if nargin < 4
    tight = false;
end
if size(vs,1) ~=3 && size(vs,2) == 3
    vs =vs';
end
if size(fs,1) ~=3 && size(fs,2) == 3
    fs =fs';
end
if islogical( subnodes(1))
  subnodes = find(subnodes);
end
vn = max(max(fs));
if ~isempty(vs)
    vs_sub = vs(:,subnodes);
else
    vs_sub = [];
end
% 쓰는거만으로 구성된 삼각형들을 가져온다. 이들을 정리한다.
if length(subnodes) > vn/2
    extnodes = 1:vn;
    extnodes(subnodes) = [];
    idx = zeros(1,size(fs,2));
    for i = 1:size(fs,1)
        idx = or(idx,ismember(fs(i,:),extnodes(:)));
    end
    idx = ~idx;
else
    idx = ones(1,size(fs,2));
    for i = 1:size(fs,1)
        idx = and(idx,ismember(fs(i,:),subnodes(:)));
    end
end

fs_sub = fs(:,idx);
if tight
    subnodes2 = unique(fs_sub);
    noisy_idx = ~ismember(subnodes,subnodes2);
    noisy_nodes = subnodes(noisy_idx);
    if ~isempty(vs)
        vs_sub = vs(:,subnodes2);
    end
    t = zeros(length(vs),1);
    t(subnodes2) = 1:length(subnodes2);
    fs_sub = t(fs_sub);
else
    t = zeros(vn,1);
    t(subnodes) = 1:length(subnodes);
    fs_sub = t(fs_sub);
end


    