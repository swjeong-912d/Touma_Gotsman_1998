function [HE, sepIdx, HE2, sepIdx2] = make_HE(fs, perform_pairing)
% [HE, sepIdx, HE2, sepIdx2] = make_HE(fs, perform_pairing)
%
% Inputs
% fs: face
% perform_paring: If it is false, we do not perform pairing operation
% 
% Outputs
%
% HE: Half edges sorted in first column. 
%     structure is [v1 v2 face_idx next_idx pair_idx]
%     if perform_paring is false, pair_idx will not exist
%
% sepIdx: seperation index. 
%         sepIdx(i):sepIdx(i+1)-1 means HE with v1 == i
% 
% HE2: Half edges sorted in second column. 
%      structure is [v1 v2 face_idx next_idx pair_idx]
%      if perform_paring is false, pair_idx will not exist
%
% sepIdx2: seperation index. 
%         sepIdx2(i):sepIdx2(i+1)-1 means HE2 with v2 == i

% HE: v1 v2 face next pair
if size(fs,1) ~= 3
    fs = fs';
end
if nargin < 2
    perform_pairing = true;
end
vn = max(max(fs));
fn = size(fs,2);
% Edges = zeros(fn*6,2);
HE = zeros(fn*3,4);
for j = 1 : 3
    a = fs(j,:);
    b = fs(mod(j,3)+1,:);
    HE(fn*(j-1)+(1:fn),1) = a;
    HE(fn*(j-1)+(1:fn),2) = b;
    HE(fn*(j-1)+(1:fn),3) = 1:fn;
    HE(fn*(j-1)+(1:fn),4) = fn*mod(j,3)+(1:fn);
end
[HE, idx1] = sortrows(HE,1);
idx1r(idx1) = 1:length(idx1);
% pointer readjustment
HE(:,4) = idx1r(HE(:,4));
sepIdx = ones(vn+1,1);
cnt = 1;
for i = 2:length(HE)
    if HE(i,1) ~= HE(i-1,1)
        cnt = cnt+1;
        sepIdx(cnt) = i;
    end
end
sepIdx(cnt+1) = length(HE)+1;

if perform_pairing
    for i = 1:length(HE)
        v1 = HE(i,1);
        v2 = HE(i,2);
        pair = find(HE(sepIdx(v2):sepIdx(v2+1)-1,2) == v1);
        if ~isempty(pair)
            HE(i,5) =  pair+sepIdx(v2)-1;
        end
    end
end

if nargout > 2
    [HE2, idx2] = sortrows(HE,2);
    idx2rev(idx2) = 1:length(idx2);
    idx2rev(end+1) = 0;
    t = 4;
    if perform_pairing
        t = [t,5];
    end
        % pair가 없으면 오류를 방지하기 위해 dummy data와 링크해서 0을 만든다.
        tmp = HE2(:,t);
        tmp(tmp == 0) = length(idx2rev);
        HE2(:,t) = idx2rev(tmp);
    sepIdx2 = ones(max(max(fs))+1,1);
    cnt = 1;
    for i = 2:length(HE2)
        if HE2(i-1,2) ~= HE2(i,2)
            cnt = cnt+1;
            sepIdx2(cnt) = i;
        end
    end
    sepIdx2(cnt+1) = length(HE2)+1;
end
end