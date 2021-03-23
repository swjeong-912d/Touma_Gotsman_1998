function [nidx nn] = compute_nrings(v_list, fs, n, A, cum_on, add_v)
% [nidx nn] = compute_nrings(v_list, fs, n, A, cum_on)
% cum_on : true/false value.
if nargin < 4
    A = [];
    cum_on = true;
    add_v = false;
elseif nargin < 5
    cum_on = true;
    add_v = false;
elseif nargin < 6
    add_v = false;
end


if isempty(fs) && isempty(A)
    error('Invalid input');
elseif isempty(A)
    if size(fs,1) < size(fs,2)
        fs = fs';
    end
    f1 = fs(:);
    f2 = fs(:,[2 3 1]); f2 = f2(:);

    A = sparse(f1,f2,ones(size(f1)));
    A = (A+A'>0);
end

nidx = cell(length(v_list),1);
nn = zeros(length(v_list),1);
for nring = 1:length(v_list)
    iring = v_list(nring);
    iring_prev = iring;
    for i = 1:n
        [iring_next, ~] = find(A(:,iring));
        if ~cum_on && i == n
            iring_prev = iring;
        end
            iring = [iring; iring_next];
        
        if i > 1
            iring = unique(iring);
        end
    end
    iring(ismember(iring,iring_prev))=[];
    if add_v
        nidx{nring} = [v_list(nring);iring];
    else
        nidx{nring} = iring;
    end
    nn(nring) = length(iring);
end
if length(nidx) == 1
    nidx = nidx{1};
end