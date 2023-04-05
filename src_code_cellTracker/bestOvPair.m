function flag = bestOvPair(movieInfo, id1, id2)
% find if id1 and id2 are the nereast neighbors to each other in the
% corresponding frames

flag = false;

f1 = movieInfo.frames(id1);
f2 = movieInfo.frames(id2);
if f1 == f2
    return;
end
if f1>f2
    id = id1;
    id1 = id2;
    id2 = id;
    
    f = f1;
    f1 = f2;
    f2 = f;
end

% [bestNei, ~] = bestOvNei(id1, movieInfo, f2);
% 
% if id2 ~= bestNei
%     return;
% end
% 
% [bestNei, ~] = bestOvNei(id2, movieInfo, f1, true);
% 
% if id1 ~= bestNei
%     return;
% end

% Because of non-rigid registration, overlap measure is not good enough.
% Here we still use CDist, and if child is the nearest neighbor of the
% parent, it's good enough to be a bestNei.
candi_nei = movieInfo.frames(movieInfo.nei{id1})==f2;
[~, min_nei] = min(movieInfo.CDist{id1}(candi_nei));
bestNei = movieInfo.nei{id1}(min_nei);

if id2 ~= bestNei
    return;
end

flag = true;


