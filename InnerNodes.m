
function [ nodes ] = InnerNodes( sp,nodeIndex )

adjSp = sp(nodeIndex,:);
nodes = find(adjSp>0);

end

