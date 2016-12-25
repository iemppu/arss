function [utrain utest] = divide2trainTest(u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

m = max(u(:,1));


u = u(randperm(size(u,1)),:);
itr = 1;
its = 1;

for i=1:m
    tr = ceil(size(find(u(:,1)==i),1)/2);
    ts = floor(size(find(u(:,1)==i),1)/2);
    if (tr~=0)
        Itrain(itr:(itr+tr-1)) = find(u(:,1)==i,tr,'first');
        Itest(its:(its+ts-1)) = find(u(:,1)==i,ts,'last');
        itr = itr + tr;
        its = its + ts;
        
    end
end


utrain = u(Itrain,:);
utest = u(Itest,:);


end

