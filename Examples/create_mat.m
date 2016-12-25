function [Itrain, Jtrain,Valtrain,Itest,Jtest,Valtest,m,n] = create_mat(train,test)
% Given the tables: data (contains all the data. The first columns lists
% all the user IDs and the second column lists all the movie IDs),
% train (contains only the revealed items with the same structure),
% test (contains the unrevealed items which are used for
% testing), we create s sparse representation of the data. m,n will 
% be the dimension of the original data. Itrain, Jtrain, Valtrain 
% refer to the row, columns indexes and values of the revealed items 
% in the train set. Similaly, Iteest,Jtest, Valtest refer to the analog items
% in the test set.

% number of users
% m = max(data(:,1));
m1 = max(train(:,1));
m2 = max(test(:,1));
m = max(m1,m2);
% number of movies
% n = max(data(:,2));
n1 = max(train(:,2));
n2 = max(test(:,2));
n = max(n1,n2);


Itrain = train(:,1);
Jtrain = train(:,2);
Valtrain = train(:,3);
mu = mean(Valtrain);
muTrain = repmat(mu,size(Valtrain));
Valtrain = Valtrain - muTrain;

% Atrain = sparse(Itrain,Jtrain,Valtrain,m,n);
% 
% Rmean = mean(Atrain,2);
% Cmean = mean(Atrain,1);

% for i=1:(size(train,1))
%     Atrain(Itrain(i),Jtrain(i)) = Atrain(Itrain(i),Jtrain(i))-((Rmean(Itrain(i))+Cmean(Jtrain(i)))/2);
% end


% [Itrain,Jtrain,Valtrain] = find(Atrain);


Itest = test(:,1);
Jtest = test(:,2);
Valtest = test(:,3);
muTest = repmat(mu,size(Valtest));
Valtest = Valtest - muTest;



% Atest = sparse(Itest,Jtest,Valtest,m,n);

% 
% for i=1:(size(test,1))
%     Atest(Itest(i),Jtest(i)) = Atest(Itest(i),Jtest(i))-((Rmean(Itrain(i))+Cmean(Jtrain(i)))/2);
% end


% [Itest,Jtest,Valtest] = find(Atest);



end

