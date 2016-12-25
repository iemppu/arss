function [] =  write_c_data(u,trainName, testName)

% write_c_data(Itrain,Jtrain,Valtrain,Itest,Jtest,Valtest,trainName,testName)

% Divides the data into train and test
[utrain utest] = divide2trainTest(u);
% Create sparse matrix which represent the data
[Itrain, Jtrain,Valtrain,Itest,Jtest,Valtest,m,n] = create_mat(utrain,utest);


% This function writes the data into the text files trainName, testName row by row.
% Itrain - the row's indexes of the training data.
% Jtrain - the columns' indexes of the training data.
% Valtrain - the values of the training data.
% Itest - the row's indexes of the test data.
% Jtest - the columns' indexes of the test data.
% Valtest - the values of the test data.



m = max(max(Itrain),max(Itest));
n = max(max(Jtrain),max(Jtest));

fid = fopen(trainName,'wt');
fprintf(fid,'%d %d %d\n',m,n,length(Valtrain));
for i=1:length(Valtrain), 
    fprintf(fid,'%d %d %d\n',Itrain(i)-1,Jtrain(i)-1,Valtrain(i)); 
end
fclose(fid);

fid = fopen(testName,'wt');
fprintf(fid,'%d %d %d\n',m,n,length(Valtest));
for i=1:length(Valtest), 
    fprintf(fid,'%d %d %d\n',Itest(i)-1,Jtest(i)-1,Valtest(i)); 
end
fclose(fid);