function [X, y, headers] = readAndreasFile(fileName)
%reads andreas file from common folder and returns 
%input features output and the headers of the features and the output


fid = fopen(fileName);
line1 = fgets(fid);
headers = strsplit(line1, ',');
A = csvread(fileName, 1);
%get input features X and output y 
X = A(:,1:end-1);
y = A(:,end);


end