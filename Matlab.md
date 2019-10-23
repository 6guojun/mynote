# Basic

```matlab
help rem;
lookfor square;
exit;
clear all; close all; clc; 
format long;
format short;
%% create section code
read file;
     pwd;
     ls;
    cd /Users/Liu/Documents/01Code/04Matlab


```

# Read and Save 

```matlab
% Open data file 
fid = fopen('Data.csv')
% Read Data
out = textscan(fid, '%f %f %f','HeaderLines',1,'Delimiter',',')
% Read numeric value
a=csvread('out2.csv');disp(a);
% Read data with string
pp = readtable('pp.csv','ReadVariableNames',1,'ReadRowNames',true,'');
a = importdata('Book1.csv');disp(a.data);
data = dlmread('data.csv','\t',6,3);
dataTable = readtable('out2.csv', 'Format', '%{MM/dd/uu HH:mm:ss}D%f%f')
% Change variable name
dataTable.Properties.VariableNames = {'date', 'col1', 'col2'};
fclose(fid)
[date, col1, col2] = deal(out{:});

% save
writetable(p,'updated_test.csv','WriteVariableNames',0)
csvwrite('save.csv',data);
```

# Plot 

```matlab	
% Extract data from readData
xData = readData{1,1}(:,1);
yData = readData{1,1}(:,2)
% Extract columns
a = xlsread(filename,'','ensemble:geneid')

% Plot data
f1 = figure(1);
cla; hold on; grid on;grid minor'
plot(xData,ydata,'k-');title('...');xlable('...')'ylable('');
end

```



# Function

```matlab
function [outupt]=name(parameters)
end

% Operator : +, -, *, /, ^
% Logic operator : >,<,>=,<=,==,~=
% Boolean operation : &, |, ~

```

# Data frame

```matlab
patient.name = 'John Doe';
matrix
tuple
```

