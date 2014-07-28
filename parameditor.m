function [out1,out2]=parameditor(simcount)
%Sam Roy 9/5/2012
%parameditor is used to edit parameters for a parameter sweep in CHILD
%the output is the matrix to be used by runchild2
%simcount=5;
cd('C:\child_n\ChildExercises\paramsweep01\');
fid=fopen('paramtemplate.txt','r');
paramtemplate=textscan(fid,'%s','delimiter', '\t','whitespace','');
load faultfile.mat;

j=1;
for i=1:2:length(paramtemplate{1,1})-1
    paramfile(j,1)=paramtemplate{1,1}(i);
    paramfile(j,2)=paramtemplate{1,1}(i+1);
    paramfile{j,2}=str2num(paramfile{j,2});
    j=j+1;
end

paramfile(:,3:simcount+1)=repmat(paramfile(:,2),1,simcount-1);

fprintf('Welcome to parameditor.\n This function allows you to assign parameters of any type to any number of CHILD simulations.\n');
fprintf('Simply add your desired parameter values for each simulation you wish to run,\n');
fprintf('then type "return". You can access parameter cell data by calling the variable paramfile{row,column}\n'); 
fprintf('faultfiles need to be stacked, one for each simulation.\n')
keyboard

out1=paramfile; out2=faultfile;
fprintf('You can now run the function "runchild3" to execute your simulations in CHILD.\n Try to contain your excitement.\n');