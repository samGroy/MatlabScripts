function product=runchild4(basenm,inputfile,sequence)

% SGR 4/2014
% SIMPLE VERSION FOR UNDERGRADS. only unfixed parameters are runtime, optintervl,x and y length, grid spacing, precip rate, uplift rate, and erodibility. 

% Indicate directories, startup files, run command
childexe = '\child_n\ChildExercises\paramsweep01\child';
basedir = '\child_n\ChildExercises\paramsweep01\';
rundir=mkdir(basedir, basenm);
rundir=([basedir basenm '\']);
mesh_assign='mesh_assign';
source='source_runchild4';
copyfile([basedir mesh_assign '.in'],rundir);
copyfile([basedir source '.in'],rundir);
cd(rundir);
outfilename=basenm;
runcommand = [childexe ' ' rundir mesh_assign '.in'];

% Save parameter cell data as .mat and convert to matrix.
save('parameters','inputfile');
if exist('faultfile','var')==1
    save('faultfile','faultfile');
end
inputfile=cell2mat(inputfile(:,2:end));
total_run=length(inputfile(1,:));
% important for sequence option
oldts=0;
olduplift=0;
oldstep=0;


% Loop runs CHILD and produces results in separate files, one for each run
% total_run=5;
for i=1:total_run %For each column of input data
    % PREPARE PARAMETERS FOR SOURCE
    % layer construction and actual model run occur here.
    intemplate=creadinfile([rundir source '.in'],1);
    in=intemplate;
    in=cchangeparam(in,'RUNTIME',inputfile(1,i));
    in=cchangeparam(in,'OPINTRVL',inputfile(2,i));
    in=cchangeparam(in,'OPTREADINPUT',0);
    in=cchangeparam(in,'GRID_SPACING',inputfile(4,i));
    in=cchangeparam(in,'X_GRID_SIZE',inputfile(3,i));
    in=cchangeparam(in,'Y_GRID_SIZE',inputfile(3,i));
    in=cchangeparam(in,'ST_PMEAN',inputfile(5,i));
    in=cchangeparam(in,'UPRATE',inputfile(6,i));
    in=cchangeparam(in,'KB',inputfile(7,i));
    in=cchangeparam(in,'KD',inputfile(8,i));
    in=cchangeparam(in,'OPTREADINPUT',10);
    in=cchangeparam(in,'ST_STDUR',inputfile(2,i)/50);
    
    if nargin==3 && i>=2 % if sequence==1, build from previous run.
        ts=inputfile(1,i-1)+oldts;
        in=cchangeparam(in,'OPTREADINPUT',1);
        in=cchangeparam(in,'INPUTDATAFILE',oldoutfilename);
        in=cchangeparam(in,'INPUTTIME',ts);
        oldts=ts;
    end
    
    outfilename=sprintf('%s_%.f',outfilename,i); outfoldername=sprintf('%.f\',i);
    in=cchangeparam(in,'OUTFILENAME',outfilename);

    fid=cwriteinfile([rundir source '.in'],in); % write out the input file for use
    
    cd(rundir);
    
    outdir=mkdir(rundir, outfoldername); outdir=([rundir outfoldername]);% make folder to hold CHILD outputs from iteration i
    
    runcommand = [childexe ' ' rundir source '.in']; 
    [s,w]=system(runcommand,'-echo'); % run the model
    if s>0
        fprintf('Error message from CHILD:\n\n%s',w);
        error('Bailing out of runflacchild.m');
    end
    
    if exist('sequence','var')==0 % organize files by model run number
        movefile([rundir outfilename '.*'],outdir);
    end
    oldoutfilename=outfilename;
    outfilename=basenm;
    
end


