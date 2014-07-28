function product=runchild5(basenm,paramfile,faultfile,sequence)

% SGR 9/2012
% Calls CHILD and runs a series of simulations (under 'basenm') with parameter inputs from
% the output from parameditor ('paramfile'). Will optionally run a model in
% sequence, with parameters from paramfile output.

ts=1e6;
% Indicate directories, startup files, run command
childexe = '\child_n\ChildExercises\paramsweep01\child';
basedir = '\child_n\ChildExercises\paramsweep01\';
rundir=mkdir(basedir, basenm);
rundir=([basedir basenm '\']);
mesh_assign='mesh_assign';
source='source';
copyfile([basedir mesh_assign '.in'],rundir);
copyfile([basedir source '.in'],rundir);
cd(rundir);
outfilename=basenm;
runcommand = [childexe ' ' rundir mesh_assign '.in'];

% Save parameter cell data as .mat and convert to matrix.
save('parameters','paramfile');
if exist('faultfile','var')==1
    save('faultfile','faultfile');
end
paramfile=cell2mat(paramfile(:,2:end));
total_run=length(paramfile(1,:));
% important for sequence option
oldts=0;
olduplift=0;
oldstep=0;

% Loop runs CHILD and produces results in separate files, one for each run
% total_run=5;
for i=1:total_run %For each column of input data
    prerunfile=([basedir 'sskp11\' sprintf('%.f',i) '\sskp11_' sprintf('%.f',i) ]);
    % read in file, save template
    % PREPARE PARAMETERS FOR SOURCE
    % layer construction and actual model run occur here.
    intemplate=creadinfile([rundir source '.in'],1);
    in=intemplate;
    in=cchangeparam(in,'RUNTIME',paramfile(1,i));
    in=cchangeparam(in,'OPINTRVL',paramfile(2,i));
    in=cchangeparam(in,'OPTREADINPUT',1);
    in=cchangeparam(in,'INPUTTIME',ts);
    in=cchangeparam(in,'INPUTDATAFILE',prerunfile);
    in=cchangeparam(in,'SEED',paramfile(3,i));
    in=cchangeparam(in,'OPTINITMESHDENS',paramfile(4,i));
    in=cchangeparam(in,'OUTFILENAME',source);
    in=cchangeparam(in,'X_GRID_SIZE',paramfile(5,i));
    in=cchangeparam(in,'Y_GRID_SIZE',paramfile(6,i));
    in=cchangeparam(in,'OPT_PT_PLACE',paramfile(7,i));
    in=cchangeparam(in,'GRID_SPACING',paramfile(8,i));
    in=cchangeparam(in,'NUM_PTS',paramfile(9,i));
    in=cchangeparam(in,'TYP_BOUND',paramfile(10,i));
    in=cchangeparam(in,'OUTLET_X_COORD',paramfile(11,i));
    in=cchangeparam(in,'OUTLET_Y_COORD',paramfile(12,i));
    in=cchangeparam(in,'MEAN_ELEV',paramfile(13,i));
    in=cchangeparam(in,'RAND_ELEV',paramfile(14,i));
    in=cchangeparam(in,'SLOPED_SURF',paramfile(15,i));
    in=cchangeparam(in,'UPPER_BOUND_Z',paramfile(16,i));
    in=cchangeparam(in,'OPTVAR',paramfile(17,i));
    in=cchangeparam(in,'ST_PMEAN',paramfile(18,i));
    in=cchangeparam(in,'ST_STDUR',paramfile(19,i));
    in=cchangeparam(in,'ST_ISTDUR',paramfile(20,i));
    in=cchangeparam(in,'OPTSINVARINFILT',paramfile(21,i));
    in=cchangeparam(in,'OPTMEANDER',paramfile(22,i));
    in=cchangeparam(in,'OPTDETACHLIM',paramfile(23,i));
    in=cchangeparam(in,'OPTREADLAYER',paramfile(24,i));
    in=cchangeparam(in,'OPTLAYEROUTPUT',paramfile(25,i));
    in=cchangeparam(in,'OPTINTERPLAYER',paramfile(26,i));
    in=cchangeparam(in,'FLOWGEN',paramfile(27,i));
    in=cchangeparam(in,'LAKEFILL',paramfile(28,i));
    in=cchangeparam(in,'TRANSMISSIVITY',paramfile(29,i));
    in=cchangeparam(in,'INFILTRATION',paramfile(30,i));
    in=cchangeparam(in,'OPTINLET',paramfile(31,i));
    in=cchangeparam(in,'OPTTSOUTPUT',paramfile(32,i));
    in=cchangeparam(in,'TSOPINTRVL',paramfile(33,i));
    in=cchangeparam(in,'OPTSTRATGRID',paramfile(34,i));
    in=cchangeparam(in,'DETACHMENT_LAW',paramfile(35,i));
    in=cchangeparam(in,'TRANSPORT_LAW',paramfile(36,i));
    in=cchangeparam(in,'KF',paramfile(37,i));
    in=cchangeparam(in,'MF',paramfile(38,i));
    in=cchangeparam(in,'NF',paramfile(39,i));
    in=cchangeparam(in,'PF',paramfile(40,i));
    in=cchangeparam(in,'KB',paramfile(41,i));
    in=cchangeparam(in,'KR',paramfile(42,i));
    in=cchangeparam(in,'KT',paramfile(43,i));
    in=cchangeparam(in,'MB',paramfile(44,i));
    in=cchangeparam(in,'NB',paramfile(45,i));
    in=cchangeparam(in,'PB',paramfile(46,i));
    in=cchangeparam(in,'TAUCB',paramfile(47,i));
    in=cchangeparam(in,'TAUCR',paramfile(48,i));
    in=cchangeparam(in,'KD',paramfile(49,i));
    in=cchangeparam(in,'DIFFUSIONTHRESHOLD',paramfile(50,i));
    in=cchangeparam(in,'OPTDIFFDEP',paramfile(51,i));
    in=cchangeparam(in,'OPT_NONLINEAR_DIFFUSION',paramfile(52,i));
    in=cchangeparam(in,'CRITICAL_SLOPE',paramfile(53,i));
    in=cchangeparam(in,'BEDROCKDEPTH',paramfile(54,i));
    in=cchangeparam(in,'REGINIT',paramfile(55,i));
    in=cchangeparam(in,'MAXREGDEPTH',paramfile(56,i));
    in=cchangeparam(in,'UPTYPE',paramfile(57,i));
    in=cchangeparam(in,'UPDUR',paramfile(58,i));
    in=cchangeparam(in,'UPRATE',paramfile(59,i));
    in=cchangeparam(in,'MINIMUM_UPRATE',paramfile(60,i));
    in=cchangeparam(in,'OPT_INCREASE_TO_FRONT',paramfile(61,i));
    in=cchangeparam(in,'NUMUPLIFTMAPS',paramfile(62,i));
    in=cchangeparam(in,'NUMGRNSIZE',paramfile(63,i)); %numg
    in=cchangeparam(in,'REGPROPORTION1',paramfile(64,i));
    in=cchangeparam(in,'BRPROPORTION1',paramfile(65,i));
    in=cchangeparam(in,'GRAINDIAM1',paramfile(66,i));
    in=cchangeparam(in,'REGPROPORTION2',paramfile(67,i));
    in=cchangeparam(in,'BRPROPORTION2',paramfile(68,i));
    in=cchangeparam(in,'GRAINDIAM2',paramfile(69,i));
    in=cchangeparam(in,'BETA',paramfile(70,i));
    in=cchangeparam(in,'HIDINGEXP',paramfile(71,i));
    in=cchangeparam(in,'CHAN_GEOM_MODEL',paramfile(72,i));
    in=cchangeparam(in,'HYDR_WID_COEFF',paramfile(73,i));
    in=cchangeparam(in,'HYDR_WID_EXP_DS',paramfile(74,i));
    in=cchangeparam(in,'HYDR_WID_EXP_STN',paramfile(75,i));
    in=cchangeparam(in,'HYDR_SLOPE_EXP',paramfile(76,i));
    in=cchangeparam(in,'HYDR_DEP_COEFF_DS',paramfile(77,i));
    in=cchangeparam(in,'HYDR_DEP_EXP_DS',paramfile(78,i));
    in=cchangeparam(in,'HYDR_DEP_EXP_STN',paramfile(79,i));
    in=cchangeparam(in,'HYDR_ROUGH_COEFF_DS',paramfile(80,i));
    in=cchangeparam(in,'HYDR_ROUGH_EXP_DS',paramfile(81,i));
    in=cchangeparam(in,'HYDR_ROUGH_COEFF_STN',paramfile(82,i));
    in=cchangeparam(in,'HYDR_ROUGH_EXP_STN',paramfile(83,i));
    in=cchangeparam(in,'BANK_ROUGH_COEFF',paramfile(84,i));
    in=cchangeparam(in,'BANKFULLEVENT',paramfile(85,i));
    in=cchangeparam(in,'OPTFLOODPLAIN',paramfile(86,i));
    in=cchangeparam(in,'OPTLOESSDEP',paramfile(87,i));
    in=cchangeparam(in,'OPTEXPOSURETIME',paramfile(88,i));
    in=cchangeparam(in,'OPTVEG',paramfile(89,i));
    in=cchangeparam(in,'OPTKINWAVE',paramfile(90,i));
    in=cchangeparam(in,'OPTMESHADAPTZ',paramfile(91,i));
    in=cchangeparam(in,'MESHADAPT_MAXNODEFLUX',paramfile(92,i));
    in=cchangeparam(in,'OPTMESHADAPTAREA',paramfile(93,i));
    in=cchangeparam(in,'OPTFOLDDENS',paramfile(94,i));
    
    fid=cwriteinfile([rundir source '.in'],in); % write out the input file for use
    
    
    cd(rundir);
    % if models are sequential and initial values are taken from previous run, any steps after 1 if nargin=4 go here
    % How many faults in the fault file??
    [s,w]=system(['cp ' prerunfile '.lay10 ' prerunfile '.lay']);
    outfilename=sprintf('%s_%.f',outfilename,i); outfoldername=sprintf('%.f\',i);
    in=cchangeparam(in,'OUTFILENAME',outfilename);
    fid=cwriteinfile([rundir source '.in'],in);% outfile name changes with each iteration
    
    outdir=mkdir(rundir, outfoldername); outdir=([rundir outfoldername]);% make folder to hold CHILD outputs from iteration i
    
    runcommand = [childexe ' ' rundir source '.in'];
    if strcmp(mysys,'unix')
        [s,w]=unix(runcommand);
    elseif strcmp(mysys,'dos')
        [s,w]=dos(runcommand,'-echo'); % run the model
    else
        error('You need to set mysys to either unix or dos, sucka');
    end
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


