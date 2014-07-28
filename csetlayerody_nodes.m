function [ld,nl]=csetlayerody_nodes( fname, ts, node_erody, numg, format_version )
% CSETLAYERODY: Sets the erodibility value in CHILD layer data. Use this
% when you want to restart a CHILD run using spatially variable
% erodibility. The function reads layer data from a CHILD run, sets the
% erodibility coefficient at each node according to the corresponding value
% in "node_erody", which must be a vector of the same length as the number
% of interior nodes, and in the same order. All layers at each node are
% assigned the same erodibility value.
%
% Inputs: fname - base name of file from which to read data
%         ts - time slice to read
%         node_erody - erodibility coefficient for each node, in order
%                      (a vector)
%         numg (optional) - number of grain-size classes (default=1)
%         format_version (optional) - version of layer format to use (0 or
%                                       1; defaults to 0)
%
% Returns: modified layer data ld, number of layers at each node nl
%
% GT, Nov 2010

% Check numg, default to 1 if not specified
if nargin<9
    numg=1;
end
if nargin<10
    format_version = 0;
end

%added for debugging 7/3/12
numg=1;
% Read in initial data (assumes the name of the file is <fname>.lay<ts-1>).
% The layer file will normally have been generated from a previous CHILD
% run.
%[ld,nl,time]=creadlayers(fname,ts,numg,format_version);
%removed last argument format_version for debugging 7/2/12
[ld,nl,time]=creadlayers(fname,ts,numg);

% Write the results to a file
cwritelayers(ld, nl, [fname '_mod'], ts, node_erody, 10000, numg, time);
