function [OPTIONS_STRUCT,OPTIONS_CELL,OPTIONS_NAMES]=kinect_get_defaults(varargin)
%
%
%
%

[pathname,~,~]=fileparts(mfilename('fullpath'));
def_options=fullfile(pathname,'kinect_analysis_defaults.config');

OPTIONS_STRUCT=kinect_read_options(def_options,varargin);
OPTIONS_CELL=kinect_map_parameters(OPTIONS_STRUCT);
OPTIONS_NAMES=fieldnames(OPTIONS_STRUCT);
