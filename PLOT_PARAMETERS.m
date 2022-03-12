%==========================================================================
%
% PLOT_PARAMETERS  Initializes plot parameters.
%
%   pp = PLOT_PARAMETERS
%
% Author: Tamas Kis
% Last Update: 2022-03-02
%
% REFERENCES:
%   [1] https://www.rapidtables.com/convert/color/rgb-to-hex.html
%
%--------------------------------------------------------------------------
%
% -------
% OUTPUT:
% -------
%   pp  - (struct) structure storing plot parameters
%
%==========================================================================
function pp = PLOT_PARAMETERS

    % plot positions [x,y,l,w]
    pp.Position = [540,300,700,500];
    pp.TwoSubplotPosition = [270,300,1200,500];
    pp.ThreeSubplotPosition = [0,300,2000,500];
    pp.FourSubplotPosition = [270,100,1200,800];
    pp.TallSubplotPosition = [540,100,700,800];

    % line width
    pp.LineWidth = 1.5;
    
    % marker size
    pp.MarkerSize = 8;

    % font sizes for standard plot
    pp.FontSize_axis = 18;
    pp.FontSize_legend = 14;
    pp.FontSize_title = 18;

    % larger font sizes
    pp.FontSize_axis_big = 24;
    pp.FontSize_legend_big = 18;
    pp.FontSize_title_big = 24;

    % basic colors
    pp.black = [0,0,0];
    pp.dark_gray = [0.25,0.25,0.25];
    pp.gray = [0.5,0.5,0.5];
    pp.light_gray = [0.75,0.75,0.75];
    pp.white = [1,1,1];
    pp.red = [1,0,0];
    pp.lime = [0,1,0];
    pp.blue = [0,0,1];
    pp.yellow = [1,1,0];
    pp.cyan = [0,1,1];
    pp.magenta = [1,0,1];
    pp.maroon = [0.5,0,0];
    pp.olive = [0.5,0.5,0];
    pp.green = [0,0.5,0];
    pp.purple = [0.5,0,0.5];
    pp.teal = [0,0.5,0.5];
    pp.navy = [0,0,0.5];

    % default MATLAB colors [rgb]
    pp.matlab_blue = [0,0.4470,0.7410];
    pp.matlab_red = [0.8500,0.3250,0.0980];
    pp.matlab_yellow = [0.9290,0.6940,0.1250];
    pp.matlab_purple = [0.4940,0.1840,0.5560];
    pp.matlab_green = [0.4660,0.6740,0.1880];
    pp.matlab_cyan = [0.3010,0.7450,0.9330];
    pp.matlab_maroon = [0.6350,0.0780,0.1840];

    % lighter default MATLAB colors [rgb]
    pp.matlab_light_blue = (1-pp.matlab_blue)*0.75+pp.matlab_blue;
    pp.matlab_light_red = (1-pp.matlab_red)*0.75+pp.matlab_red;
    pp.matlab_light_yellow = (1-pp.matlab_yellow)*0.75+pp.matlab_yellow;
    pp.matlab_light_purple = (1-pp.matlab_purple)*0.75+pp.matlab_purple;
    pp.matlab_light_green = (1-pp.matlab_green)*0.75+pp.matlab_green;
    pp.matlab_light_cyan = (1-pp.matlab_cyan)*0.75+pp.matlab_cyan;
    pp.matlab_light_maroon = (1-pp.matlab_maroon)*0.75+pp.matlab_maroon;

    % Stanford colors [rgb]
    pp.cardinal_red = [140,21,21]/255;
    pp.light_cardinal_red = (1-pp.cardinal_red)*0.75+pp.cardinal_red;

end