function han = m_cline(varargin)

% Draw a colored line on the map with the color dependenent on a third variable
%
% M_CLINE Plot objects on an M_MAP plot.   All of the normal
% cline options are available.  NOTE - this isn't exactly
% like cline as only the first two arguments are actually converted
% to map coords, i.e.:
%
% USAGE: M_CLINE(LON,LAT,[OPTIONS]) 
%
% Greg Pelletier 7/4/2024
%
% This software is provided "as is" without warranty of any kind.
%
% - - -
% documentation from cline.m is as follows:
%
% Draw a color-coded line by using the edge of a patch with no facecolor
%
% SYNTAX
% ======
% h = cline(x, y [, z, cdata])
%
% INPUT
% =====
% x                     vector with x-values
% y                     vector with y-values
% z (opt.)              vector with z-values
% cdata (opt.)          vector with color-data
%
% 2 input arguments =>  cdata = y; z=0      % s. Example 1
% 3 input arguments =>  cdata = z           % s. Example 2
% 4 i.a. & z = []   =>  cdata = y; z=0      % s. Example 4
%
% OUPUT
% =====
% h                 Handle to line (i.e. patch-object !!!)
%
%
% Author & Version of cline.m
% ================
% S. Hölz, TU-Berlin, seppel_mit_ppATweb.de
% V 1.0, 16.4.2007
% Created using Matlab 7.0.4 (SP2)
%
%
% Info
% ====
% This function uses the edges of a patch to represent the colored 2D/3D-line. The marker-related
% properties (i.e. 'maker','markersize','markeredgecolor','markerfacecolor') can be used as with a
% regular line. 
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% The line-related properties (i.e. 'linestyle','linewidth') WILL HAVE NO EFFECT 
% while displaying the line on screen, but will change the output when printing to file !!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!

global MAP_PROJECTION MAP_VAR_LIST

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

if nargin < 2
  help m_cline
  return
end

[x,y] = m_ll2xy(varargin{1},varargin{2});
varargin = varargin(:);
s = size(varargin,1);
h=cline(x,y,varargin{3:s});

if nargout == 1
  han = h;
end

return
