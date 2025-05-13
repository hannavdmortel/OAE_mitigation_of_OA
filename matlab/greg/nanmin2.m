function [y,idx] = nanmin(a,dim,b)
% FORMAT: [Y,IDX] = NANMIN(A,DIM,[B])
% 
%    Minimum ignoring NaNs
%
%    This function enhances the functionality of NANMIN as distributed in
%    the MATLAB Statistics Toolbox and is meant as a replacement (hence the
%    identical name).
%
%    If fact NANMIN simply rearranges the input arguments to MIN because
%    MIN already ignores NaNs.
%
%    NANMIN(A,DIM) calculates the minimum of A along the dimension DIM of
%    the N-D array X. If DIM is omitted NANMIN calculates the minimum along
%    the first non-singleton dimension of X.
%
%    NANMIN(A,[],B) returns the minimum of the N-D arrays A and B.  A and
%    B must be of the same size.
%
%    Comparing two matrices in a particular dimension is not supported,
%    e.g. NANMIN(A,2,B) is invalid.
%    
%    [Y,IDX] = NANMIN(X,DIM) returns the index to the minimum in IDX.
%    
%    Similar replacements exist for NANMAX, NANMEAN, NANSTD, NANMEDIAN and
%    NANSUM which are all part of the NaN-suite.
%
%    See also MIN

% -------------------------------------------------------------------------
%    author:      Jan Gläscher
%    affiliation: Neuroimage Nord, University of Hamburg, Germany
%    email:       glaescher@uke.uni-hamburg.de
%    
%    $Revision: 1.1 $ $Date: 2004/07/15 22:42:14 $

% Copyright (c) 2004, Jan Gläscher
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% * Redistributions of source code must retain the above copyright notice, this
  % list of conditions and the following disclaimer.

% * Redistributions in binary form must reproduce the above copyright notice,
  % this list of conditions and the following disclaimer in the documentation
  % and/or other materials provided with the distribution

% * Neither the name of Universitätsklinikum Hamburg-Eppendorf nor the names of its
  % contributors may be used to endorse or promote products derived from this
  % software without specific prior written permission.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if nargin < 1
	error('Requires at least one input argument')
end

if nargin == 1
	if nargout > 1
		[y,idx] = min(a);
	else
		y = min(a);
	end
elseif nargin == 2
	if nargout > 1
		[y,idx] = min(a,[],dim);
	else
		y = min(a,[],dim);
	end
elseif nargin == 3
	if ~isempty(dim)
		error('Comparing two matrices along a particular dimension is not supported')
	else
		if nargout > 1
			[y,idx] = min(a,b);
		else
			y = min(a,b);
		end
	end
elseif nargin > 3
	error('Too many input arguments.')
end

% $Id: nanmin.m,v 1.1 2004/07/15 22:42:14 glaescher Exp glaescher $
