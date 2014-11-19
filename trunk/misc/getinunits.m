function out = getinunits(h, prop, unts)

%GETINUNITS   Get object properties in specified units
%   V = GETINUNITS(H, PROP, UNITS) returns the object property
%   in the specified UNITS. It will leave the 'Units' and 'FontUnits'
%   property unchanged afterwards.
%
%   H is the handle of the object. If it is an M-element array of handles,
%   the function will return an M-by-1 cell array. PROP can be a string or
%   a cell array of strings. If it is a 1-by-N or N-by-1 cell array, the
%   function will return an M-by-N cell array of values. UNITS can be a
%   string or a cell array. If it is a cell array, then PROP must also be a
%   cell array with the same size as UNITS, and each cell element of UNITS
%   corresponds to a cell element of PROP.
%
%   V = GETINUNITS(H, PROP) is the same as GET(H, PROP)
%
%   Examples:
%     V = GETINUNITS(H, 'Position', 'Pixels')
%     V = GETINUNITS(H, {'FontSize', 'Position'}, 'Normalized')
%     V = GETINUNITS(H, {'FontSize', 'Position'}, {'Points', 'Pixels'})
%
%   See also GET, SET

%   VERSIONS:
%     v1.0 - first version
%     v1.1 - check to see if FontUnits and Units properties exist or not.
%            This addresses the bug that gives an error with figures.
%
% Jiro Doke
% Copyright 2005-2010 The MathWorks, Inc.

error(nargchk(2, 3, nargin));

if nargin == 2
  % use the original units
  out = get(h, prop);
  return;
  
elseif ischar(unts)
  % convert to cell arrays 
  if ischar(prop)
    prop = {prop};
    unts = {unts};
  elseif iscell(prop) && min(size(prop)) == 1
    unts = repmat({unts}, size(prop));
  else
    error('PROP must be a string or a cell array of strings');
  end
  
elseif iscell(unts) && iscell(prop)
  if ~(isequal(size(unts), size(prop)) && min(size(prop)) == 1)
    error('Cell arrays PROP and UNITS must have the same number of elements');
  end
  
else
  error('UNITS must be a string, or both PROP and UNITS must be a cell array of strings');
end


% Now, PROP and UNTS both are cell array of strings with the same size
% Process each property individually, in case each property is asked in
% different units.
for iProp = 1:length(prop)
  try
    curFontUnits = cellstr(get(h, 'FontUnits'));
  catch % FontUnits property doesn't exist
    curFontUnits = [];
  end
  try
    curUnits     = cellstr(get(h, 'Units'));
  catch % Units property doesn't exist
    curUnits     = [];
  end
  
  % attempt to change units
  try
    set(h, 'FontUnits', unts{iProp});
  end
  try
    set(h, 'Units', unts{iProp});
  end
  
  if length(h) == 1
    out{1, iProp} = get(h, prop{iProp});
  else
    out(:, iProp) = get(h, prop{iProp});
  end
  
  % restore units
  if ~isempty(curFontUnits)
    set(h, {'FontUnits'}, curFontUnits);
  end
  if ~isempty(curUnits)
    set(h, {'Units'}, curUnits);
  end
end

% if only one object and one property, extract from cell
if numel(out) == 1
  out = out{1};
end




%#ok<*AGROW>
%#ok<*TRYNC>
%#ok<*CTCH>