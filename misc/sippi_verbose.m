% sippi_verbose : displays verbose information to the console
%
% Call:
%  sippi_verbose(txt,verbose)
%
% txt [string] : text to be displayed
% verbose [integer] (def=0) : increase to see more information
%
% 'vlevel' must be set in the sippi_verbose.m m-file.
%
% All entries with vebose>vlevel are displayed
%

%
% entries with a higher verbose value has a higher chance of being displayed
% that entries with lower verbose values
% verbose [0] : normal (default)
%         [-1] : little info 
%         [-2] : no info
%
% Default level is set using an environmental variable
% e.g. setenv('SIPPI_VERBOSE_LEVEL','1')

% 


%
function varargout=sippi_verbose(txt,verbose)
  

  vlevel=0; % SHOW ALL VERBOSE INFO WITH 'VERBOSE' ABOVE 'VLEVEL'
  try
      % GET VERBOSE LEVEL FROM SYSTEM VARIABLE IF SET
      tmp=str2num(getenv('SIPPI_VERBOSE_LEVEL'));
      if ~isempty(tmp)
          vlevel=tmp;
      end
  end
  
  if nargout;
      varargout{1} = vlevel;    
  end
  if nargin==0;
      return
  end
  
  if nargin==1,
    verbose=0;
  end
 
  
  if (verbose<=vlevel),
  %if (verbose>=vlevel),
    txt1=mfilename;
    disp(sprintf('%s',txt));
  end