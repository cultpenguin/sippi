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
% The verbose level can be set either using a an environmental variable
% e.g. setenv('SIPPI_VERBOSE_LEVEL','1')
% or one can call sippi_verbose with a third argument, which will set the
% verbose level. To set the verbose level to 2, use:
% sippi_verbose('',0,2);
%
% The verbose level can also be set using an environmental variable:
%    %% VERBOSITY
%    The amount of text info displayed at the prompt, can be controlled by
%    setenv('SIPPI_VERBOSE_LEVEL','2') % all: information on chain swapping
%    setenv('SIPPI_VERBOSE_LEVEL','1') % information about seq-gibbs step update
%    setenv('SIPPI_VERBOSE_LEVEL','0'); % [def] frequent update
%    setenv('SIPPI_VERBOSE_LEVEL','-1'); % rare update om finish time
%    setenv('SIPPI_VERBOSE_LEVEL','-2'); % indication of stop and start
%    setenv('SIPPI_VERBOSE_LEVEL','-3'); % none
%
%

function varargout=sippi_verbose(txt,verbose,set_verbose_level)
  
  if nargin==3;
      setenv('SIPPI_VERBOSE_LEVEL',num2str(set_verbose_level));
      sippi_verbose(sprintf('setting verbose level at %d',set_verbose_level));
  end

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