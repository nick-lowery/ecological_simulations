%Image Stack Saver
%Tristan Ursell, June 2014
%
% image_save(Im1,basename)
% image_save(Im1,basename,fmax)
% image_save(Im1,basename,fmax,'jpeg')
%
% For use with appended image stacks, e.g. most commonly stacked TIF files.
% Some operating systems throw an error indicating that the file is
% inaccessible at the time of write, which leads to a code fault and can be
% very frustrating.  This script fixes that issue, mainly on the Windows
% operating system.
%
% Im1 = the current image (matrix) you want to add to the stack
% 
% basename = is a string that specifices the file name, and potentially the
% path of the stack to save to.  Best practice is the put '.tif' at the end
% of the file name.  
%
% fmax = maximum number of allowed failures, i.e. if the script attempts to
% write to the file 'fmax' times and fails, it will give up.  This prevents
% entering an infinite loop.  The default value is fmax = 10.
%
% Any of the 'imwrite' parameters can be modified below on line 38.
%

function image_save(Im1,basename,varargin)

if ~isempty(varargin)
    fmax=varargin{1};
    if length(varargin)==2
        if strcmp(varargin{2},'jpeg')
            jpg_q=1;
        end
    else
        jpg_q=0;
    end
else
    fmax=10;
    jpg_q=0;
end

er1=0;
f=0;
while and(er1==0,f<fmax)
    try
        if jpg_q
            imwrite(Im1,basename,'Quality',100)
        else
            imwrite(Im1,basename,'writemode','append','compression','none');
        end
        er1=1;
    catch
        f=f+1;
        pause(0.1*rand)
    end
end

if f==fmax
    error(['File: ' basename ', remained inaccessible after ' num2str(fmax) ' attempts.'])
end


