function saveNii(nii, filename, tmpFolder, untouch)
% saveNii write a nifti or gzipped nifti files
%   saveNii(nii, filename) write nii nifti structure to filename file. 
%   filename should have a .nii or .gz extension. Uses the system's tempdir to temporarily
%   nifti files if gz filename given
%
%   nii = saveNii(nii, filename, tmpFolder) uses tmpFolder as the temporary folder needed
%
%   nii = saveNii(nii, filename, tmpFolder, untouch) determines whether to use the save_untouch_nii
%   version of loader from the NIFTI library. untouch is a logical, defaulting to false;
%
%   Requires:
%   save_nii as part of NIFTI library by Jimmy Shen.
% 
%   Example:
%   nii = make_nii(zeros(10, 10, 10));
%   saveNii(nii, 'niftyfile.nii.gz'); 
%
%   Contact: Adrian V. Dalca, www.mit.edu/~adalca
%   Last Update: December, 2013.

    if nargin < 3
        untouch = false;
    end
    
    if untouch
        save_fcn = @save_untouch_nii;
    else
        save_fcn = @save_nii;
    end

    [pathstr, name, ext] = fileparts(filename);
    assert(strcmp(ext, '.nii') || strcmp(ext, '.gz'), ...
        sprintf('unknown extension <%s>. Should be .nii or .gz', ext));

    % if tmpFolder is not provided, use the system's temporary folder
    giventmp = exist('tmpFolder', 'var');
    if ~giventmp
        % note just using tempdir can fail when doing operations in parallel using the same file.
        % thus we make a new folder with tempname; 
        tmpFolder = tempname; 
    end
    
    fs = filesep;
    if strcmp(tmpFolder(end), fs)
    
        % Strip the last filesep to be consistent with fileparts behavior
        tmpFolder = tmpFolder(1:end-1); 
    end
    
    % make folder if necessary
    if ~exist(tmpFolder, 'dir')
        mkdir(tmpFolder);
    end
    
    % if the filename is a nifti, just save it
    if strcmp(ext, '.nii')
        try
            save_fcn(nii, filename);
        catch e
            warning('Caught error: "%s". Using save_untouch_nii', e.identifier);
            save_untouch_nii(nii, filename);
        end

    % else, it should be a gz extension. In this case, extract the nifti name first and save it as
    % a nifti, then gzip and erase the original file. Make sure nifti file doesn't get
    % overwritten
    else
        % set the nifti filename
        niiname = [tmpFolder, filesep, name];
        [~, ~, newext] = fileparts(niiname);
        assert(strcmp(newext, '.nii'), 'GZ filename must have a nifti pre-extention (file.nii.gz)');
        
        % if the output dir is not the temp dir, don't overwrite files. 
        if giventmp && ~strcmp(strtrimchr(tmpFolder, filesep), strtrimchr(tempdir, filesep))
            msg = sprintf('%s exists as a file or folder', niiname);
            assert(~(exist(niiname, 'file') > 1), msg);
        end

        % save the nifti
        try
            save_fcn(nii, niiname);
        catch e
            warning('Caught error: "%s". Using save_untouch_nii', e.identifier);
            if isempty(e.identifier)
                warning(e.message);
            end
            save_untouch_nii(nii, niiname);
        end

        % gzip the file (which creates a .gz file) and delete the .nii file
        if strcmp(pathstr, '')
            pathstr = pwd;
        end
        sys.fastgzip(niiname, pathstr);
        delete(niiname);
        if ~giventmp
            rmdir(tmpFolder);
        end
    end
