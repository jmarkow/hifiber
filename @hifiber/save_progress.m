function save_progress(OBJ, fname)
	%
	%
	%
	if nargin==1
		fname = 'hifiber_object.mat';
	end
	
	savefun=@(save_path,hifiber_object) save(save_path,'hifiber_object','-v7.3');
	
	if ~exist(OBJ.options.save_dir,'dir')
		mkdir(OBJ.options.save_dir);
	end
	
	savefun(fullfile(OBJ.options.save_dir,fname),OBJ);
	
end % function
