function [DATA]=kinect_read_legacy_fit(H5FILE)
%
%
%
%
%
%


tmp=hdf5info(H5FILE);
var_names={tmp.GroupHierarchy.Groups(:).Name};

DATA=struct();

for i=1:length(var_names)
	fprintf('Reading %s\n',var_names{i});
	datasets=tmp.GroupHierarchy.Groups(i).Datasets;

	for j=1:length(datasets)
		tokens=regexp(datasets(j).Name,'/','split');

		try
			DATA.(tokens{2}).(tokens{3})=h5read(H5FILE,datasets(j).Name);
		catch
			fprintf('Couldn''t read %s\n',datasets(j).Name);
			continue;
		end
	end
end
