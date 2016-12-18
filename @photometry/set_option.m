function set_option(OBJ,CATEGORY,OPTION,OPTION_VAL);
% Makes sure the user sets options correctly
%
%
%

for i=1:length(OBJ(i))
	if isprop(OBJ(i).options,CATEGORY)
		if isfield(OBJ(i).options.(CATEGORY),OPTION)
			use_class=class(OBJ(i).options.(CATEGORY).(OPTION));
			if strcmp(use_class,class(OPTION_VAL))
				OBJ(i).options.(CATEGORY).(OPTION)=OPTION_VAL;
			else
				fprintf('Option %s must be %s\n',OPTION,use_class);
			end
		else
			fprintf('No option %s in category %s\n',OPTION,CATEGORY);
		end
	else
		fprintf('No category %s\n',CATEGORY);
	end
end
