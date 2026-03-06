current_gcp = gcp('nocreate');
if isempty(current_gcp)
	parpool('local', 4); % 启动一个包含4个工作者的本地并行池
else
	disp('Parallel pool already exists.');
end
