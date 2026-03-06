current_gcp = gcp('nocreate');
if ~isempty(current_gcp)
    delete(gcp('nocreate')); % 关闭当前的并行池
else
    disp('No parallel pool to close.');
end