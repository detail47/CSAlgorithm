% 运行所有测试
test_files = dir(fullfile(pwd, '+Test', '*.m'));
for i = 1:length(test_files)
    name = test_files(i).name;
    names = split(name, '.');
    test_name = names{1};
    relative_name = ['Test.' test_name];
    test_func = str2func(relative_name);
    test_func();
    disp(['Test ' test_name ' passed.']);
end
