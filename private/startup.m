set(groot, "DefaultFigureWindowStyle", "normal");
set(groot, "DefaultAxesFontSize", 16);

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','Default');
    set(groot, default_name, 'latex');
end

clc; clear; close all;
