function [] = adpath(pathname)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
addpath('subfunction_code');
example_name=['source_code/',pathname];
addpath(example_name); 
end

