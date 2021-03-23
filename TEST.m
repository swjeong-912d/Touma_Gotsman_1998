addpath(genpath('./utils'));
mesh_dir = './meshes/';
% %% add only, simple
% [vs fs] = read_off([mesh_dir,'tetrahedron1.off']);
% [success,print_list,vs_e,fs_e] = connectivity_encoding(fs,vs,true);
% [vs_d,fs_d] = connectivity_decoding(print_list,true,vs_e,fs_e);

%% add and split, paper example
[vs fs] = read_off([mesh_dir,'tg_example.off']);
vs = vs';
fs = fs';

[success,print_list,vs_e,fs_e] = connectivity_encoding(fs,vs,true);
[vs_d,fs_d] = connectivity_decoding(print_list,true,vs_e,fs_e);
% 
% %% add and split, hard case
% [vs fs] = read_off([mesh_dir,'tetrahedron2.off']);
% [success,print_list,vs_e, fs_e] = connectivity_encoding(fs,vs,true);
% [vs_d,fs_d] = connectivity_decoding(print_list,true,vs_e,fs_e);

%% add,split, and merge, torus case
% [vs fs] = read_off([mesh_dir,'tri_torus.off']);
% [success,print_list,vs_e,fs_e] = connectivity_encoding(fs,vs,true);
% [vs_d,fs_d] = connectivity_decoding(print_list,true,vs_e,fs_e);

