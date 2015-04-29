#! /usr/bin/env python

from tempfile import mkdtemp
import subprocess
from os import chdir,getcwd,listdir,path
from shutil import rmtree,copy2

orig_dir = getcwd()
 
## Syncying small_ntuple files
temp_dir = mkdtemp()
chdir(temp_dir)
subprocess.call(['git','clone','git@github.com:manuelfs/susy_cfa'])

susy_cfa = temp_dir+'/susy_cfa'

copy2(susy_cfa+'/src/generate_small_tree.cxx',orig_dir+'/src')
copy2(susy_cfa+'/src/timer.cpp',orig_dir+'/src')
copy2(susy_cfa+'/src/utilities.cpp',orig_dir+'/src')
copy2(susy_cfa+'/inc/generate_small_tree.hpp',orig_dir+'/inc')
copy2(susy_cfa+'/inc/timer.hpp',orig_dir+'/inc')
copy2(susy_cfa+'/inc/utilities.hpp',orig_dir+'/inc')

cfg_path = susy_cfa+'/txt/small_tree_cfg'
cfg_files = [ file for file in listdir(cfg_path) if path.isfile(cfg_path+'/'+file) ]

for file in cfg_files:
    src = susy_cfa+'/txt/small_tree_cfg/'+file
    dst = orig_dir+'/txt/small_tree_cfg'
    copy2(src,dst)

chdir(orig_dir)
rmtree (temp_dir)

## Syncying plotting files
temp_dir = mkdtemp()
chdir(temp_dir)
subprocess.call(['git','clone','git@github.com:ald77/ra4_macros'])

ra4_macros = temp_dir+'/ra4_macros'

copy2(ra4_macros+'/txt/plot_styles.txt',orig_dir+'/txt')
copy2(ra4_macros+'/src/styles.cpp',orig_dir+'/src')
copy2(ra4_macros+'/src/utilities_macros.cpp',orig_dir+'/src')
copy2(ra4_macros+'/inc/utilities_macros.hpp',orig_dir+'/inc')
copy2(ra4_macros+'/inc/styles.hpp',orig_dir+'/inc')

chdir(orig_dir)
rmtree (temp_dir)
