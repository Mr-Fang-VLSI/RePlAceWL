Root_dir = `pwd`
Lib_file = "$Root_dir/cds.lib"
source /opt/coe/synopsys/lc/P-2019.03-SP2

read_lib Lib_file
write_lib cds -format db -output cds.db

quit

