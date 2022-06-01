# 
# Examples for Non Timing-driven RePlAce with TCL usage
#

set design gcd
set lib_dir ./library/nangate45/
set bench_dir ./library/
set design_dir ./design/nangate45/${design}
set num   1
set fname ispd18_test${num}

replace_external rep

# Import LEF/DEF files
rep set_fastWL 1
rep set_isPlot 0
#rep import_lef ${lib_dir}/NangateOpenCellLibrary.lef
#rep import_def ${design_dir}/${design}.def

rep import_lef ${bench_dir}/${fname}/${fname}.input.lef
rep import_def ${bench_dir}/${fname}/${fname}.input.def

rep set_output ./output/

rep set_verbose_level 0

# Initialize RePlAce
rep init_replace

# place_cell with BiCGSTAB 
rep place_cell_init_place


# print out instances' x/y coordinates
#rep print_instances

# place_cell with Nesterov method
rep place_cell_nesterov_place

# print out instances' x/y coordinates
#rep print_instances

# Export DEF file
rep export_def ./${design}_nan45_nontd.def
puts "Final HPWL: [rep get_hpwl]"

