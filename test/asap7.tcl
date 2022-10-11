# 
# Examples for Non Timing-driven RePlAce with TCL usage
#

set design WS16FP8
set lib_dir ./library/asap7/
set bench_dir ./library/
set design_dir ./


replace_external rep

# Import LEF/DEF files

rep set_isPlot 0
rep import_lef ${lib_dir}/asap7_tech_4x_170803.lef
rep import_lef ${lib_dir}/asap7sc7p5t_24_L_4x_170912.lef
rep import_lef ${lib_dir}/asap7sc7p5t_24_R_4x_170912.lef
rep import_lef ${lib_dir}/asap7sc7p5t_24_SL_4x_170912.lef
rep import_lef ${lib_dir}/asap7sc7p5t_24_SRAM_4x_170912.lef

rep import_def ${design_dir}/${design}.def



rep set_output ./output/

rep set_verbose_level 0

# Initialize RePlAce
rep init_replace

# place_cell with BiCGSTAB 
#rep place_cell_init_place


# print out instances' x/y coordinates
#rep print_instances

# place_cell with Nesterov method
#rep place_cell_nesterov_place

# print out instances' x/y coordinates
#rep print_instances

# Export DEF file
#rep export_def ./${design}_ASAP7_gp.def
#puts "Final HPWL: [rep get_hpwl]"

