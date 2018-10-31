# Assign inputs/outputs to actual pins on the FPGA
# need one pin per bit, so have to assign it accordingly
set_property PACKAGE_PIN F22 [get_ports clk]
set_property CLOCK_DEDICATED_ROUTE FALSE [get_nets {clk}];
set_property PACKAGE_PIN G22 [get_ports reset] 
set_property PACKAGE_PIN H22 [get_ports enable]
#make sure you assign both bits to different pins otherwise you'll get an error!
set_property PACKAGE_PIN T22 [get_ports out[0]]
set_property PACKAGE_PIN T21 [get_ports out[1]]  
set_property PACKAGE_PIN U21 [get_ports out[2]]  
 
# Define voltage levels (3.3 for LEDs and 1.8 for Switches)
set_property IOSTANDARD LVCMOS18 [get_ports enable] 
set_property IOSTANDARD LVCMOS18 [get_ports reset] 
set_property IOSTANDARD LVCMOS18 [get_ports clk] 
set_property IOSTANDARD LVCMOS33 [get_ports out[0]]
set_property IOSTANDARD LVCMOS33 [get_ports out[1]] 
set_property IOSTANDARD LVCMOS33 [get_ports out[2]] 
