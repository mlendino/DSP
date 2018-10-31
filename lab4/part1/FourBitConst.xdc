# Assign inputs/outputs to actual pins on the FPGA
# need one pin per bit, so have to assign it accordingly
set_property PACKAGE_PIN F22 [get_ports Clk]
set_property CLOCK_DEDICATED_ROUTE FALSE [get_nets {Clk}];
set_property PACKAGE_PIN G22 [get_ports load] 
set_property PACKAGE_PIN H22 [get_ports reset]
#make sure you assign both bits to different pins otherwise you'll get an error!
set_property PACKAGE_PIN F21 [get_ports D[0]]
set_property PACKAGE_PIN H19 [get_ports D[1]]  
set_property PACKAGE_PIN H18 [get_ports D[2]]  
set_property PACKAGE_PIN H17 [get_ports D[3]]  
set_property PACKAGE_PIN V22 [get_ports Q[0]]  
set_property PACKAGE_PIN W22 [get_ports Q[1]]  
set_property PACKAGE_PIN U19 [get_ports Q[2]]  
set_property PACKAGE_PIN U14 [get_ports Q[3]]  

 
# Define voltage levels (3.3 for LEDs and 1.8 for Switches)
set_property IOSTANDARD LVCMOS18 [get_ports reset] 
set_property IOSTANDARD LVCMOS18 [get_ports load] 
set_property IOSTANDARD LVCMOS18 [get_ports Clk] 
set_property IOSTANDARD LVCMOS18 [get_ports D[0]]
set_property IOSTANDARD LVCMOS18 [get_ports D[1]] 
set_property IOSTANDARD LVCMOS18 [get_ports D[2]] 
set_property IOSTANDARD LVCMOS18 [get_ports D[3]] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[0]] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[1]] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[2]] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[3]]