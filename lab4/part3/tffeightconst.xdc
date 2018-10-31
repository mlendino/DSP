# Assign inputs/outputs to actual pins on the FPGA
# need one pin per bit, so have to assign it accordingly
set_property PACKAGE_PIN F22 [get_ports Clk]
set_property CLOCK_DEDICATED_ROUTE FALSE [get_nets {Clk}];
set_property PACKAGE_PIN G22 [get_ports T] 
set_property PACKAGE_PIN H22 [get_ports reset_n]
#make sure you assign both bits to different pins otherwise you'll get an error!
set_property PACKAGE_PIN T22 [get_ports Q[0]]
set_property PACKAGE_PIN T21 [get_ports Q[1]]  
set_property PACKAGE_PIN U21 [get_ports Q[2]]  
set_property PACKAGE_PIN U22 [get_ports Q[3]]  
set_property PACKAGE_PIN V22 [get_ports Q[4]]  
set_property PACKAGE_PIN W22 [get_ports Q[5]]  
set_property PACKAGE_PIN U19 [get_ports Q[6]]  
set_property PACKAGE_PIN U14 [get_ports Q[7]]  

 
# Define voltage levels (3.3 for LEDs and 1.8 for Switches)
set_property IOSTANDARD LVCMOS18 [get_ports reset_n] 
set_property IOSTANDARD LVCMOS18 [get_ports T] 
set_property IOSTANDARD LVCMOS18 [get_ports Clk] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[0]]
set_property IOSTANDARD LVCMOS33 [get_ports Q[1]] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[2]] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[3]] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[4]] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[5]] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[6]] 
set_property IOSTANDARD LVCMOS33 [get_ports Q[7]]