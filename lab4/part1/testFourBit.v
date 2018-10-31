`timescale 1ns / 1ps

module testFourBit();
   //testbench contains regs that feed into the input of the device under test 
   reg [3:0] D;
   reg Clk;
   reg reset;
   reg load;
   integer i;
   //device under test outputs wires
   wire[3:0]Q;
   
   FourBitReg DUT(.D(D),.Clk(Clk),.load(load),.reset(reset),.Q(Q));
   
    always #5 Clk = !Clk;
    
    initial
        begin
        Clk = 0;
        //when reset is 1 it should set everything to 0
        reset = 0;
        D = 0;
        //load is 1 to set Q to D
        load = 1; 
        //then we need to set out input, 4 bit shift register means we need to go from 0 to 15
        //and at the output we should see that Q matches D
        for(i = 0; i<16; i = i + 1)
            #5 D = i;
        #10 
        Clk = 0;
        reset = 1;
        D = 0;
        load = 1;
        //this next for loop verifies that that reset should take precendence and that Q should always be 0
        for(i=0;i<16;i=i+1)
            #5 D = i;
        #10;
     end
     
endmodule
