`timescale 1ns / 1ps


module counterTest();

   reg clk;
   reg reset;
   reg enable;
   wire [2:0] out;
   integer j; 
   
   counter1 DUT(.clk(clk), .reset(reset), .enable(enable), .out(out));
   
   always #2 clk = !clk;
   
   initial
        begin
        clk = 0; reset = 0; enable = 0;
        for(j=0; j<8; j=j+1)
        begin
            #4 enable = j;
        end
        clk = 0; reset = 0; enable = 1;
        for(j=0; j<8; j=j+1)
        begin 
            #4 enable = j;
        end
        clk = 0; reset = 1; enable = 0;
        for(j=0; j<8; j=j+1)
        begin
            #4 enable = j;
        end
        clk = 0; reset = 1; enable = 1;
        for(j=0; j<8; j=j+1)
        begin
            #4 enable = j;
        end
     end
endmodule
