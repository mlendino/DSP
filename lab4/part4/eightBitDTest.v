`timescale 1ns / 1ps

module eightBitDTest();

reg D;
reg Clk;
reg clear;
wire [7:0] Q;
integer j;

EightBitD DUT(.D(D), .Clk(Clk), .clear(clear), .Q(Q));

always #2 Clk = !Clk;

initial
    begin
    Clk = 0; D = 0; clear = 0;
    for(j=0; j<10; j=j+1)
    begin
        #4;
    end
    Clk = 0; D = 0; clear = 1;
    for(j=0; j<10; j=j+1)
    begin
        #4;
    end
    Clk = 0; D = 1; clear = 0;
    for(j=0; j<256; j=j+1)
    begin
        #4;
    end
    Clk = 0; D = 1; clear = 1;
    for(j=0; j<256; j=j+1)
    begin
        #4;
    end  
end
endmodule
