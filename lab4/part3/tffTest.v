`timescale 1ns / 1ps

module tffTest();

reg Clk;
reg reset_n;
reg T;
wire [7:0] Q;
integer j;

EightBitT DUT(.Clk(Clk), .reset_n(reset_n), .T(T), .Q(Q));

always #2 Clk = !Clk;

initial
    begin
    Clk = 0; reset_n = 0; T = 0;
    for(j=0; j<10; j=j+1)
    begin
        #4;
    end  
    Clk = 0; reset_n = 0; T = 1;
    for(j=0; j<260; j=j+1)
    begin
        #4;
    end
    Clk = 0; reset_n = 1; T = 0;
    for(j=0; j<10; j=j+1)
    begin
        #4;
    end
    Clk = 0; reset_n = 1; T = 1;
    for(j=0; j<260; j=j+1)
    begin
        #4;
    end
end 
endmodule
