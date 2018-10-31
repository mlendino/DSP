`timescale 1ns / 1ps


module tff(
    Clk,
    reset_n,
    T,
    Q);
    
    input Clk;
    input reset_n;
    input T;
    output reg Q = 0;
    
    always @(negedge Clk)
        if(reset_n)
            Q <= 1'b0;
        else if (T)
            Q <= ~Q;
endmodule
