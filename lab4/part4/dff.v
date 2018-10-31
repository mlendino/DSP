`timescale 1ns / 1ps

module dff(
    D,
    Clk,
    clear,
    Q);
    
    input D;
    input Clk;
    input clear;
    output reg Q;
    
    always@(posedge Clk or posedge clear)
        if (clear)
        begin
            Q <= 1'b0;
        end else
        begin
            Q <= D;
        end
endmodule
