`timescale 1ns / 1ps

//Model a 4-bit register with synchronous reset and load.
//Assign Clk, D input, reset, load, and output Q
module FourBitReg(
    D,
    Clk,
    load,
    reset,
    Q);
    
    input [3:0] D;
    input Clk;
    input reset;
    input load;
    output reg [3:0]Q;
    //Controlled register clocking with a synchronous reset and enable/load
    always @(posedge Clk)
        if(reset)
        begin
            Q <= 0;
        end else if(load)
        begin
            Q <= D;
        end 
endmodule
