`timescale 1ns / 1ps

module EightBitD(
    D,
    Clk,
    clear,
    Q);
    
  input D;
  input Clk;
  input clear;
  output [7:0] Q;
  
  dff d1(.D(~Q[0]), .Clk(Clk), .clear(clear), .Q(Q[0]));
  dff d2(.D(~Q[1]), .Clk(~Q[0]), .clear(clear), .Q(Q[1]));
  dff d3(.D(~Q[2]), .Clk(~Q[1]), .clear(clear), .Q(Q[2]));
  dff d4(.D(~Q[3]), .Clk(~Q[2]), .clear(clear), .Q(Q[3]));
  dff d5(.D(~Q[4]), .Clk(~Q[3]), .clear(clear), .Q(Q[4]));
  dff d6(.D(~Q[5]), .Clk(~Q[4]), .clear(clear), .Q(Q[5]));
  dff d7(.D(~Q[6]), .Clk(~Q[5]), .clear(clear), .Q(Q[6]));
  dff d8(.D(~Q[7]), .Clk(~Q[6]), .clear(clear), .Q(Q[7]));
  
endmodule
