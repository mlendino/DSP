`timescale 1ns / 1ps

module EightBitT(
    Clk,
    reset_n,
    T,
    Q);
    
     input Clk;
     input reset_n;
     input T;
     output [7:0] Q;
     
     tff t(.Clk(Clk), .reset_n(reset_n), .T(T), .Q(Q[0]));
     tff t1(.Clk(Clk), .reset_n(reset_n), .T(T && Q[0]), .Q(Q[1]));
     tff t2(.Clk(Clk), .reset_n(reset_n), .T(T && Q[0] && Q[1]), .Q(Q[2]));
     tff t3(.Clk(Clk), .reset_n(reset_n), .T(T && Q[0] && Q[1] && Q[2]), .Q(Q[3]));
     tff t4(.Clk(Clk), .reset_n(reset_n), .T(T && Q[0] && Q[1] && Q[2] && Q[3]), .Q(Q[4]));
     tff t5(.Clk(Clk), .reset_n(reset_n), .T(T && Q[0] && Q[1] && Q[2] && Q[3] && Q[4]), .Q(Q[5]));
     tff t6(.Clk(Clk), .reset_n(reset_n), .T(T && Q[0] && Q[1] && Q[2] && Q[3] && Q[4] && Q[5]), .Q(Q[6]));
     tff t7(.Clk(Clk), .reset_n(reset_n), .T(T && Q[0] && Q[1] && Q[2] && Q[3] && Q[4] && Q[5] && Q[6]), .Q(Q[7]));
     
    
    
endmodule
