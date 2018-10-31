`timescale 1ns / 1ps

module counter1(
    clk,
    reset,
    enable,
    out);
    
    input clk;
    input reset;
    input enable;
    output [2:0] out;
    reg [2:0] current_state;
    reg [2:0] next_state;
    localparam  s0 = 3'b000, s1 = 3'b001, s2 = 3'b011, s3 = 3'b101, s4 = 3'b111, s5 = 3'b010;
    //just guarantees it changes on the positive edge of the clock
    assign out=current_state;
    always@(posedge clk)
        begin
            if(reset)
                current_state = s0;
            if(enable && ~reset)
                current_state = next_state;
        end
        
    always@(*)
        begin
            next_state = current_state;
            case(current_state)
            s0: if(enable)
                    next_state = s1;
            s1: if(enable)
                    next_state = s2;
            s2: if(enable)
                    next_state = s3;
            s3: if(enable)
                    next_state = s4;
            s4: if(enable)
                    next_state = s5;
            s5: if(enable)
                    next_state = s0;
            default:
                    next_state = s0;
            endcase 
        end                  
endmodule
