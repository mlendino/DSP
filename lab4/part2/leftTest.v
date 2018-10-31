`timescale 1ns / 1ps


module leftTest();

    reg Clk;
    reg ShiftIn;
    reg [3:0] ParallelIn;
    reg load;
    reg ShiftEn;
    wire ShiftOut;
    wire [3:0] RegContent;
    integer j;
    
    leftShift DUT(.Clk(Clk), .ShiftIn(ShiftIn), .ParallelIn(ParallelIn), .load(load), .ShiftEn(ShiftEn), .ShiftOut(ShiftOut), .RegContent(RegContent));
    
    always #2 Clk = !Clk;
    
    initial
        begin
        ShiftEn = 0; load = 0; Clk = 0;
        for(j=0; j<16; j=j+1)
        begin
            #4 ParallelIn = j;
            ShiftIn = j%2;       
        end
        ShiftEn = 0; load = 1; Clk = 0;
        for(j=0; j<16; j=j+1)
        begin
            #4 ParallelIn = j;
            ShiftIn = j%2;       
        end
        ShiftEn = 1; load = 0; Clk = 0;
        for(j=0; j<16; j=j+1)
        begin
            #4 ParallelIn = j;
             ShiftIn = j%2;       
        end
        ShiftEn = 1; load = 1; Clk = 0;
        for(j=0; j<16; j=j+1)
        begin
            #4 ParallelIn = j;
            ShiftIn = j%2;       
        end
     end
endmodule
