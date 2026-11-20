//*************************************************************************************************
//*************************************************************************************************
//* The MIT License (MIT)                                                                         *
//* Copyright (C) 2026                                                                            *
//* Ivan Matic and Rados Radoicic                                                                 *
//*                                                                                               *
//* Permission is hereby granted, free of charge, to any person obtaining a copy of this          *
//* software and associated documentation files (the "Software"), to deal in the Software         *
//* without restriction, including without limitation the rights to use, copy, modify, merge,     *
//* publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons    *
//* to whom the Software is furnished to do so, subject to the following conditions:              *
//*                                                                                               *
//* The above copyright notice and this permission notice shall be included in all copies or      *
//* substantial portions of the Software.                                                         *
//*                                                                                               *
//* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,           *
//* INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR      *
//* PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE     *
//* FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR          *
//* OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        *
//* DEALINGS IN THE SOFTWARE.                                                                     *
//*************************************************************************************************

#ifndef _INCL_HSP_CPP
#define _INCL_HSP_CPP
namespace HSP{
std::string helpString(){
    std::string res;
    res+="Usage: ./piercing -i W14.txt\n";
    res+="\n";
    res+="CONTRACT. The input file contains the records\n";
    res+="    bUpperNum            integers bNum_ij for 0<=i<=j<=14, row i listing bNum_ii,...,bNum_i14\n";
    res+="    bDen                 a positive integer\n";
    res+="    betaNum, betaDen     positive integers\n";
    res+="    massNum, massDen     integers, massDen positive.\n";
    res+="Let b_ij = b_ji = bNum_ij/bDen and\n";
    res+="    W_14(x,y) = sum_{i,j=0}^{14} b_ij B^14_i((2x-1)^2) B^14_j((2y-1)^2).\n";
    res+="If the last line printed is\n";
    res+="    "+GL_PASS+"\n";
    res+="then (V1), (V2), (V3) of Section 4 of the paper hold for W_14 with\n";
    res+="beta = betaNum/betaDen, and Theorem 3.2 of the paper holds. Hence\n";
    res+="    liminf p_n/n >= (massNum/massDen)/(betaNum/betaDen),\n";
    res+="where p_n is the least number of lines that pierce every cell of an n x n board.\n";
    res+="No modification of the input can make the program print the line above for a\n";
    res+="false statement. Tampering only makes the proof fail.\n";
    res+="Text outside the [key]...[/key] [value]...[/value] records is ignored.\n";
    return res;
}
}
#endif
