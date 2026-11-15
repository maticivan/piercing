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

#include "src/fundamental/SHF_standardHeadersMinimal.cpp"
#include <thread>
#include "src/fundamental/ssm.cpp"
#include "src/fundamental/fst.cpp"
#include "src/fundamental/MFRF_massiveFindReplace.cpp"
#include "src/fundamental/BF_bijFunction.cpp"
#include "src/fundamental/TMF_timer.cpp"
#include "src/fundamental/SF_strings.cpp"
#include "src/fundamental/GF_globalInitsMinimal.cpp"
#include "src/fundamental/RNDF_randomCodesAndPermutations.cpp"
#include "src/fundamental/IOF_io.cpp"
std::string GL_PASS="THE THEOREM WAS PROVED";
#include "src/tools/CLAT_commandLineArguments.cpp"
#include "src/piercing/BSP_bernsteinSubdivision.cpp"
#include "src/piercing/WCP_wCertificate.cpp"
#include "src/piercing/HSP_helpString.cpp"
int main(int argc, char **argv){
    TMF::Timer tm;
    CLAT::CommandLineArgs inputArgs=CLAT::getArguments(CLAT::argsInStrings(argc,argv));
    if(inputArgs.inputFile=="notFound"){
        std::cout<<HSP::helpString()<<"\n";
        return 0;
    }
    std::string sourceTxt=IOF::fileToString(inputArgs.inputFile);
    if(sourceTxt=="fileNotFound"){
        std::cout<<"Input file not found\n";
        return 0;
    }
    int success=0;
    tm.start();
    try{
        success=WCP::mainWork<Int1024>(sourceTxt);
    }
    catch(const std::overflow_error&){
        success=0;
        std::cout<<"WARNING: an integer did not fit in 1024 bits. Results are not trustworthy.\n";
    }
    tm.end();
    std::cout<<"Total time: "<<tm.getTimeSec()<<" seconds.\n";
    if(success){
        std::cout<<GL_PASS<<"\n";
    }
    return 0;
}
