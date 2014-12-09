/****************************************************************
 normalize-wig.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/WigBinary.H"
using namespace std;
using namespace BOOM;


class Application
{
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])  {
    try  {Application app; return app.main(argc,argv); }
    catch(const char *p) { cerr << p << endl; }
    catch(const string &msg) { cerr << msg.c_str() << endl; }
    catch(const exception &e) 
      { cerr << "STL exception caught in main:\n" << e.what() << endl; }
    catch(...)
      {	cerr << "Unknown exception caught in main" << endl; }
    return -1;
  }



Application::Application()
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("normalize-wig <in-out.wig> <pseudocount>");
  String filename=cmd.arg(0);
  float pseudocount=cmd.arg(1).asFloat();

  WigBinary wig(filename,true);
  const int L=wig.getLength();
  float min, max;
  wig.getExtrema(min,max);
  max+=pseudocount;
  for(int i=0 ; i<L ; ++i) 
    wig.write(i,wig.read(i)/max);
  
  return 0;
}

