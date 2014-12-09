/****************************************************************
 dump-junctions.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Pipe.H"
#include "genezilla/RnaJunction.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  int main(int argc,char *argv[]);
};

int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=1)
    throw String("\n\
dump-junctions <in.junctions>\n\
");
  String infile=cmd.arg(0);
  File in(infile);
  const int recSize=RnaJunction::getSize();
  const int N=in.getSize()/recSize;
  RnaJunction J;
  for(int i=0 ; i<N ; ++i) {
    J.read(in);
    cout<<J<<endl;
  }
  return 0;
}


