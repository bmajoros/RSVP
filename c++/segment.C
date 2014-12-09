/****************************************************************
 segment.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "DiskMatrix.H"
#include "BOOM/Constants.H"
#include "BOOM/CommandLine.H"
#include "BOOM/WigBinary.H"
#include "genezilla/RnaJunctions.H"

class Application {
public:
  Application();
  void main(int argc,char *argv[]);
protected:
};

int main(int argc,char *argv[]) {
  try {	Application app;  app.main(argc,argv); }
  catch(const char *p) { cerr << p << endl; return -1; }
  catch(const string &msg) { cerr << msg.c_str() << endl; return -1; }
  catch(const exception &e) {
    cerr << "STL exception caught in main:\n" << e.what() << endl;
    return -1;
  }
  catch(...) {
    cerr << "Unknown exception caught in main" << endl;
    return -1;
  }
  return 0;
}



Application::Application()
{
}



void Application::main(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw BOOM::String("\n\
segment <in.pileup> <in.junctions> <threshold>\n\
");
  String pileupFile=cmd.arg(0);
  String junctionFile=cmd.arg(1);
  const float threshold=cmd.arg(2).asFloat();

  // Open pileup and junction files
  WigBinary wig(pileupFile);
  RnaJunctions junctions(junctionFile);

  // Create a vector of binary indicators for coveredness
  const char COVERED=(char)1;
  const char NOT_COVERED=(char)0;
  const int L=wig.getLength();
  File covered(tmpfile());
  for(int i=0 ; i<L ; ++i) {
    float x=wig.read(i);
    char y=x>=threshold ? COVERED : NOT_COVERED;
    covered.write(y);
  }
  int N=getNumJunctions() const;
  for(int i=0 ; i<N ; ++i) {
    const RnaJunction &junction=junctions[i];
    if(junction.getDepth()<threshold) continue;
    int begin=junction.getBegin(), end=junction.getEnd();
    covered.seek(begin);
    for(int j=begin ; j<end ; ++j) covered.write(COVERED);
  }

  // Perform segmentation



  //FILE *f=tmpfile();
  //DiskMatrix M(f,L,3);

  //  int getLength(); // in residues, not bytes
  //virtual float read(int pos); // pos is in residue units (not bytes!)

  //int getNumJunctions() const;
  //const RnaJunction &operator[](int) const;

}


