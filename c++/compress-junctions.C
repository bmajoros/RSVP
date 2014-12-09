/****************************************************************
 compress-junctions.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Pipe.H"
#include "BOOM/VectorSorter.H"
#include "genezilla/RnaJunction.H"
using namespace std;
using namespace BOOM;


// first run:
// python BED2Intron.py tophat/junctions.bed > junctions.bed.intron

class Application {
  String outDir;
  JunctionComparator cmp;
  void dump(const String &currentSubstrate,Vector<RnaJunction> &);
public:
  Application();
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



Application::Application()
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("\n\
compress-junctions <junctions.bed> <outdir>\n\
  First run this: \n\
    python BED2Intron.py tophat/junctions.bed > junctions.bed\n\
");
  String infile=cmd.arg(0);
  outDir=cmd.arg(1);
  File::mkdir(outDir);
  String currentSubstrate;
  File in(infile);
  Vector<RnaJunction> V;
  while(!in.eof()) {
    String line=in.getline();
    line.trimWhitespace();
    BOOM::Vector<BOOM::String> &fields=*line.getFields();
    int numFields=fields.size();
    if(numFields==0) { delete &fields; continue; }
    if(numFields<5) 
      throw String("wrong number of fields in input record: ")+line;
    const String substrate=fields[0];
    const int begin=fields[1].asInt(), end=fields[2].asInt();
    const float depth=fields[3].asFloat();
    const char strand=fields[4][0];
    delete &fields;
    if(currentSubstrate.isEmpty()) currentSubstrate=substrate;
    if(substrate!=currentSubstrate) {
      dump(currentSubstrate,V);
      currentSubstrate=substrate;
      V.clear();
    }
    RnaJunction J(begin-1,end,depth,strand);
    V.push_back(J);
  }
  if(V.size()>0) dump(currentSubstrate,V);
  in.close();
  return 0;
}



void Application::dump(const String &currentSubstrate,
		       Vector<RnaJunction> &V) {
  String outfile=outDir+"/"+currentSubstrate+".junctions";
  File out(outfile,"w");
  VectorSorter<RnaJunction> sorter(V,cmp);
  sorter.sortAscendInPlace();
  Vector<RnaJunction>::iterator cur=V.begin(), end=V.end();
  for(; cur!=end ; ++cur) {
    const RnaJunction &J=*cur;
    out.write(J.getBegin());
    out.write(J.getEnd());
    out.write(J.getDepth());
    out.write(J.getStrand());
  }
}
