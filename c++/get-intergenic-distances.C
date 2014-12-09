/****************************************************************
 get-intergenic-distances.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "BOOM/Constants.H"
#include "BOOM/CommandLine.H"
#include "BOOM/WigBinary.H"
#include "BOOM/Regex.H"
#include "BOOM/Array1D.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/GffReader.H"
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
  if(cmd.numArgs()!=2)
    throw BOOM::String("\n\
get-intergenic-distances <transcripts.gff> <pileup-dir>\n\
");
  String gffFile=cmd.arg(0);
  String dir=cmd.arg(1);

  GffReader reader(gffFile);
  Map<String,TranscriptList*> &hash=*reader.loadByContig();
  Set<String> substrates;
  hash.getKeys(substrates);
  const int numSubstrates=substrates.size();
  for(Set<String>::iterator cur=substrates.begin(), end=substrates.end() ; 
      cur!=end ; ++cur) {
    const String substrate=*cur;
    String pileupFile=dir+"/"+substrate+".pileup";
    WigBinary wig(pileupFile);
    String junctionFile=dir+"/"+substrate+".junctions";
    RnaJunctions junctions(junctionFile);
    TranscriptList &transcripts=*hash[substrate];
    const int numTranscripts=transcripts.size();
    for(int i=0 ; i<numTranscripts-1 ; ++i) {
      GffTranscript *thisTranscript=transcripts[i];
      GffTranscript *nextTranscript=transcripts[i+1];
      int begin=thisTranscript->getEnd();
      int end=nextTranscript->getBegin();
      while(begin<end && wig.read(begin)>0) ++begin;
      while(end>begin && wig.read(end)>0) --end;

      // ### what about junctions???

      int len=end-begin;
      if(len>0) cout<<len<<endl;
    }

  }

}
