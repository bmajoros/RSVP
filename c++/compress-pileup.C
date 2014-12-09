/****************************************************************
 compress-pileup.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Pipe.H"
using namespace std;
using namespace BOOM;


// first, run:
// samtools mpileup sorted.bam | simplify-pileup.pl | gzip > pileup.txt.gz

class Application
{
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
compress-pileup <pileup.txt.gz> <outdir>\n\
    where pileup.txt.gz is the gzipped output of simplify-pileup.pl\n\
");
  String infile=cmd.arg(0), outDir=cmd.arg(1);
  File::mkdir(outDir);
  String currentSubstrate;
  Pipe in(String("cat ")+infile+" | gunzip","r");
  File *out=NULL;
  long currentPos=0;
  const float f0=0;
  while(!in.eof()) {
    String line=in.getline();
    line.trimWhitespace();
    BOOM::Vector<BOOM::String> &fields=*line.getFields();
    int numFields=fields.size();
    if(numFields==0) { delete &fields; continue; }
    if(numFields!=3) 
      throw String("wrong number of fields in input record: ")+line;
    const String substrate=fields[0];
    const int pos=fields[1].asInt(); // ### no "-1" here
    float depth=fields[2].asFloat();
    delete &fields;
    if(substrate!=currentSubstrate) {
      if(out) delete out;
      String outfile=outDir+"/"+substrate+".pileup";
      out=new File(outfile,"w");
      currentSubstrate=substrate;
      currentPos=0;
    }
    if(pos!=currentPos) //{ out->seek(pos); currentPos=pos; }
      for(; currentPos<pos ; ++currentPos) out->write(f0);
    out->write(depth);
    ++currentPos;
  }
  in.close();
  if(out) out->close();
  
  return 0;
}

