/****************************************************************
 extract-pileup-interval.C
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


// first, run:
// samtools mpileup sorted.bam | simplify-pileup.pl | gzip > pileup.txt.gz

class Application {
  void processPileup(const String &infile,int begin,int end,
		     const String &outfile);
  void processJunctions(const String &infile,int begin,int end,
			const String &outfile);
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
  if(cmd.numArgs()!=5)
    throw String("\n\
extract-pileup-interval <*.pileup> <*.junctions> <begin> <end> <filestem-out>\n\
");
  String pileupInfile=cmd.arg(0);
  String junctionInfile=cmd.arg(1);
  int begin=cmd.arg(2).asInt();
  int end=cmd.arg(3).asInt();
  String filestem=cmd.arg(4);
  processPileup(pileupInfile,begin,end,filestem+".pileup");
  processJunctions(junctionInfile,begin,end,filestem+".junctions");
  return 0;
}



void Application::processPileup(const String &infile,int begin,int end,
				const String &outfile)
{
  File in(infile,"r");
  in.seek(begin*sizeof(float));
  File out(outfile,"w");
  for(int pos=begin ; pos<end ; ++pos) {
    float f=in.readFloat();
    out.write(f);
  }
  in.close();
  out.close();
}



void Application::processJunctions(const String &infile,int begin,int end,
				   const String &outfile)
{
  const int recSize=RnaJunction::getSize();
  File in(infile,"r"), out(outfile,"w");
  long fileSize=in.getSize();
  int N=fileSize/recSize;
  if(N==0) return;
  // Start with a binary-search like procedure to find the last
  // record having the right "begin" value
  int b=0, e=N;
  RnaJunction J;
  while(b+1<e) {
    int m=(b+e)/2;
    in.seek(m*recSize);
    J.read(in);
    if(begin<J.getBegin()) e=m;
    else b=m;
  }
  // Now if there are any records with the appropriate "begin" value,
  // they must lie immediately to the left of e
  int i;
  for(i=e-1 ; i>=0  ; --i) {
    in.seek(i*recSize);
    J.read(in);
    if(J.getBegin()<begin) break;
  }
  for(; i<N ; ++i) {
    in.seek(i*recSize);
    J.read(in);
    if(J.getBegin()>=begin && J.getEnd()<end) {
      J.shift(-begin);
      J.write(out);
    }
    if(J.getBegin()>=end) break;
  }
  in.close();
  out.close();
}



