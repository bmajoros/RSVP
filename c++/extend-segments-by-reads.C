/****************************************************************
 extend-segments-by-reads.C
 Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/WigBinary.H"
#include "BOOM/GffReader.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  void getStats(GffFeature &segment,int begin,int end,WigBinary &wig,
		int &minReads,int &maxReads,float &aveReads);
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
  if(cmd.numArgs()!=3)
    throw String("extend-segments-by-reads <segments.gff> <pileup-dir> <package-dir>");
  const String infile=cmd.arg(0);
  const String pileupDir=cmd.arg(1);
  const String packageDir=cmd.arg(2);

  GffReader reader(infile);
  Map<String,Vector<GffFeature*> > &hash=*reader.featuresBySubstrate();
  Set<String> substrates;
  hash.getKeys(substrates);
  for(Set<String>::iterator cur=substrates.begin(), end=substrates.end() ;
      cur!=end ; ++cur) {
    String substrate=*cur;
    cerr<<"processing "<<substrate<<"..."<<endl;
    String pileupFile=pileupDir+"/"+substrate+".pileup";
    if(!File::exists(pileupFile)) {
      cerr<<pileupFile<<" : not found"<<endl;
      cerr<<"SKIPPING "<<substrate<<" - no pileup file found"<<endl; 
      continue; }
    WigBinary wig(pileupFile);
    const int chrLen=wig.getLength();
    Vector<GffFeature*> &segments=hash[substrate];
    int numSegments=segments.size();
    for(int i=0 ; i<numSegments ; ++i) {
      GffFeature &segment=*segments[i];
      int begin=segment.getBegin(), end=segment.getEnd();
      int minReads, maxReads; float aveReads;
      getStats(segment,begin,end,wig,minReads,maxReads,aveReads);
      while(begin>=0 && wig.read(begin)>0) --begin;
      while(end<chrLen && wig.read(end)>0) ++end;
      segment.setEnd(end);
      segment.setBegin(begin);
      Vector<String> &extra=segment.getExtraFields();
      extra.push_back(String("minReads=")+minReads);
      extra.push_back(String("maxReads=")+maxReads);
      extra.push_back(String("aveReads=")+aveReads);
      cout<<segment.toGff();
      /*
      if(begin<segment.getBegin()) {
	int extend=segment.getBegin()-begin;
	cout<<"extended begin by "<<extend<<endl;
      }
      if(segment.getEnd()<end) {
	int extend=end-segment.getEnd();
	cout<<"extended end by "<<extend<<endl;
      }
      */
    }
  }
  return 0;
}



void Application::getStats(GffFeature &segment,int begin,int end,
			   WigBinary &wig,int &minReads,int &maxReads,
			   float &aveReads)
{
  int N=0, sum=0;
  minReads=maxReads=wig.read(begin);
  for(int i=begin ; i<end ; ++i) {
    int reads=wig.read(i);
    if(reads<minReads) minReads=reads;
    else if(reads>maxReads) maxReads=reads;
    sum+=reads;
    ++N;
  }
  aveReads=float(sum)/N;
}



