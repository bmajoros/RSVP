/****************************************************************
 train-segmenter.C
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
#include "genezilla/RnaJunctions.H"

struct Region {
  Region(int b,int e,float m) : begin(b), end(e), meanValue(m) {}
  int begin, end;
  float meanValue;
};


class Application {
public:
  Application();
  void main(int argc,char *argv[]);
protected:
  void fillArray(Array1D<float> &,WigBinary &,RnaJunctions &);
  void getRegions(Array1D<float> &,Vector<WigInterval> &regions);
  void addCounts(WigBinary &wig,WigInterval w,Vector<float> &counts);
  void addCounts(Array1D<float> &,WigInterval w,Vector<float> &counts);
  Regex filePattern;
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
  : filePattern("([^\/]+).pileup")
{
}



void Application::main(int argc,char *argv[])
{
  // Process command line
  BOOM::CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw BOOM::String("\n\
train-segmenter <chunks-dir> <bin-size>\n\
");
  String dir=cmd.arg(0);
  const int binSize=cmd.arg(1).asInt();

  Vector<String> files;
  File::getFileList(dir,files);
  const int numBins=1000/binSize;
  Array1D< Vector<float> > hist(numBins);
  const int numFiles=files.size();
  for(int i=0 ; i<numFiles ; ++i) {
    const String &filename=files[i];
    if(filePattern.match(filename)) {
      String chunkID=filePattern[1];
      WigBinary wig(dir+"/"+filename);
      RnaJunctions junctions(dir+"/"+chunkID+".junctions");
      Vector<WigInterval> regions;
      const int L=wig.getLength();
      Array1D<float> A(L);
      fillArray(A,wig,junctions);
      getRegions(A,regions);
      //if(regions.size()<2) continue;
      if(regions.size()==0) continue;
      if(regions.size()==1) {
	Vector<float> counts;
	addCounts(A,regions[0],counts);
	SummaryStats stats(counts);
	float meanCoverage=stats.getMean();
	int bin=int(ceil(meanCoverage)/binSize);
	if(bin<numBins) { hist[bin].push_back(0); }
	continue;
      }
      int N=regions.size();
      for(int j=0 ; j<N-1 ; ++j) {
	WigInterval w1=regions[j], w2=regions[j+1];
	int len=w2.begin-w1.end;
	Vector<float> counts;
	addCounts(A,w1,counts);
	addCounts(A,w2,counts);
	SummaryStats stats(counts);
	//cout<<stats.getMean()<<"\t"<<len<<endl;
	float meanCoverage=stats.getMean();
	int bin=int(ceil(meanCoverage)/binSize);
	if(bin<numBins) { hist[bin].push_back(len); }
      }
    }
  }
  Vector<float> meanCurve, sdCurve;
  for(int i=0 ; i<numBins ; ++i) {
    if(hist[i].size()<10) continue;
    DirectComparator<float> cmp;
    VectorSorter<float> sorter(hist[i],cmp);
    sorter.sortDescendInPlace();
    int N=hist[i].size()/2;
    Vector<float> half(N);
    for(int j=0 ; j<N ; ++j) half[j]=hist[i][j];
    //SummaryStats stats(hist[i]);
    SummaryStats stats(half);
    float coverage=i*binSize;
    //cout<<coverage<<"\t"<<stats.getMean()<<"\t"<<stats.getStdDev()<<endl;
    meanCurve.push_back(stats.getMean());
    sdCurve.push_back(stats.getStdDev());
  }
  
  // Smooth curves
  const int numIterations=20;
  const int RADIUS=2;
  Vector<float> smoothMeanCurve(numBins), smoothSdCurve(numBins);
  smoothMeanCurve=meanCurve;
  smoothSdCurve=sdCurve;
  for(int i=0 ; i<numIterations ; ++i) {
    for(int i=1 ; i<numBins-1 ; ++i) {
      smoothMeanCurve[i]=0;
      int radius=RADIUS;
      if(radius>i) radius=i;
      if(i+radius>=numBins) radius=numBins-1-i;
      for(int k=i-radius ; k<i+radius ; ++k) 
	smoothMeanCurve[i]+=meanCurve[k];
      smoothMeanCurve[i]/=radius*2+1;
      smoothSdCurve[i]=0;
      for(int k=i-radius ; k<i+radius ; ++k) 
	smoothSdCurve[i]+=meanCurve[k];
      smoothSdCurve[i]/=radius*2+1;
    }
    meanCurve=smoothMeanCurve;
    sdCurve=smoothSdCurve;
  }
  for(int i=0 ; i<numBins ; ++i) {
    cout<<i*binSize<<"\t"<<smoothMeanCurve[i]<<"\t"<<smoothSdCurve[i]<<endl;
  }
}


void Application::addCounts(WigBinary &wig,WigInterval w,Vector<float> &counts)
{
  for(int i=w.begin ; i<w.end ; ++i) {
    //cout<<"====> "<<wig.read(i)<<endl;
    counts.push_back(wig.read(i));
  }
}



void Application::addCounts(Array1D<float> &A,WigInterval w,
			    Vector<float> &counts)
{
  for(int i=w.begin ; i<w.end ; ++i) counts.push_back(A[i]);
}



void Application::getRegions(Array1D<float> &A,Vector<WigInterval> &regions)
{
  const int L=A.size();
  float prev=0;
  int begin=-1;
  for(int i=0 ; i<L ; ++i) {
    //cout<<"XXXXXXX "<<i<<" -> "<<A[i]<<endl;
    if(A[i]>0 && prev==0) begin=i;
    else if(A[i]==0 && prev>0) regions.push_back(WigInterval(begin,i));
    prev=A[i];
  }
  if(A[L-1]>0) regions.push_back(WigInterval(begin,L));
}



void Application::fillArray(Array1D<float> &A,WigBinary &wig,
			    RnaJunctions &junctions) 
{
  const int L=A.size();
  for(int i=0 ; i<L ; ++i) A[i]=wig.read(i);
  const int numJunctions=junctions.getNumJunctions();
  for(int i=0 ; i<numJunctions ; ++i) {
    const RnaJunction &junction=junctions[i];
    const float depth=junction.getDepth();
    for(int j=junction.getBegin() ; j<junction.getEnd() ; ++j) {
      if(A[j]<depth) A[j]=depth;
    }
  }
}




