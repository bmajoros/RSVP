/****************************************************************
 evaluate.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
#include "BOOM/File.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/FastaReader.H"
#include "BOOM/ProteinTrans.H"
using namespace std;
using namespace BOOM;

const char *donorConsensuses[]={"GT","GC","AT"};
const char *acceptorConsensuses[]={"AG","AC"};
int numDonorConsensuses=sizeof(donorConsensuses)/sizeof(char*);
int numAcceptorConsensuses=sizeof(acceptorConsensuses)/sizeof(char*);


/****************************************************************
                       class LightExon
 ****************************************************************/
struct LightExon {
  int begin, end;
  float score;
  ExonType exonType;
  LightExon(int b,int e,float s,ExonType t) : begin(b), end(e), score(s),
					      exonType(t) {}
  bool operator==(const LightExon &b) const
  { return begin==b.begin && end==b.end; }
  bool operator>(const LightExon &b) const
  { return begin>b.begin || begin==b.begin && end>b.end; }
  bool operator<(const LightExon &b) const
  { return begin<b.begin || begin==b.begin && end<b.end; }
};

/****************************************************************
                       class ExonComparer
 ****************************************************************/
class ExonComparer : public Comparator<LightExon> {
public:
  bool equal(LightExon &a,LightExon &b)   { return a==b; }
  bool greater(LightExon &a,LightExon &b) { return a>b; }
  bool less(LightExon &a,LightExon &b)    { return a<b; }
};

/****************************************************************
                      class LightTranscript
 ****************************************************************/
struct LightTranscript {
  Vector<LightExon> exons;
  float score;
  void sortExons();
  bool operator==(const LightTranscript &);
  bool internallyEqual(const LightTranscript &); // only internal exons
};

/****************************************************************
                   class LightTransComparator
 ****************************************************************/
struct LightTransComparator : Comparator<LightTranscript> {
  bool equal(LightTranscript &a,LightTranscript &b) 
    { return a.score==b.score; }
  bool greater(LightTranscript &a,LightTranscript &b)
    { return a.score>b.score; }
  bool less(LightTranscript &a,LightTranscript &b)
    { return a.score<b.score; }
};

/****************************************************************
                      class PredictionIndex
 ****************************************************************/
struct PredictionIndex {
  int chunkID;
  int predictionID;
  PredictionIndex(int c,int p) : chunkID(c), predictionID(p) {}
  PredictionIndex() {}
};

/****************************************************************
                        class Substrate
 ****************************************************************/
struct Substrate {
  Vector<LightTranscript> annotations;
  Vector<LightTranscript> predictions;
};

/****************************************************************
                       class Application
 ****************************************************************/
class Application {
  bool internalOnly;
  bool sameNumbers; // use same #'s of transcripts as in another GFF file
  Vector<Substrate> chunks;
  int numIdenticalExons;// all exon types
  int numUniqPredExons; // all exon types
  int numUniqAnnoExons; // all exon types
  int numIdentExonsIntern;    // internal exons only
  int numUniqPredInternExons; // internal exons only
  int numUniqAnnoInternExons; // internal exons only
  int numIdentInternTrans;    // internal exons only
  int numIdenticalTranscripts;// full gene
  int numAnnoTranscripts; // full transcripts
  int numAnnoInternTrans; // internal exons only
  int numPredTranscripts; // full transcripts
  int numPredInternTrans; // internal exons only
  int maxTranscripts;      // full transcripts
  int maxInternTranscripts;//internal exons only
  bool allChunksInMemory;
  bool ensureCanonical;
  Map<String,int> chunkCounts; // chunkID -> #transcripts
  int chunkLimit; // max #transcripts for any chunk
  void applySameNumbers(Vector<LightTranscript> &,const String &chunkID);
  void load(const String &filename,Vector<LightTranscript> &into,
	    const String &sequence);
  void load(const String &filename,Vector<LightTranscript> &into);
  void updateStats(Vector<LightTranscript> &predictions,
		   Vector<LightTranscript> &annotations);
  void uniqExons(Vector<LightTranscript> &,Vector<LightExon> &into,
		 Set< pair<int,int> > &seen);
  void uniqExons(LightTranscript &transcript,Vector<LightExon> &into,
		 Set< pair<int,int> > &seen);
  void sortExons(Vector<LightExon> &exons);
  void updateExonStats(Vector<LightExon> &predictions,
		       Vector<LightExon> &annotations);
  int countExons(Vector<LightExon> &,ExonType);
  void getInternalTranscripts(Vector<LightTranscript> &from,
			      Vector<LightTranscript> &to);
  void getInternalExons(Vector<LightExon> &from,Vector<LightExon> &to);
  String getKey(Vector<LightExon> &);
  bool canonicalSplicing(GffExon &,const String &seq);
  void limitTranscripts(int);
  void limitInternTranscripts(int);
  void getChunkCounts(const String &gffFilename);
public:
  Application();
  int main(int argc,char *argv[]);
};

/****************************************************************
                       class PredIndexCmp
 ****************************************************************/
struct PredIndexCmp : Comparator<PredictionIndex> {
  PredIndexCmp(Vector<Substrate> &chunks) : chunks(chunks) {}
  bool equal(PredictionIndex &a,PredictionIndex &b) 
  { return chunks[a.chunkID].predictions[a.predictionID].score==
      chunks[b.chunkID].predictions[b.predictionID].score; }
  bool greater(PredictionIndex &a,PredictionIndex &b)
  { return chunks[a.chunkID].predictions[a.predictionID].score>
      chunks[b.chunkID].predictions[b.predictionID].score; }
  bool less(PredictionIndex &a,PredictionIndex &b)
  { return chunks[a.chunkID].predictions[a.predictionID].score<
      chunks[b.chunkID].predictions[b.predictionID].score; }
private:
  Vector<Substrate> &chunks;
};



/****************************************************************
                             main()
 ****************************************************************/
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



/****************************************************************
                     Application methods
 ****************************************************************/


Application::Application()
  : numIdenticalExons(0), numUniqPredExons(0), numUniqAnnoExons(0),
    numIdentExonsIntern(0), numUniqPredInternExons(0), 
    numUniqAnnoInternExons(0), numIdentInternTrans(0),
    numIdenticalTranscripts(0), numAnnoTranscripts(0),
    numAnnoInternTrans(0), numPredTranscripts(0),
    numPredInternTrans(0)
{
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"it:T:S:C:c");
  if(cmd.numArgs()!=3)
    throw String("\n\
evaluate <predictions-dir> <chunks-dir> <ext>\n\
  -i = internal exons only\n\
  -T n = keep top n transcripts\n\
  -t n = keep top n internal transcripts\n\
  -S *.gff = keep same # of transcripts per chunk as in other gff file\n\
  -C n = keep no more than n transcripts for each chunk\n\
  -c = don't enforce splice consensus\n\
");
  String predictionDir=cmd.arg(0);
  String chunksDir=cmd.arg(1);
  String ext=cmd.arg(2);
  internalOnly=cmd.option('i');
  maxTranscripts=cmd.option('T') ? cmd.optParm('T').asInt() : INT_MAX;
  maxInternTranscripts=cmd.option('t') ? cmd.optParm('t').asInt() : INT_MAX;
  allChunksInMemory=cmd.option('T') || cmd.option('t');
  if(cmd.option('T') && cmd.option('t'))
    throw "Can't use both -T and -t at the same time";
  sameNumbers=cmd.option('S');
  if(sameNumbers) getChunkCounts(cmd.optParm('S'));
  chunkLimit=cmd.option('C') ? cmd.optParm('C').asInt() : INT_MAX;
  ensureCanonical=!cmd.option('c');

  Vector<String> predictionFiles;
  File::getFileList(predictionDir,predictionFiles);
  int numFiles=predictionFiles.size();
  for(int i=0 ; i<numFiles ; ++i) {
    String predFile=predictionFiles[i];
    Vector<String> &fields=*predFile.getFields("/.");
    int numFields=fields.size();
    Substrate chunk;
    if(numFields>1 && fields[numFields-1]==ext) {
      const String &substrate=fields[numFields-2];
      String fastaFile=chunksDir+'/'+substrate+".fasta";
      FastaReader reader(fastaFile);
      String defline, seq;
      reader.nextSequence(defline,seq);
      load(predictionDir+'/'+predFile,chunk.predictions,seq);
      //cout<<chunk.predictions.size()<<" transcripts"<<endl;
      String annoPath=chunksDir+'/'+substrate+".gff";
      load(annoPath,chunk.annotations,seq);
      if(sameNumbers) applySameNumbers(chunk.predictions,substrate);
      if(chunk.predictions.size()>chunkLimit)
	chunk.predictions.resize(chunkLimit);
      if(allChunksInMemory) chunks.push_back(chunk);
    }
    delete &fields;
    if(!allChunksInMemory) updateStats(chunk.predictions,chunk.annotations);
  }

  if(allChunksInMemory) {
    int numChunks=chunks.size();
    if(cmd.option('T')) limitTranscripts(maxTranscripts);
    else if(cmd.option('t')) limitInternTranscripts(maxInternTranscripts);
    for(int i=0 ; i<numChunks ; ++i) 
      updateStats(chunks[i].predictions,chunks[i].annotations);
  }

  cout<<"INTERNAL EXONS: Sn="
      <<int(1000*numIdentExonsIntern/float(numUniqAnnoInternExons)+5/9.0)/10.0
      <<"% #preds="<<numUniqPredInternExons<<" (internal exons only)"<<endl;
  cout<<"ALL EXONS: Sn="<<int(1000*numIdenticalExons/float(numUniqAnnoExons)
			      +5/9.0)/10.0
      <<"% #preds="<<numUniqPredExons<<endl;
  cout<<"IDENTICAL FULL TRANSCRIPTS: Sn="
      <<int(1000*numIdenticalTranscripts/float(numAnnoTranscripts)+5/9.0)/10.0
      <<"%  "<<numIdenticalTranscripts<<"/"<<numAnnoTranscripts<<endl;
  cout<<"IDENTICAL INTERNAL TRANSCRIPTS: Sn="
      <<int(1000*numIdentInternTrans/float(numAnnoInternTrans)+5/9.0)/10.0
      <<"%  "<<numIdentInternTrans<<"/"<<numAnnoInternTrans<<endl;
  cout<<"predicted transcripts: "<<numPredTranscripts
      <<" ("<<numPredInternTrans<<" internal)"<<endl;
  return 0;
}



void Application::load(const String &filename,Vector<LightTranscript> &into,
		       const String &seq)
{
  GffReader reader(filename);
  Vector<GffTranscript*> &transcripts=*reader.loadTranscripts();
  int numTranscripts=transcripts.size();
  Set<String> seen;
  for(int i=0 ; i<numTranscripts ; ++i) {
    GffTranscript *transcript=transcripts[i];
    LightTranscript trans;
    int numExons=transcript->getNumExons();
    for(int j=0 ; j<numExons ; ++j) {
      BOOM::GffExon &exon=transcript->getIthExon(j);
      if(canonicalSplicing(exon,seq)) {
	LightExon le(exon.getBegin(),exon.getEnd(),exon.getScore(),
		     exon.getExonType());
	trans.exons.push_back(le);
      }	
    }
    trans.sortExons();
    String key=getKey(trans.exons);
    if(seen.isMember(key) || trans.exons.size()==0) 
      {delete transcript; continue;}
    seen.insert(key);
    trans.score=transcript->getScore();
    into.push_back(trans);
    delete transcript;
  }
  delete &transcripts;
}



void Application::load(const String &filename,Vector<LightTranscript> &into)
{
  GffReader reader(filename);
  Vector<GffTranscript*> &transcripts=*reader.loadTranscripts();
  int numTranscripts=transcripts.size();
  Set<String> seen;
  for(int i=0 ; i<numTranscripts ; ++i) {
    GffTranscript *transcript=transcripts[i];
    LightTranscript trans;
    int numExons=transcript->getNumExons();
    for(int j=0 ; j<numExons ; ++j) {
      BOOM::GffExon &exon=transcript->getIthExon(j);
      LightExon le(exon.getBegin(),exon.getEnd(),exon.getScore(),
		   exon.getExonType());
      trans.exons.push_back(le);
    }
    trans.sortExons();
    String key=getKey(trans.exons);
    if(seen.isMember(key) || trans.exons.size()==0) 
      {delete transcript; continue;}
    seen.insert(key);
    trans.score=transcript->getScore();
    into.push_back(trans);
    delete transcript;
  }
  delete &transcripts;
}



void Application::uniqExons(LightTranscript &transcript,
			    Vector<LightExon> &into,
			    Set< pair<int,int> > &seen)
{
  for(Vector<LightExon>::iterator cur=transcript.exons.begin(),
	end=transcript.exons.end() ; cur!=end ; ++cur) {
    const LightExon &exon=*cur;
    pair<int,int> p(exon.begin,exon.end);
    if(seen.isMember(p)) continue;
    seen.insert(p);
    into.push_back(exon);
  }
}



void Application::uniqExons(Vector<LightTranscript> &transcripts,
			    Vector<LightExon> &into,
			    Set< pair<int,int> > &seen)
{
  for(Vector<LightTranscript>::iterator cur=transcripts.begin(), end=
	transcripts.end() ; cur!=end ; ++cur) 
    uniqExons(*cur,into,seen);
}
  
  

void Application::sortExons(Vector<LightExon> &exons)
{
  ExonComparer cmp;
  VectorSorter<LightExon> sorter(exons,cmp);
  sorter.sortAscendInPlace();
}



void Application::updateExonStats(Vector<LightExon> &predictions,
				  Vector<LightExon> &annotations)
{
  int numPreds=predictions.size(), numAnno=annotations.size();
  for(int p=0, a=0 ; p<numPreds && a<numAnno ; ) {
    const LightExon &pExon=predictions[p], &aExon=annotations[a];
    if(pExon==aExon) {
      ++numIdenticalExons;
      if(aExon.exonType==ET_INTERNAL_EXON) ++numIdentExonsIntern;
      ++p; ++a;
      continue;
    }
    if(pExon<aExon) ++p;
    else ++a;
  }
  numUniqPredExons+=numPreds;
  numUniqAnnoExons+=numAnno;
}



void Application::updateStats(Vector<LightTranscript> &predictions,
			      Vector<LightTranscript> &annotations) 
{
  // Exons:
  Vector<LightExon> uniqPredExons, uniqAnnoExons;
  Set< pair<int,int> > predSeen, annoSeen;
  uniqExons(predictions,uniqPredExons,predSeen); sortExons(uniqPredExons);
  uniqExons(annotations,uniqAnnoExons,annoSeen); sortExons(uniqAnnoExons);
  updateExonStats(uniqPredExons,uniqAnnoExons);
  numUniqPredInternExons+=countExons(uniqPredExons,ET_INTERNAL_EXON);
  numUniqAnnoInternExons+=countExons(uniqAnnoExons,ET_INTERNAL_EXON);

  // Whole genes:
  int numPredTrans=predictions.size(), numAnnoTrans=annotations.size();
  for(int i=0 ; i<numPredTrans ; ++i)
    for(int j=0 ; j<numAnnoTrans ; ++j)
      if(predictions[i]==annotations[j]) ++numIdenticalTranscripts;
  Vector<LightTranscript> internalPred, internalAnno;
  getInternalTranscripts(predictions,internalPred);
  getInternalTranscripts(annotations,internalAnno);
  int numInternPredTrans=internalPred.size();
  int numInternAnnoTrans=internalAnno.size();
  for(int i=0 ; i<numInternPredTrans ; ++i)
    for(int j=0 ; j<numInternAnnoTrans ; ++j)
      if(internalPred[i]==internalAnno[j]) ++numIdentInternTrans;
  numAnnoTranscripts+=numAnnoTrans;
  numAnnoInternTrans+=numInternAnnoTrans;
  numPredTranscripts+=numPredTrans;
  numPredInternTrans+=numInternPredTrans;
}



int Application::countExons(Vector<LightExon> &exons,ExonType type)
{
  int count=0;
  for(Vector<LightExon>::iterator cur=exons.begin(), end=exons.end() ; 
      cur!=end ; ++cur)
    if((*cur).exonType==type) ++count;
  return count;
}



String Application::getKey(Vector<LightExon> &exons) 
{
  String key;
  for(Vector<LightExon>::iterator cur=exons.begin(), end=exons.end() ; 
      cur!=end ; ++cur) {
    const LightExon &exon=*cur;
    key+=String(exon.begin)+","+exon.end+" ";
  }
  return key;
}



void Application::getInternalTranscripts(Vector<LightTranscript> &from,
					 Vector<LightTranscript> &to)
{
  Set<String> seen;
  for(Vector<LightTranscript>::iterator cur=from.begin(), end=from.end() ; 
      cur!=end ; ++cur) {
    LightTranscript t;
    getInternalExons((*cur).exons,t.exons);
    if(t.exons.isEmpty()) continue;
    String key=getKey(t.exons);
    if(seen.isMember(key)) continue;
    seen.insert(key);
    to.push_back(t);
  }
}	



void Application::getInternalExons(Vector<LightExon> &from,
				   Vector<LightExon> &to)
{
  for(Vector<LightExon>::iterator cur=from.begin(), end=from.end() ; 
      cur!=end ; ++cur) {
    const LightExon &exon=*cur;
    if(exon.exonType==ET_INTERNAL_EXON) to.push_back(exon);
  }	
}



bool Application::canonicalSplicing(GffExon &exon,const String &seq)
{
  if(!ensureCanonical) return true;
  //return true; // ###

  const int L=seq.length();
  Vector<String> donors, acceptors;
  int b=exon.getBegin(), e=exon.getEnd();
  if(exon.getStrand()=='+') 
    switch(exon.getExonType())
      {
      case ET_INITIAL_EXON:
	if(e>=0 && e<L-2) donors.push_back(seq.substring(e,2));
	break;
      case ET_INTERNAL_EXON:
	if(b-2>=0 && b-2<L-2) acceptors.push_back(seq.substring(b-2,2));
	if(e>=0 && e<L-2) donors.push_back(seq.substring(e,2));
	break;
      case ET_FINAL_EXON:
	if(b-2>=0 && b-2<L-2) acceptors.push_back(seq.substring(b-2,2));
	break;
      case ET_SINGLE_EXON:
	break;
      }
  else switch(exon.getExonType())
    {
    case ET_INITIAL_EXON:
      if(b-2>=0 && b-2<L-2) donors.push_back(ProteinTrans::
		       reverseComplement(seq.substring(b-2,2)));
      break;
    case ET_INTERNAL_EXON:
      if(e>=0 && e<L-2) acceptors.push_back(ProteinTrans::
			  reverseComplement(seq.substring(e,2)));
      if(b-2>=0 && b-2<L-2) donors.push_back(ProteinTrans::
		       reverseComplement(seq.substring(b-2,2)));
      break;
    case ET_FINAL_EXON:
      if(e>=0 && e<L-2) acceptors.push_back(ProteinTrans::
			  reverseComplement(seq.substring(e,2)));
      break;
    case ET_SINGLE_EXON:
      break;
    }
  for(Vector<String>::iterator cur=donors.begin(), end=donors.end() ; 
      cur!=end ; ++cur) {
    const String &donor=*cur;
    bool found=false;
    for(int i=0 ; i<numDonorConsensuses ; ++i)
      if(donor==donorConsensuses[i]) found=true;
    if(!found) return false;
  }
  for(Vector<String>::iterator cur=acceptors.begin(), end=acceptors.end() ; 
      cur!=end ; ++cur) {
    const String &acceptor=*cur;
    bool found=false;
    for(int i=0 ; i<numAcceptorConsensuses ; ++i)
      if(acceptor==acceptorConsensuses[i]) found=true;
    if(!found) return false;
  }
  return true;
}



void Application::limitTranscripts(int N)
{
  Vector<PredictionIndex> V;
  const int numChunks=chunks.size();
  for(int i=0 ; i<numChunks ; ++i) {
    Vector<LightTranscript> &predictions=chunks[i].predictions;
    int numPred=predictions.size();
    for(int j=0 ; j<numPred ; ++j)
      V.push_back(PredictionIndex(i,j));
  }
  if(V.size()<=N) return;
  PredIndexCmp cmp(chunks);
  VectorSorter<PredictionIndex> sorter(V,cmp);
  sorter.sortDescendInPlace();
  V.resize(N);
  Array1D< Vector<int> > keep(numChunks);
  for(int i=0 ; i<N ; ++i) {
    PredictionIndex index=V[i];
    keep[index.chunkID].push_back(index.predictionID);
  }
  V.purge();
  for(int i=0 ; i<numChunks ; ++i) {
    Substrate &chunk=chunks[i];
    Vector<LightTranscript> &predictions=chunk.predictions;
    Vector<LightTranscript> kept;
    Vector<int> toKeep=keep[i];
    for(Vector<int>::iterator cur=toKeep.begin(), end=toKeep.end() ; 
	cur!=end ; ++cur)
      kept.push_back(predictions[*cur]);
    predictions=kept;
  }
}



void Application::limitInternTranscripts(int N)
{
}



void Application::getChunkCounts(const String &gffFilename)
{
  Vector<LightTranscript> transcripts;
  //load(gffFilename,transcripts);
  GffReader reader(gffFilename);
  Map<String,TranscriptList*> &byContig=*reader.loadByContig();
  Set<String> &keys=*byContig.getKeys();
  for(Set<String>::iterator cur=keys.begin(),end=keys.end();cur!=end;++cur) {
    const String &key=*cur;
    chunkCounts[key]=byContig[key]->size();
  }
  delete &keys;
  delete &byContig;

}



void Application::applySameNumbers(Vector<LightTranscript> &array,
				   const String &chunkID)
{
  int chunkCount=chunkCounts[chunkID];
  if(chunkCount>0 && array.size()<=chunkCount) return;

  LightTransComparator cmp;
  VectorSorter<LightTranscript> sorter(array,cmp);
  sorter.sortDescendInPlace();

  array.resize(chunkCount);
}



/****************************************************************
                     LightTranscript methods
 ****************************************************************/


void LightTranscript::sortExons()
{
  ExonComparer cmp;
  VectorSorter<LightExon> sorter(exons,cmp);
  sorter.sortAscendInPlace();
}



bool LightTranscript::operator==(const LightTranscript &other)
{
  int n=exons.size(), n2=other.exons.size();
  if(n!=n2) return false;
  for(int i=0 ; i<n ; ++i) if(!(exons[i]==other.exons[i])) return false;
  return true;
}


bool LightTranscript::internallyEqual(const LightTranscript &other)
{
  int n=exons.size(), n2=other.exons.size();
  if(n!=n2) return false;
  for(int i=0 ; i<n ; ++i) {
    LightExon &exon=exons[i];
    if(exon.exonType!=ET_INTERNAL_EXON) continue;
    if(!(exons[i]==other.exons[i])) return false;
  }	
  return true;
}



