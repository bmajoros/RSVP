/****************************************************************
 DiskMatrix.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_RnaJunctions_H
#define INCL_RnaJunctions_H
#include <iostream>
#include "BOOM/File.H"

#define FLOATSIZE sizeof(float)

class DiskMatrix {
public:
  DiskMatrix(const String &filename,int m,int n) 
    : m(m), n(n), file(filename) {}
  DiskMatrix(FILE *f,int m,int n) : m(m), n(n), file(f) {}
  void set(int i,int j,const float &value)
    { file.seek(FLOATSIZE*(m*j+i)); file.write(value); }
  float get(int i,int j)
    { file.seek(FLOATSIZE*(m*j+i)); return file.readFloat(); }
  void close() { file.close(); }
protected:
  File file;
  const int m, n;
};

class DiskArray {
public:
  DiskArray(const String &filename) : file(filename) {}
  DiskArray(FILE *f) : file(f) {}
  void set(int i,const float &value)
    { file.seek(FLOATSIZE*i); file.write(value); }
  float get(int i)
    { file.seek(FLOATSIZE*i); return file.readFloat(); }
  void close() { file.close(); }
protected:
  File file;
}

#endif
