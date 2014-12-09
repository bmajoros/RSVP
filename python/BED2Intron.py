#!/home/sysbio/sl228/Softwares/Python-2.7.1/bin/python2.7
'''
Created on Mar 17, 2011

@author: sl228

given a BED file, which usually has one header line,
and the following structured entries: 

Chr2    4587    6211    JUNC00000002    1       +       4587    6211    255,0,0 2       23,27   0,1597

write a *.BED.intron file of tab delimited file

Chr2    leftofintron    rightofintron
 ||        ||           ||
Chr2    4587+23    6211-27 #need to test it's 4587+23+/-1 and 6211-27+/-1

@tested on May 09, 2011

#for the left of the location:
4587 is 0 based genomic location, 
4587+1 is translated to 1 based genomic location is 4588
23 is the length of the thick area.
4587 + 1 + 23 - 1 is the end of the thick area
4587 + 1 + 23 -1 + 1 is the start of the thin area (intron) in 1 based genomic location
simplify to 4587 + 1 + 23 is the start of the thin area(intron) in 1 based genomic location
#above explain the line of code:
left = int(ct[1])+int(hangingend[0]) + 1

#for the right of the location:
6211 is position that is just outside the thick area in 0 based genomic location
6211 + 1 is the same position in 1 base genomic location
6211 + 1 -1 is the last base location of the thick area on 1 based genomic location
6211 + 1 -1 -27 + 1 is the first base location of the thick area in 1 based genomic location
6211 + 1 -1 -27 + 1 - 1 is the last base of thin area in 1 based genomic location
simplify to 6211 - 27 is the last base of thin area in 1 based genomic location

#above explain the line of code:
right = int(ct[2])-int(hangingend[1])

@ Mar 24, 2011
To keep the number of reads that cover certain junction, which is column 5(4 in Python).


'''
import sys
def BED2intron(infilename):
    infile = open(infilename,'r')
    outfile = open(infilename+'.intron','w')
    for line in infile:
        if line.find('name=')>0:
            continue
        ct = line.split('\t')
        try:
            hangingend = ct[10].split(',')
            chrname = ct[0]
            left = int(ct[1])+int(hangingend[0]) + 1
            right = int(ct[2])-int(hangingend[1])
            count = int(ct[4])
            strand = ct[5]
            outfile.write(chrname+'\t'+`left`+'\t'+`right`+'\t'+`count`+'\t'+strand+'\n')
        except IndexError:
            print ct
            break
    
    outfile.close()
    infile.close()


if __name__ == '__main__':
    BED2intron(sys.argv[1])
    pass
