#!usr/bin/env python3
"""
The script will compare given sam file with gff3 file. 
Then it will print out the labels for each read.
Each label contains readID, rRNA label, the number of overlapped base pairs, and speciesID
"""
import sys
import readFile

def cont_overlap(read_start,rRNA_start,rRNA_end):
    """
    input: start position of read
           start position of rRNA
           end position of rRNA
    output: the number of overlapped base pairs
    """
    if read_start+250 >= rRNA_start and read_start <= rRNA_end:
        return min(rRNA_end-read_start,
                   read_start+250-rRNA_start,
                   250,
                   rRNA_end-rRNA_start)
    else:
        return 0
    
def create_labels(gff3_path,sam_path,ignore_strand,bp):
    """
    input: path of gff3 file
           path of sam path
    output: A generator contains labels for each read.
            Each label has  readID
                           ,rRNA label
                           ,the number of overlapped base pairs
                           ,and speciesID
    """
    gff3 = readFile.read_gff3(gff3_path)
    sam = readFile.read_sam(sam_path)
    li =[]
    for rRNA in gff3:
        li.append(rRNA)

    length = len(li)

    for read in sam:
        index = 0
        isntRibo = True
        while True:                
            overlap = cont_overlap(int(read[3]),int(li[index][3]),int(li[index][4]))
            if read[2] == li[index][0] and (ignore_strand or read[1] == li[index][6]) and overlap >= int(bp):
                isntRibo = False
                yield (
                    read[0]+'\t'+read[1]
                    +'\t'+str(overlap)
                    +'\t'+ li[index][8].split(';')[0].split('=')[1].split('_')[0]
                    +'\t'+read[-1]
                )
            if index == length -1:
                if isntRibo:
                     yield (read[0]+'\t'+read[1]+'\t0\tnone+\t'+read[-1])
                break
            index += 1
if __name__ == '__main__':
    for label in create_labels(sys.argv[1],sys.argv[2],True,sys.argv[3]):
        print(label)                    
