#! usr/bin/env python3

import numpy as np
import datetime
import sklearn.metrics as skm
import sys
import seaborn as sb
import pandas as pd
import matplotlib.pyplot as plt

def build_meta(fp,overlap):
     """
     input: path of .coord file
     output: a generator contains id, strand, and gene
     """
     with open(fp,'rb') as coord:
          for line in coord:
               if line[0] != 35:
                    item = line.split(b'\t')
                    if item[6] == b'-':
                         item[6] = 16
                    else:
                         item[6] = 0
                    if int(item[4]) - int(item[3]) >= overlap:
                         yield (item[0],item[6],
                                int(bytes.decode(item[7].split(b':')[1].split(b'S')[0])))
     coord.close()
def build_ribo(p5s,p16s,p23s):
     """
     input: path of 5s
            path of 16s
            path of 23s
     output: a numpy array contains seqid and sequences of gene
     """
     with open(p5s,'rb') as S5:
          S5.readline()
          for line in S5:
               yield (line.split(b'\t')[0],-1,5)
     S5.closed

     with open(p16s,'rb') as S16:
          S16.readline()
          for line in S16:
               yield (line.split(b'\t')[0],-1,16)
     S16.closed

     with open(p23s,'rb') as S23:
          S23.readline()
          for line in S23:
               yield (line.split(b'\t')[0],-1,23)
     S23.closed

     
     
def build_sortme(psam,p5s,p16s,p23s):
     """
     input: path of sam file
     output: a generator contains seqid,strand, and seqences of gene
     """
     
     table = {}
     with open(p5s,'r') as S5:
          for line in S5:
               if line[0] == '>':
                    table[line[1:].strip()] = 5
     S5.close()
     with open(p16s,'r') as S16:
          for line in S16:
               if line[0] == '>':
                    table[line.split(' ')[0][1:]] = 16
     S16.close()
     with open(p23s,'r') as S23:
          for line in S23:
               if line[0] == '>':
                    table[line.split(' ')[0][1:]] = 23
     S23.close()
     with open(psam,'rb') as sam:
          for line in sam:
               if line[0] != 64:
                    item = line.split(b'\t')
                    yield(item[0],int(bytes.decode(item[1])),table[bytes.decode(item[2])])
     sam.close()

def build_labels(fp):
     """
     input: path of labels file
     output: a numpy array with seqid,strand,overlap,type,species attributes
     """
     labels = np.loadtxt(fp,delimiter='\t',dtype=[('seqid','S30'),('strand','i'),('overlap','i')
                                                ,('type','S10'),('species','S100')])
     labels.sort()
     return labels

def binary_search(start,end,seq,labels):
     """
     input(first call): start index of searching
                        end index of searching
                        source
                        destination
     output: index of source in destination or none
     """
     if start>end:
          return None
     label = labels[int((start+end)/2)]['seqid']
     if seq == label:
          return int((start+end)/2)
     elif seq < label:
          return binary_search(start,int((start+end)/2)-1,seq,labels)
     else:
          return binary_search(int((start+end)/2)+1,end,seq,labels)

def find_pred_rna(target,labels,ignore_strand):
     """
     input: RNA filter
            labels
            if is ignore strand
     output: a binary numpy array of predicting rna
     """
     len_labels = len(labels)
     pred_rna = np.zeros(len_labels)
     index = 0

     for item in target:
          index = binary_search(0,len_labels-1,item[0],labels)
          if index != None:
               for i in range(-3,3):
                    if labels[index+i]['seqid'] == item[0] and (ignore_strand or labels[index+i]['strand'] == item[1]):
                         pred_rna[index+i] = item[2]
          else:
               print('can\'t find it')
     return pred_rna

def find_true_rna(labels):
     """
     input: labels
     output: a dictionary of true RNA which contains total,5S,16S,23S
     """
     len_labels = len(labels)
     true_rna = np.zeros(len_labels)

     for index in range(0,len_labels):
          if labels[index]['type'] == b'5S':
               true_rna[index] = 5
          elif labels[index]['type'] == b'16S':
               true_rna[index] = 16
          elif labels[index]['type'] == b'23S':
               true_rna[index] = 23

     return true_rna

def build_data_frame2(tool,str):
     """
     input: the list of sklearn.classification_report
     output: dataframe divided by tool, precision and recall
     """
     return pd.DataFrame({ 'tool': [str+',5S',str+',16S',str+',23S'],
                          'precision' : pd.Series([tool[10],tool[15],tool[21]],dtype='d'),
                          'recall' : pd.Series([tool[11],tool[16],tool[22]],dtype='d')},
                         index=[0,1,2])

def build_data_frame_by_hit(true,pred,labels,legend):
     """
     input: true_rna,pred_rna,labels
     output: precision of pred_rna and recall of pred_rna
     """
     count = 0
     species_tp = []
     species_tn = []
     species_fp = []
     species_fn = []
     fl = open(legend,'w')
     for i in range(0,len(labels)):
          if true[i] == pred[i]:
               if true[i] == 0:
                    species_tn.append(bytes.decode(b" ".join(labels[i]['species'].split()[0:2])))
               else:
                    fl.write(bytes.decode(labels[i][0])+'\t'+str(labels[i][1])+'\t'+str(labels[i][2])+'\t'+bytes.decode(labels[i][3])+'\t'+bytes.decode(labels[i][4])+'\n')
                    species_tp.append(bytes.decode(b" ".join(labels[i]['species'].split()[0:2])))
          else:
               if true[i] == 0:
                    count += 1
                    species_fp.append(bytes.decode(b" ".join(labels[i]['species'].split()[0:2])))
               else:
                    species_fn.append(bytes.decode(b" ".join(labels[i]['species'].split()[0:2])))
     tp = pd.DataFrame({'species' : species_tp,
                        'count' : np.zeros(len(species_tp))}).groupby('species').count()
     fp = pd.DataFrame({'species' : species_fp,
                        'count' : np.zeros(len(species_fp))}).groupby('species').count()
     fn = pd.DataFrame({'species' : species_fn,
                        'count' : np.zeros(len(species_fn))}).groupby('species').count()
     precision = tp.divide(tp.add(fp,fill_value=0.0),fill_value=0.0).reset_index(level=0).sort_values('count')
     recall = tp.divide(tp.add(fn,fill_value=0.0),fill_value=0.0).reset_index(level=0).sort_values('count')

     fl.close()
     print(count)
     return precision,recall

def display_confusion_matrix(true,pred):
     """
     input: true_rna,pred_rna
     output: a confusion matrix contains the number of hits or miss splited by types of ribosome 
     """
     cm = skm.confusion_matrix(true,pred)
     print('Confusion Matrix:\n\tnone\t\t5S\t\t16S\t\t23S')
     print('none\t'+str(cm[0][0])+'\t\t'+str(cm[0][1])+'\t\t'+str(cm[0][2])+'\t\t'+str(cm[0][3]))
     print('5S\t'+str(cm[1][0])+'\t\t'+str(cm[1][1])+'\t\t'+str(cm[1][2])+'\t\t'+str(cm[1][3]))
     print('16S\t'+str(cm[2][0])+'\t\t'+str(cm[2][1])+'\t\t'+str(cm[2][2])+'\t\t'+str(cm[2][3]))
     print('23S\t'+str(cm[3][0])+'\t\t'+str(cm[3][1])+'\t\t'+str(cm[3][2])+'\t\t'+str(cm[3][3]))

if __name__ == '__main__':
     tf = 'test_files/'
     
     labels_20 = build_labels(tf+'labels_20')
     labels_100 = build_labels(tf+'labels_100')
     
     pred_meta_20 = find_pred_rna(build_meta(tf+'mt.coord',20),labels_20,True)
     pred_meta_100 = find_pred_rna(build_meta(tf+'mt.coord',100),labels_100,True)
     pred_ribo_100 = find_pred_rna(build_ribo(tf+'rb_5s.tsv',tf+'rb_16s.tsv',tf+'rb_23s.tsv'),labels_100,True) 
     pred_sortme_100 = find_pred_rna(build_sortme(tf+'sm.sam',tf+'sm_5s.fasta',tf+'sm_16s.fasta',tf+'sm_23s.fasta'),labels_100,True)

     true_100 = find_true_rna(labels_100)
     true_20 = find_true_rna(labels_20)
     
     pmt_100,rmt_100 = build_data_frame_by_hit(true_100,pred_meta_100,labels_100,'MetaRNA')
     prb_100,rrb_100 = build_data_frame_by_hit(true_100,pred_ribo_100,labels_100,'Ribopicker')
     psm_100,rsm_100 = build_data_frame_by_hit(true_100,pred_sortme_100,labels_100,'SortmeRNA')

     recall = pd.concat([rmt_100,rrb_100,rsm_100],keys=['Meta','Ribo','Sortme']).reset_index(level=0)
     precision = pd.concat([pmt_100,prb_100,psm_100],keys=['Meta','Ribo','Sortme']).reset_index(level=0)

     

     

