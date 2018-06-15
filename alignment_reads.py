#!usr/bin/env ipython3

import numpy as np
import readFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.pipeline import make_pipeline
from sklearn import preprocessing
from sklearn import svm
import sklearn as sk
import sklearn.metrics as skm
import sklearn.model_selection as skmod
import os

class state(object):
    def __init__(self):
        self.next = [None,None,None,None,None,None,None,None]
        self.count = 0
        self.kmer = ''
class dataset(object):
    def __init__(self,fp_labels,fp_fna):
        self.data = load_traning_data(fp_labels,fp_fna)
        self.species = self.data['species']
        self.seq = self.data['seq']
    def build_metrics(self,num):
        return construct_kmer(num,self.species,self.seq)
    
def load_traning_data(fp_labels,fp_fna):
    fna = {}
    for items in readFile.read_fna(fp_fna):
        fna[items[0]] = items[1]

    labels = []
    for items in readFile.read_labels(fp_labels):
        labels.append((' '.join(items[4].split()[0:2]),items[3],int(items[2]),items[0],fna[items[0]]))
    data = np.array(labels,[('species','<U30'),('type','<U30'),('overlap','i'),('reads','<U30'),('seq','<U400')])
    
    return data

def print_confusion_matrix(confusion_matrix, class_names, figsize = (10,7), fontsize=10):
    """Prints a confusion matrix, as returned by sklearn.metrics.confusion_matrix, as a heatmap.
    
    Arguments
    ---------
    confusion_matrix: numpy.ndarray
        The numpy.ndarray object returned from a call to sklearn.metrics.confusion_matrix. 
        Similarly constructed ndarrays can also be used.
    class_names: list
         An ordered list of class names, in the order they index the given confusion matrix.
    figsize: tuple
        A 2-long tuple, the first value determining the horizontal size of the ouputted figure,
        the second determining the vertical size. Defaults to (10,7).
    fontsize: int
        Font size for axes labels. Defaults to 14.
        
    Returns
    -------
    matplotlib.figure.Figure
        The resulting confusion matrix figure
    """
    df_cm = pd.DataFrame(
        confusion_matrix, index=class_names, columns=class_names, 
    )
    fig = plt.figure(figsize=figsize)
    try:
        heatmap = sns.heatmap(df_cm, annot=True, fmt="d")
    except ValueError:
        raise ValueError("Confusion matrix values must be integers.")
    heatmap.yaxis.set_ticklabels(heatmap.yaxis.get_ticklabels(), rotation=0, ha='right', fontsize=fontsize)
    heatmap.xaxis.set_ticklabels(heatmap.xaxis.get_ticklabels(), rotation=45, ha='right', fontsize=fontsize)
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    return fig

def print_prob_distribution(prob_dtb,label,title,figsize=(10,7), fontsize = 10):
    """
    input: probability distribution of seq
           label of species
           title of figure
    
    output: a matplot barchart in /test_figures
    """
    df = pd.DataFrame(prob_dtb,index=label,columns = ['value'])
    df = df.reset_index(level=0)
    fig = plt.figure(figsize=figsize)
    barplot = sns.barplot(x='value',y='index',data=df)
    plt.title(title)
    plt.subplots_adjust(left=0.4)
    plt.savefig('./test_figures/'+title+'.png')
    plt.close(fig)
    return fig
    
def build_confusion_matrix(true,pred,label,isProba,scale,save = False):
    """
    input: a True matrix,
           a pred matrix or probability matrix,
           a list of label whose length is equal to the length of unique items in matrix,
           a boolean isProba show whether it is a probability matrix,
           the allowed distance between the first prob and second,
           
           
    """
    if isProba:
        pred_indices = find_pred_matrix(pred,scale)
    
    new_pred = np.zeros(len(true))
    
    for i in range(0,len(true)):
        if true[i] in pred_indices[i]:
            if true[i] != pred_indices[i][0]:
                if save:
                    title = 'Ambiguous_index:'+str(i)+'_species:'+str(label[true[i]])
                    print_prob_distribution(pred[i],label,title)
                new_pred[i] = 0
            else:
                new_pred[i] = true[i]+1
        else:
            if save:
                title = 'Misspredict_index:'+str(i)+'_species:'+str(label[true[i]])
                print_prob_distribution(pred[i],label,title)
            new_pred[i] = pred_indices[i][0]+1
    
    return new_pred,print_confusion_matrix(skm.confusion_matrix([i+1 for i in true],new_pred),['Ambiguous']+label)

def find_pred_matrix(pred_prob,scale):
    """
    input: prediction of probabilities for each genomes
           scale of allowed error
    
    output: prediction of each sequence by indice 
    """
    pred_indices = []
    for i in range(0,len(pred_prob)):
        max_index = pred_prob[i].argmax()
        curr=[max_index]
        for j in range(0,len(pred_prob[i])):
            if pred_prob[i][j] >curr[0]*scale and j != max_index:
                curr.append(j)
        pred_indices.append(curr)
    return pred_indices
    
def construct_kmer(num,species,sequences):
    """
    input: num of kmer
           array of species
           array of sequences
    output: array of kmer
            array of label for species
    """
    index = 0
    dict_species = {}
    for i in np.unique(species):
        dict_species[i] = index
        index += 1
    list_species = [dict_species[i] for i in species]
    
    dict_kmer = {'A':0,'C':1,'G':2,'T':3,
                 'N':4,'W':5,'Y':6,'M':7}
    array = []
    root = state()
    root.seq = '1'
    states = []
    for seq in sequences:
        #reset list
        #for item in list:
        for item in states:
            item.count = 0
        #contruct initial state
        if len(seq) < num:
            print('error! out of index')
            break
        for index in range(0,len(seq)-num):
            curr = root
            isNew = False
            kmer = seq[index:index+num]
            for char in kmer:
                if not curr.next[dict_kmer[char]]:
                    curr.next[dict_kmer[char]] = state()
                    isNew = True
                else:
                    isNew = False
                curr = curr.next[dict_kmer[char]]
            curr.count = curr.count + 1
            if isNew:
                states.append(curr)
                curr.kmer = kmer
        array.append([i.count for i in states])
    list_kmers = np.zeros((len(sequences),len(states)))
    for i in range(0,len(array)):
        for j in range(0,len(array[i])):
            list_kmers[i,j] = array[i][j]
    return pd.DataFrame(data=list_kmers), list_species

def build_index(fp_fastam,prefix):
    os.system('bwa index '+fp_fasta+'-p '+prefix)

def bwa_sequences(prefix,reads):
    os.system('echo "" > curr_reads')
    for i in reads:
        output = os.system('echo "%s\\n%s" >> curr_reads'%(i[0],i[1]))
        if output != 0:
            print('errors ocurr in curr reads')
    os.system('bwa mem bwa_index/'+prefix+' curr_reads > alignments/bwa.sam')
    
def build_db(fp_genomes):
    output = os.system('bash ./build_index.sh '+fp_genomes)
    if output != 0:
        print('errors ocurr in file \'build_index.sh\'')
        return
    
def blast_sequences(preds,reads):
    """
    input: predictions of sequence
           reads
    
    output: files have alignments       
    """
    
    for index in range(0,len(preds)):
        output = os.system('cat '+' '.join(['db/'+str(i) for i in preds[index]])+' > curr_db')
        if output != 0:
            print('errors ocurr in curr database')
            return
        output = os.system('echo "%s\\n%s" > curr_reads'%(reads[index][0],reads[index][1]))
        if output != 0:
            print('errors ocurr in curr reads')
            return
        output = os.system('blastn -query curr_reads -outfmt 6 -subject curr_db -out "alignments/%s"'%(reads[index][0]))
        print('blastn -query curr_reads -outfmt 6 -subject curr_db -out "%s"'%(reads[index][0]),output)

if __name__ == '__main__': 
    dataset_meta_100 = dataset('test_files/traning_labels_MetaRNA'
                               ,'test_files/simulated.fna') 
    kmers, species = dataset_meta_100.build_metrics(4)
    X_train, X_test, y_train, y_test = skmod.train_test_split(kmers,species,test_size = 0.3,random_state =10086)
    label = []
    for i in np.unique(dataset_meta_100.data['species']):
        label.append(i)

    clf = make_pipeline(preprocessing.StandardScaler(),sk.neighbors.KNeighborsClassifier(n_neighbors=10,weights='distance',leaf_size=120,p=3))
    cv = skmod.ShuffleSplit(n_splits=10, test_size=0.1,random_state=0)
    
    fit = clf.fit(X_train,y_train)
    
    predict_proba = fit.predict_proba(X_test)
    #pred, cm = build_confusion_matrix(y_test,predict_proba,label,True,0.1,save = True)
    #predict = fit.predict(X_test)
    #print_confusion_matrix(skm.confusion_matrix([i for i in y_test],predidct),label)
    
    reads = [('>'+dataset_meta_100.data[i]['reads'],dataset_meta_100.data[i]['seq']) for i in X_test.index]
    #build_db('../genomes')
    #blast_sequences(predict_indices,reads)
    bwa_sequences('synthetic_metagenome.fna',reads)

    


