import numpy as np
import readFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.pipeline import make_pipeline
from sklearn import preprocessing
from sklearn import svm
import sklearn.metrics as skm
import sklearn.model_selection as skmod

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
    df = pd.DataFrame(prob_dtb,index=label,columns = ['value'])
    df = df.reset_index(level=0)
    fig = plt.figure(figsize=figsize)
    barplot = sns.barplot(x='value',y='index',data=df)
    plt.title(title)
    plt.subplots_adjust(left=0.4)
    plt.savefig('./prob_distribution_figures/'+title)
    plt.clf()
    return fig
    
def build_confusion_matrix(true,pred,label,isProba,distance):
    if isProba:
        new_pred = []
        for i in range(0,len(pred)):
            curr = [pred[i].argmax()]
            max = pred[i][0]
            max_index = 0
            for j in range(0,len(pred[i])):
                 if pred[i][j] > max and j != curr[0]:
                      max_index = j
                      max = pred[i][j]
            if pred[i][curr[0]] - max < distance:
                 curr.append(max_index)
            new_pred.append(curr)
        pred_index = new_pred
    
    new_pred = np.zeros(len(true))
    
    for i in range(0,len(true)):
        print('true:'+str(true[i])+' pred:'+str(pred_index[i]))
        if true[i] != pred_index[i][0]:
            title = 'index:'+str(i)+' species:'+''+str(label[true[i]])
            print_prob_distribution(pred[i],label,title)
        if true[i] in pred_index[i]:
            new_pred[i] = true[i]
        else:
            new_pred[i] = pred_index[i][0]
    
    return new_pred,print_confusion_matrix(skm.confusion_matrix(true,new_pred),label)

def construct_kmer(num,species,sequences):
    index = 0
    dict_species = {}
    for i in np.unique(species):
        dict_species[i] = index
        index += 1
    list_species = [dict_species[i] for i in species]
    
    dict_kmer = {'A':0,'C':1,'G':2,'T':3,
                 'N':4,'W':4,'Y':4,'M':4}
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
    return list_kmers, list_species


def training(X_train,y_train):
    clf = make_pipeline(preprocessing.StandardScaler(),svm.SVC(C=1,probability=True))
    cv = skmod.ShuffleSplit(n_splits=10, test_size=0.1,random_state=0)

    fit = clf.fit(X_train,y_train)
    
    return fit

if __name__ == '__main__':
    dataset_meta_100 = dataset('test_files/traning_labels_MetaRNA'
                               ,'test_files/simulated.fna') 
    kmers, species = dataset_meta_100.build_metrics(8)
    X_train, X_test, y_train, y_test = skmod.train_test_split(kmers,species,test_size = 0.3,random_state =10086)
    label = np.unique(dataset_meta_100.data['species'])
    
    fit = training(X_train,y_train)
    
    predict_proba = fit.predict_proba(X_test)
    pred, cm = build_confusion_matrix(y_test,predict_proba,label,True,0.9)
