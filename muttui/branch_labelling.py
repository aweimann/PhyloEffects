#Functions for labelling branches into spectrum categories

import operator
from Bio import Phylo

#Converts taxon labels to a dictionary with labels as keys and taxa as values
def getLabelDict(tree, labels):
    #Import the labels
    taxonLabels = open(labels).readlines()

    #Labels and keys, taxa as values
    labelDict = {}

    #Taxa that have a label, others will be assigned a label of OTHER
    analysedTaxa = []

    #Iterate through the taxa and add to labelDict
    for taxon in taxonLabels[1:]:
        if taxon.strip().split(",")[1] in labelDict:
            labelDict[taxon.strip().split(",")[1]].append(taxon.strip().split(",")[0])
        else:
            labelDict[taxon.strip().split(",")[1]] = [taxon.strip().split(",")[0]]
        analysedTaxa.append(taxon.strip().split(",")[0])
    
    #Check if OTHER has been used as a label, if it has use another OTHER label
    if "OTHER" in labelDict:
        otherLabel = "OTHER_A"
        print("OTHER has been used as a label, using OTHER_A as the label for taxa not included in the label file")
    else:
        otherLabel = "OTHER"
    
    labelDict[otherLabel] = []

    #Iterate through the taxa and add to labelDict if not already included
    for eachTaxon in tree.get_terminals():
        if eachTaxon.name not in analysedTaxa:
            labelDict[otherLabel].append(eachTaxon.name)
    
    return(labelDict)

#Adds the same label to all branches, used when there are no taxon labels provided
def labelAllBranches(tree):
    #Will contain the labels in the tree
    treeLabels = ["A"]

    #Iterate through the branches and add the label A to them
    for clade in tree.find_clades():
        clade.clade_label = "A"
    
    return(tree, treeLabels)

#Labels branches in a given tree with names as in treetime
#Branches are labelled from NODE_0000000 onwards
def labelBranchesTreetime(tree):
    #Appended with each internal node
    iterator = 0

    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = "NODE_" + str(iterator).zfill(7)
            iterator += 1
    
    return(tree)

#Writes the label dictionary to a file that can be used by treetime mugration
def writeLabels(labelDict, out_dir):
    labelFile = open(out_dir + "all_taxon_labels.csv", "w")
    labelFile.write("name,group\n")
    for eachLabel in labelDict:
        for sequence in labelDict[eachLabel]:
            labelFile.write(sequence + "," + eachLabel + "\n")
    labelFile.close()



#Label all branches in the tree with their branch name
def labelBranchesNames(tree):
    #Will contain the labels in the tree
    treeLabels = []

    #Iterate through the branches and add the label A to them
    for clade in tree.find_clades():
        clade.clade_label = clade.name
    
    return(tree, treeLabels)
