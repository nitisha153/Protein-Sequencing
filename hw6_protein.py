"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f = open(filename,'r')
    text = f.read()
    f.close()

    text = text.splitlines()
    #print(text)
    text_1 = ""
    for line in text:
        text_1 += line
    
    
    return text_1
#print(readFile(("data/human_p53.txt"))


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    rna = ""
    for char in dna:
        if char == "T":
            rna += "U"
        else:
            rna += char
    stop_codons = ["UAA","UAG","UGA"]
    rna_list = []
    i = (len(dna)-startIndex)%3
    j = len(dna) - i
    for index in  range(startIndex,j,3):
        codon = rna[index]+rna[index+1]+rna[index+2]
        rna_list.append(codon)
        if codon in stop_codons:
            break



        

    return rna_list


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    data = open(filename,'r')
    dictionary = json.load(data)
    data.close()
    new_dict ={}
    for element in dictionary:
        codon_list = dictionary[element]
        for codon in codon_list:
            rna =""
            for char in codon:
                if char == "T":
                    rna += "U"
                else:
                    rna += char
            new_dict[rna] = element


    #print(new_dict)
    return new_dict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    proteins = []
    for index in range(len(codons)):
        if index == 0:
            proteins.append("Start")
        else:
            proteins.append(codonD[codons[index]])
    
    #print(proteins)
    return proteins


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna = readFile(dnaFilename)
    codonD = makeCodonDictionary(codonFilename)
    proteins = []
    i = 0
    while i < len(dna):
        dna_string = dna[i:]
        start_index = dna_string.find("ATG")
        if start_index < 0:
            break
        rna_list = dnaToRna(dna_string,start_index)
        protein = generateProtein(rna_list,codonD)
        proteins.append(protein)
        i += start_index + len(protein)*3
    #print(proteins) 
    return proteins


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    common_proteins = []
    for ele_1 in proteinList1:
        if ele_1 in proteinList2:
            if ele_1 not in common_proteins:
                common_proteins.append(ele_1)
            

    #print(common_proteins)

    return common_proteins


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    amino_acids = []
    for element in proteinList:
        for sub_element in element:
            amino_acids.append(sub_element)
    
    #print(amino_acids)

    return amino_acids


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aa_dict = {}
    for aa in aaList:
        if aa not in aa_dict:
            aa_dict[aa] = 0
        aa_dict[aa] += 1
        
        
    #print(aa_dict)
    return aa_dict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    a_acid1 = combineProteins(proteinList1)
    a_acid2 = combineProteins(proteinList2)
    aa_dict1 = aminoAcidDictionary(a_acid1)
    aa_dict2 = aminoAcidDictionary(a_acid2)
    freq1 = {}
    freq2 = {}
    for acid in aa_dict1:
        freq1[acid] = aa_dict1[acid]/len(a_acid1)
    for acid in aa_dict2:
        freq2[acid] = aa_dict2[acid]/len(a_acid2)
    amino_acids = [] #without Start and Stop amino acids in the list
    for aa in aa_dict1:
        #if aa != "Start" and aa != "Stop":
            if aa not in amino_acids and aa != "Start" and aa != "Stop":
                amino_acids.append(aa)
    for aa in aa_dict2:
        #if aa != "Start" and aa != "Stop":
            if aa not in amino_acids and aa != "Start" and aa != "Stop":
                amino_acids.append(aa)
    aa_difference = []
    for acid in amino_acids:
        f1 = 0
        f2 = 0
        if acid in freq1:
            f1 = freq1[acid]
        if acid in freq2:
            f2 = freq2[acid]
        diff = f2 - f1
        if diff > cutoff or diff < -cutoff:
            aa_difference.append([acid,f1,f2])
        
        


    

    return aa_difference


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The following proteins occurred in both DNA Sequences:")
    for proteins in commonalities:
        common_protein = ""
        protein = proteins[1:len(proteins)-1]
        count = 0
        for element in protein:
            common_protein += element
            count += 1
            if count != len(protein):
                common_protein += "-"
        if len(common_protein) != 0:
            print(common_protein)
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for amino_acid in differences:
        acid = amino_acid[0]
        f1 = round(amino_acid[1]*100,2)
        f2 = round(amino_acid[2]*100,2)

        print(acid + " : " + str(f1) + "% in Seq1," + str(f2) + "% in Seq2")

    return None


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    a_acid1 = combineProteins(proteinList1)
    a_acid2 = combineProteins(proteinList2)
    aa_dict1 = aminoAcidDictionary(a_acid1)
    aa_dict2 = aminoAcidDictionary(a_acid2)
    labels = []
    for amino_acid in aa_dict1:
        if amino_acid not in labels:
            labels.append(amino_acid)
    for amino_acid in aa_dict2:
        if amino_acid not in labels:
            labels.append(amino_acid)
    labels.sort()
    return labels


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    proteins = combineProteins(proteinList)
    amino_acid_dict = aminoAcidDictionary(proteins)
    frequency = []
    for label in labels:
        f = 0
        if label not in amino_acid_dict:
            frequency.append(f)
        else:
            f = amino_acid_dict[label]/len(proteins)
            frequency.append(f)



    

    return frequency


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    import numpy as np
    # reference link: https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py
    x = np.arange(len(freqList1))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, freqList1, width, label=label1, edgecolor = edgeList)
    rects2 = ax.bar(x + width/2, freqList2, width, label=label2, edgecolor = edgeList)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Frequency')
    ax.set_title('Frequency comparison of Genes')
    ax.set_xticks(x)
    ax.set_xticklabels(xLabels)
    ax.legend()

    #ax.bar_label(rects1, padding=3)
    #ax.bar_label(rects2, padding=3)

    fig.tight_layout()

    plt.show()
    return None


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    new_list = []
    for element in biggestDiffs:
        #print(element)
        new_list.append(element[0])
    edge_colour = []
    for label in labels:
        if label in new_list:
            edge_colour.append("black")
        else:
            edge_colour.append("white")
    

    return edge_colour


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    human_dna = "data/human_p53.txt"
    elephant_dna = "data/elephant_p53.txt"

    human_proteins = synthesizeProteins(human_dna,"data/codon_table.json")
    elephant_proteins = synthesizeProteins(elephant_dna,"data/codon_table.json")
    
    commons = commonProteins(human_proteins,elephant_proteins)

    differences = findAminoAcidDifferences(human_proteins,elephant_proteins,0.005)

    displayTextResults(commons,differences)

    labels = makeAminoAcidLabels(human_proteins,elephant_proteins)
    f1 = setupChartData(labels,human_proteins)
    f2 = setupChartData(labels,elephant_proteins)
    edges = makeEdgeList(labels,differences)
    createChart(labels,f1,"Human",f2,"Elephant",edgeList=edges)
    return None



### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()

    ## Uncomment these for Week 2 ##

    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    

    ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
