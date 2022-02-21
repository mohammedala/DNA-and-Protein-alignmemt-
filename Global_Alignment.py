import numpy as np
import pandas as pd
from Bio.SubsMat import MatrixInfo as mf


#Get score form blossom
def blossom62Score(row , col):
        Blossom62 = mf.blosum62
        try:
            val = Blossom62[row , col]
        except KeyError:
            val = Blossom62[col , row]
        return val


def GA(sequence_1, sequence_2, type = 0):   
    ## initializing the matrices using numpy
    main_matrix = np.zeros((len(sequence_1)+1,len(sequence_2)+1))
    matching_matrix = np.zeros((len(sequence_1),len(sequence_2)))

    #initializing the scoring penalties
    match = 1
    mismatch= -2
    gap = -1
    
    ## filling the matchnig matrix with comparing elmetns
    for i in range(len(sequence_1)):
        for j in range(len(sequence_2)):
            # for protien
            if type == 1:
                if sequence_1[i] == sequence_2[j]:
                    matching_matrix[i][j]= blossom62Score(sequence_1[i],sequence_2[j])
                else:
                    matching_matrix[i][j]= blossom62Score(sequence_1[i],sequence_2[j])
            elif type == 0:
                if sequence_1[i] == sequence_2[j]:
                    matching_matrix[i][j]= match
                else:
                    matching_matrix[i][j]= mismatch
                
   ## creating the first row and coloumn of the main matrix
    for i in range(len(sequence_1)+1):
        main_matrix[i][0] = i*gap
    for j in range(len(sequence_2)+1):
        main_matrix[0][j] = j*gap

    ## filling the main matrix based on the input sequence
    for i in range(1,len(sequence_1)+1):
        for j in range(1,len(sequence_2)+1):
            main_matrix[i][j] = max(main_matrix[i-1][j-1]+matching_matrix[i-1][j-1],
                                main_matrix[i-1][j]+gap,
                                main_matrix[i][j-1]+ gap)
            
   ## Creating the traceback
    aligned_1 = ""
    aligned_2 = ""

    row = len(sequence_1)
    col = len(sequence_2)

    while(row >0 and col > 0):

        if (main_matrix[row][col] == main_matrix[row-1][col-1]+ matching_matrix[row-1][col-1]):

            aligned_1 = sequence_1[row-1] + aligned_1
            aligned_2 = sequence_2[col-1] + aligned_2

            row = row - 1
            col = col - 1
    
        elif(main_matrix[row][col] == main_matrix[row-1][col] + gap):
            aligned_1 = sequence_1[row-1] + aligned_1
            aligned_2 = "-" + aligned_2

            row = row -1
        else:
            aligned_1 = "-" + aligned_1
            aligned_2 = sequence_2[col-1] + aligned_2

            col = col - 1
    print(aligned_1)
    print(aligned_2)



GA("ATCGT", "ACGT")

