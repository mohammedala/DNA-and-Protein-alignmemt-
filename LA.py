from Bio.SubsMat import MatrixInfo as mf
import numpy as np 

#Get score form blossom
def blossom62Score(row , col):
        Blossom62 = mf.blosum62
        try:
            val = Blossom62[row , col]
        except KeyError:
            val = Blossom62[col , row]
        return val

def LA(sequence_1 ,sequence_2, type = 0):
    # deault type is DNA
    ## initializing the matrices using numpy
    main_matrix = np.zeros((len(sequence_1)+1,len(sequence_2)+1))
    matching_matrix = np.zeros((len(sequence_1),len(sequence_2)))

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

    # filling the main matrix based on the input sequence
    for i in range(1,len(sequence_1)+1):
        for j in range(1,len(sequence_2)+1):
            main_matrix[i][j] = max(main_matrix[i-1][j-1]+matching_matrix[i-1][j-1],
                            main_matrix[i-1][j]+gap,
                            main_matrix[i][j-1]+ gap)

            if main_matrix[i][j] < 0:
                main_matrix[i][j] = 0

    #finding the largest cell
    biggest = main_matrix[0][0]
    for i in range(1,len(sequence_1)+1):
        for j in range(1,len(sequence_2)+1):
            if biggest <= main_matrix[i][j]:
                biggest = main_matrix[i][j]
                row = i
                column = j

    aligned_1 = ""
    aligned_2 = ""
    #traceback
    while(row >0 and column > 0):

        #diagonal 
        if (main_matrix[row][column] == main_matrix[row-1][column-1]+ matching_matrix[row-1][column-1]):

            aligned_1 = sequence_1[row-1] + aligned_1
            aligned_2 = sequence_2[column-1] + aligned_2

            row = row - 1
            column = column - 1
        #left
        elif(main_matrix[row][column] == main_matrix[row-1][column] + gap):
            aligned_1 = sequence_1[row-1] + aligned_1
            aligned_2 = "-" + aligned_2

            row = row -1
        #Up
        else:
            aligned_1 = "-" + aligned_1
            aligned_2 = sequence_2[column-1] + aligned_2

            column = column - 1
    
    print(aligned_1)
    print(aligned_2)

sequence_1 = "CGTGAATTCAT"
sequence_2 = "GACTTAC"
LA (sequence_2,sequence_1)



