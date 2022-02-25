#This script is from the pure GA and are used for the hybrid GA so that I dont have
#to write them again.
import random
from europeanCities import*
import itertools
import csv
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from IPython.display import display
import pandas as pd
import numpy as np
with open("european_cities.csv", "r") as f:
    data = list(csv.reader(f, delimiter=';'))
    cityNames = data[0]
    cityTravelLengths=data[1:]

#creating a random route.
def randomSearch(n):
    randomConfig=list(range(n))
    randomConfig=random.sample(randomConfig, len(randomConfig)) 
    return randomConfig  

#Calculating the distance of a route
def distanceTour(cityDistance):
    totalDistance=0
    for i in range(len(cityDistance)-1):
        totalDistance=totalDistance+float(cityTravelLengths[cityDistance[i]][cityDistance[i+1]])
    return round(totalDistance,3)

#Swaping neighbors and returning the best.   
def swapingNeighbor(randomSearch):  
    crom=randomSearch
    distance=distanceTour(crom)
    for i in range(10):
        gene = random.sample(crom,2)
        crom[gene[0]],crom[gene[1]]=crom[gene[1]],crom[gene[0]]
        tempDistance=distanceTour(crom)
        if tempDistance>distance:
            crom[gene[1]],crom[gene[0]]=crom[gene[0]],crom[gene[1]]
        else:
            distance=tempDistance
    return round(distance,5),crom

#Mutating a parent, exchaning 2 genes.
def swapMutation(parent):
    Chromosome=parent    
    n=len(parent)
    gene1=random.randrange(n)
    gene2=random.randrange(n)
    Chromosome[gene1],Chromosome[gene2]=Chromosome[gene2],Chromosome[gene1]
    return Chromosome

#Randomizing genes in the parent.
def scrambleMutation(parent):
    Chromosome=parent
    n=len(parent)
    nGenes=random.randrange(2,n)
    genePositions=random.sample(range(0, n), nGenes)
    changePositions=random.sample(genePositions, len(genePositions))
    while len(genePositions)>1:
        allele1=genePositions[0]
        allele2=changePositions[0]
        Chromosome[allele1],Chromosome[allele2]=Chromosome[allele2],Chromosome[allele1]
        if genePositions[0]==changePositions[0]:
            del(genePositions[0])
            del(changePositions[0])
        else:
            del(genePositions[genePositions.index(allele1)])
            del(genePositions[genePositions.index(allele2)])
            del(changePositions[changePositions.index(allele2)])
            del(changePositions[changePositions.index(allele1)])
    return Chromosome


def inversionMutation(parent):
    Chromosome=parent
    n=len(parent)
    p=random.sample(range(0, n), 2)
    p.sort()
    tempList=Chromosome[p[0]:p[1]+1]
    mutatedChromosome=[]
    for i in range(len(tempList)):
        mutatedChromosome.append(tempList[len(tempList)-i-1])
    mutatedChromosome=Chromosome[0:p[0]]+mutatedChromosome+Chromosome[p[1]+1:]
    return mutatedChromosome

def insertMutation(parent):
    Chromosome=parent
    n=len(parent)
    genes=random.sample(range(0, n), 2)
    if genes[0]==n-1:
        Chromosome.insert(genes[0],Chromosome[genes[1]])
        del(Chromosome[genes[1]])
    elif genes[0]>genes[1]:
        Chromosome.insert(genes[0]+1,Chromosome[genes[1]])
        del(Chromosome[genes[1]])
    else:
        Chromosome.insert(genes[0]+1,Chromosome[genes[1]])
        del(Chromosome[genes[1]+1]) 
    return Chromosome

    
    
def orderCrossover(parent1,parent2):
    n=len(parent1)
    child1=[-1]*n
    p=random.sample(range(0, n), 2)
    p.sort()
    child1[p[0]:p[1]+1]=parent1[p[0]:p[1]+1]
    remainingGenes=[]
    for i in range(len(parent2)):
        if parent2[i] not in child1:
            remainingGenes.append(parent2[i])
    child1=remainingGenes[n-p[1]-1:]+parent1[p[0]:p[1]+1]+remainingGenes[0:n-p[1]-1]
    parent1,parent2=parent2,parent1
    child2=[-1]*n
    p=random.sample(range(0, n), 2)
    p.sort()
    child2[p[0]:p[1]+1]=parent1[p[0]:p[1]+1]
    remainingGenes=[]
    for i in range(len(parent2)):
        if parent2[i] not in child2:
            remainingGenes.append(parent2[i])
    child2=remainingGenes[n-p[1]-1:]+parent1[p[0]:p[1]+1]+remainingGenes[0:n-p[1]-1]
    return child1,child2


def cycleCrossOver(parent1,parent2):
    p1,p2=parent1,parent2
    child1=[-1]*len(parent1)
    while -1 in child1:
        element = child1.index(-1)
        l1 = []
        l2 = []
        while element not in l1:
            certainElement = parent1[element]
            l1.append(element)
            l2.append(certainElement)
            element = parent1.index(parent2[element])
        for element,certainElement in zip(l1, l2):
            child1[element] = certainElement
        parent1,parent2 = parent2,parent1
    parent1,parent2=p2,p1
    child2=[-1]*len(parent1)
    while -1 in child2:
        element = child2.index(-1)
        l1 = []
        l2 = []
        while element not in l1:
            certainElement = parent1[element]
            l1.append(element)
            l2.append(certainElement)
            element = parent1.index(parent2[element])
        for element,certainElement in zip(l1, l2):
            child2[element] = certainElement
        parent1,parent2 = parent2,parent1
    return child1,child2



def PMX(parent1,parent2):
    n=len(parent1)
    p=random.sample(range(0, n), 2)
    p.sort()
    child = [-1]*len(parent1)
    child[p[0]:p[1]] = parent1[p[0]:p[1]]
    for allele in range(len(parent2[p[0]:p[1]])):
        failSafe=0
        gene=(parent2[p[0]:p[1]])[allele]
        if gene not in child:
            while child[allele] != -1:
                failSafe+=1
                allele = parent2.index(parent1[allele])
                if failSafe==len(parent1):
                    break
            child[allele] = gene
    for i in range(len(child)):
        if child[i] == -1:
            for j in range(len(parent2)):
                if parent2[j] not in child:
                    child[i] = parent2[j]
                    
    child1=child
    parent1,parent2=parent2,parent1
    p=random.sample(range(0, n), 2)
    p.sort()
    child = [-1]*len(parent1)
    child[p[0]:p[1]] = parent1[p[0]:p[1]]
    for allele in range(len(parent2[p[0]:p[1]])):
        failSafe=0
        gene=(parent2[p[0]:p[1]])[allele]
        if gene not in child:
            while child[allele] != -1:
                failSafe+=1
                allele = parent2.index(parent1[allele])
                if failSafe==len(parent1):
                    break
            child[allele] = gene    
    for i in range(len(child)):
        if child[i] == -1:
            for j in range(len(parent2)):
                if parent2[j] not in child:
                    child[i] = parent2[j]
    child2=child
    return child1, child2

#choosing a certain mutation based on the probabilty.
def mutationOperators(parent):
    p=random.random()   
    if 0<p<0.95:
        return inversionMutation(parent)   #best
    if 0.94<p<0.97:
        return insertMutation(parent)      #ok
    if 0.97<p<0.99:
        return swapMutation(parent)        #bad
    if 0.99<p:
        return scrambleMutation(parent)    #worst
    
#choosing a recombination operator.    
def recombinationOperators(parent1,parent2):
    p=random.random()  
    if 0<p<0.1:    
        return cycleCrossOver(parent1,parent2)  #worst    
    if 0.1<p<0.2:
        return PMX(parent1,parent2)             #ok      
    if 0.2<p<1:        
        return orderCrossover(parent1,parent2)  #best                         


#choosing what sort of mutation the parents goes through.
def operators(nPop,popList,mode):
    counter =0
    while counter <nPop:
        mutChance=random.random()
        if mutChance <0.01:
            popList[nPop+counter]=mutationOperators(hybridSelection(popList[0:nPop],mode))
            popList[nPop+counter+1]=mutationOperators(hybridSelection(popList[0:nPop],mode))
            counter=counter+2
        else:
            mutChild=recombinationOperators(hybridSelection(popList[0:nPop],mode),hybridSelection(popList[0:nPop],mode))
            popList[nPop+counter]=mutationOperators(mutChild[0])
            popList[nPop+counter+1]=mutationOperators(mutChild[1])
            counter=counter+2
    return popList

#Gives back best route.
def currentBestGen(nPop,pop):
    s=[]
    for i in range(nPop):
        s.append(distanceTour(pop[i]))
    return min(s)
