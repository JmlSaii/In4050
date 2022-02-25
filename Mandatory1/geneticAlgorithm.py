from europeanCities import*
import itertools
import csv
import time
import matplotlib.pyplot as plt
import random
from scipy.optimize import curve_fit
from IPython.display import display
import pandas as pd
import numpy as np


with open("european_cities.csv", "r") as f:
        data = list(csv.reader(f, delimiter=';'))
        cityNames = data[0]
        cityTravelLengths=data[1:]

def randomSearch(n):
    randomConfig=list(range(n))
    randomConfig=random.sample(randomConfig, len(randomConfig)) 
    return randomConfig

def distanceTour(cityDistance):
    totalDistance=0
    for i in range(len(cityDistance)-1):
        totalDistance=totalDistance+float(cityTravelLengths[cityDistance[i]][cityDistance[i+1]])
    return round(totalDistance,3)
    

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

def swapMutation(parent):
    Chromosome=parent    
    n=len(parent)
    gene1=random.randrange(n)
    gene2=random.randrange(n)
    Chromosome[gene1],Chromosome[gene2]=Chromosome[gene2],Chromosome[gene1]
    return Chromosome

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
    n=len(parent1)
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
                
#this is an incomplete crossover operator            
def edgeCrossover(parents):
    n=10
    parent1=randomSearch(n)
    parent1.insert(0,parent1[len(parent1)-1])
    parent1.insert(len(parent1),parent1[1])
    
    parent2=randomSearch(n)
    
    parent2.insert(0,parent2[len(parent2)-1])
    parent2.insert(len(parent2),parent2[1])
    
    edgeList=[0]*(len(parent1)-2)

    for i in range(1,(len(parent1)-1)):
        neighbors=[]
        neighbors.append(parent1[i-1])
        neighbors.append(parent1[i+1])
        edgeList[parent1[i]]=(neighbors)
        pos=parent2[1:(len(parent1)-1)].index(parent1[1:(len(parent1)-1)][i-1])+1
        
        if parent2[pos-1] not in edgeList[parent1[i]]:
            neighbors.append(parent2[pos-1])
        if parent2[pos+1] not in edgeList[parent1[i]]:
            neighbors.append(parent2[pos+1])
    
    print(edgeList)
    child=[]
    for i in range(len(edgeList)):

        lenList=[]
        if i==0:
            checkIndex=0
        
        
        for j in range(len(edgeList[checkIndex])):
            lenList.append(len(edgeList[edgeList[checkIndex][j]]))
        
       
        child.append(checkIndex)
        
        tempIndex=checkIndex
        print(edgeList)
        print(edgeList[tempIndex])
        checkIndex=edgeList[tempIndex][lenList.index(max(lenList))]
        print(edgeList[checkIndex])
        
     
        k=0
        for k in range(len(edgeList)):
            trueFalse=tempIndex in edgeList[k]
            
            if trueFalse==True:
                edgeList[k].remove(tempIndex)
                
            k+=1

    return edgeList

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
    
    
def recombinationOperators(parent1,parent2):
    p=random.random()  
    if 0<p<0.1:    
        return cycleCrossOver(parent1,parent2)  #worst    
    if 0.1<p<0.2:
        return PMX(parent1,parent2)             #ok      
    if 0.2<p<1:        
        return orderCrossover(parent1,parent2)  #best                         
    if p<0:
        return edgeCrossover(parent1,parent2)


    
def selection(nPop,pop_list):
    selectPop=[-1]*2*nPop
    for k in range(nPop):
        if distanceTour(pop_list[nPop+k])<distanceTour(pop_list[k]):
            selectPop[k]=pop_list[nPop+k]
        else:
            selectPop[k]=pop_list[k]   
    
  
    return selectPop    
    
    
#choosing which operators to use
def operators(nPop,pop_list):
    counter =0
    while counter <nPop:
        mutChance=random.random()
        if mutChance <0.01:
            pop_list[nPop+counter]=mutationOperators(pop_list[counter])
            pop_list[nPop+counter+1]=mutationOperators(pop_list[counter+1])
            counter=counter+2
        else:
            mutChild=recombinationOperators(pop_list[counter],pop_list[counter+1])
            pop_list[nPop+counter]=mutationOperators(mutChild[0])
            pop_list[nPop+counter+1]=mutationOperators(mutChild[1])
            counter=counter+2
    
    return pop_list

#choosing inital population based on number of cities and population
def initialPop(nPop,nCities):
    pop_list=[-1]*nPop*2
    for j in range(nPop):
        pop_list[j]=randomSearch(nCities)
    
    return pop_list

#For plotting and printing
def plotPrintStuff(Gen,Fitness,nCities,generation,diffTime,nPop,listOfCities,nMap):
    s=[]
    for i in range(nPop):
        s.append(distanceTour(listOfCities[i]))
    
    routePlan=listOfCities[s.index(min(s))]
    
        
    plan=[]
    for i in range(len(routePlan)):
        plan.append(cityNames[routePlan[i]])
    plot_plan(plan)
    
    if nMap==0:
        plt.savefig('geneticMap4.png', dpi=300, bbox_inches='tight')
    if nMap==1:
        plt.savefig('geneticMap5.png', dpi=300, bbox_inches='tight')
    if nMap==2:
        plt.savefig('geneticMap6.png', dpi=300, bbox_inches='tight')
    print(plan)
    print("------------------------------")
    print("Population        |",nPop)
    print("Number of runs    |",generation)
    print("Tours inspected   |",2*generation*nPop)
    print("Number of cities  |",nCities)
    print("Time used         |",round(diffTime,5))
    print("Best run          |",round(np.min(s),1))
    print("Worst run         |",round(np.max(s),1))
    print("Average           |",round(np.mean(s),1))
    print("Standard deviation|",round(np.std(s),1))
    print("------------------------------")

#Gives best parent in the population
def currentBestGen(nPop,pop):
    s=[]
    for i in range(nPop):
        s.append(distanceTour(pop[i]))
    return min(s)
    
    
    
    
def GA(nPop,nCities,generation):
    stime=time.time()
    pop_list=initialPop(nPop,nCities)
    Gen=[]
    Fitness=[]
    iteration=0
    while iteration<generation:
        Gen.append(iteration)
        Fitness.append(currentBestGen(nPop,pop_list[0:nPop]))
        pop_list=selection(nPop,operators(nPop,pop_list))
        iteration=iteration+1
    etime=time.time()
    diffTime=etime-stime
    return Gen,Fitness,nCities,generation,diffTime,nPop,pop_list[0:nPop]
    

    

    
def mainGA():
    nPop=[10,100,1000]
    nCities=24
    generation=20
    par=[]
    coords=[]
    for i in range(3):
        par.append(GA(nPop[i],nCities,generation))
        plotPrintStuff(par[i][0],par[i][1],par[i][2],par[i][3],par[i][4],par[i][5],par[i][6],i)
        coords.append([par[i][0],par[i][1]])
        
    plt.figure()
    plt.xlabel('Generations')
    plt.ylabel('Distance')

    
    
    
    c=["green","blue","red"]
    
    for i in range(len(coords)):
        x=coords[i][0]
        y=coords[i][1]
        plt.plot(x,y,color=c[i])
        
        
        
        
    plt.legend(['population=10', 'population=100', 'population=1000'])
    plt.savefig('genetic.png', dpi=300, bbox_inches='tight')
    

mainGA()

#This function plots time it takes to find best solution for a given number of cities. 
def findingOptimaltime():
    nPop=20
    nCities=[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
    bestDist=[3168,3637,4817,4828,5273,6326,6334,6820,7375,8132,
              8386,8635,8697,9262,9888,9964,10031,10055,10838]
    Time=[]
    for j in range(3):
        for i in range(len(nCities)):
            pop_list=initialPop(nPop,nCities[i])
            stime=time.time()
            fitness=currentBestGen(nPop,pop_list)
            while bestDist[i]<fitness:
                fitness=currentBestGen(nPop,pop_list[0:nPop])
                pop_list=selection(nPop,operators(nPop,pop_list))
           
            
            etime=time.time()
            Time.append(etime-stime)
    
    averageTime=[]
    
    for i in range(len(nCities)):
        averageTime.append((Time[i]+Time[len(nCities)+i]+Time[2*len(nCities)+i])/3)
        
    plt.figure()
    plt.title(('Estimated runtime - pure GA'))
    plt.scatter(nCities,averageTime,s=30,marker ="x",color="red")
    plt.xlabel('City (n)')
    plt.ylabel('Time (sec)')
    x=nCities
    y=np.log(averageTime)
    C=np.polyfit(x,y,1)
    x=np.linspace(0,24,1000)
    y=np.exp(C[0]*x+C[1])
    plt.plot(x,y,color="blue")
    
    
    plt.savefig('geneticTime.png', dpi=300, bbox_inches='tight')

    
#This may take some time.    
#findingOptimaltime()    


def mainGA():
    nCities=[10,24]
    nPop=[20,20]
    generation=[200,10000]
    par=[]
    textTime=["  / 10 sec < t","   / 1 mill year < t"]
    textDis=["/ 5272.7","/ 10837.1 "]
    for i in range(2):
        par.append(GA(nPop[i],nCities[i],generation[i]))
        print("------------------------------")
        print("Tours inspected   |",generation[i]*nPop[i],"  /",np.math.factorial(nCities[i]))
        print("number of runs    |",generation[i],"  /")
        print("Number of cities  |",nCities[i], "     /",nCities[i])
        print("Time used         |",round(par[i][4],3),textTime[i])
        print("Best run          |",round(min(par[i][1]),1),textDis[i])
        print("------------------------------")
    print("About 0.02 seconds to find the best solution for 10 cities (at least 100 times faster), and around 10 seconds for all 24 cities. ") 

    s=[]
    for i in range(20):
        s.append(distanceTour(par[1][6][i]))
    
    routePlan=par[1][6][s.index(min(s))]
    
    
    plan=[]
    for i in range(len(routePlan)):
        plan.append(cityNames[routePlan[i]])
    print("Best Plan: ") 
    print(plan) 
    plot_plan(plan)
    plt.savefig('geneticMap7.png', dpi=300, bbox_inches='tight')
    
    
mainGA()