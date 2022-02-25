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
    
    for i in range(1000):
        gene = random.sample(crom,2)
        crom[gene[0]],crom[gene[1]]=crom[gene[1]],crom[gene[0]]
        tempDistance=distanceTour(crom)
        
        if tempDistance>distance:
            crom[gene[1]],crom[gene[0]]=crom[gene[0]],crom[gene[1]]
        else:
            distance=tempDistance
    
    

    return round(distance,5),crom


def printStuff(numberOfCities,numberOfRuns,picN):
    sTime=time.time()
    currentBest=swapingNeighbor(randomSearch(numberOfCities))
    overallBest=currentBest
    allRoutes=[]
    
    bestRoute=[]
    for i in range(numberOfRuns):
        allRoutes.append(currentBest[0])
        if currentBest[0]<overallBest[0]:
            overallBest=currentBest
        currentBest=swapingNeighbor(randomSearch(numberOfCities))
    eTime=time.time()
    diffTime=eTime-sTime
    
    for i in range(len(overallBest[1])):
        bestRoute.append(cityNames[overallBest[1][i]])
    
    print("------------------------------")
    print("Number of runs    |",numberOfRuns)
    print("Number of cities  |",numberOfCities)
    print("Time used         |",round(diffTime,5))
    print("Best run          |",round(np.min(allRoutes),1))
    print("Worst run         |",round(np.max(allRoutes),1))
    print("Average           |",round(np.mean(allRoutes),1))
    print("Standard deviation|",round(np.std(allRoutes),1))
    print("------------------------------")
    plot_plan(bestRoute)
    if picN==0:
        print(bestRoute)
        plt.savefig('Hillclimber1.png', dpi=300, bbox_inches='tight')
    else:
        plt.savefig('Hillclimber2.png', dpi=300, bbox_inches='tight')
        print(bestRoute)
    
def performance(optimaDistance,nCities):
    nCities=[19,20,21,22]
    optimaDistance=[9262,9888,9964,10031]
    repeat=[1,1,1,1]
    
    timeArray=[]
    cityArray=[]
    
    for j in range(len(nCities)):
        for i in range(repeat[j]):
            sTime=time.time()
            currentBest=optimaDistance[j]+1
            while optimaDistance[j]<currentBest:
                currentBest=round(swapingNeighbor(randomSearch(nCities[j]))[0],1)
            eTime=time.time()
            timeArray.append(eTime-sTime)
            cityArray.append(nCities[j])
    

    for i in range(len(timeArray)):
        if timeArray[i]>0:
            timeArray[i]=np.log(timeArray[i])
        else:
            timeArray[i]=np.log(0.00001)
        
    
    plt.figure()
    C=np.polyfit(cityArray, timeArray, 1)
    cityCoord=np.linspace(0,24,1000)
    timeCoord=np.exp(C[0]*cityCoord+C[1])
    plt.title(('Estimated run time - Hill Climbing'))
    plt.xlabel('City (n)')
    plt.ylabel('time (sec)')
    
    
    
    
    plt.plot(cityCoord,timeCoord,color="blue")
    plt.xlim([0, 24.5])
    plt.ylim([0, timeCoord[len(timeCoord)-1]])

    print("Time to find global solution for 24 cities using a simple hill climbing algorithm: ",round((np.exp(C[0]*24+C[1])),2)," seconds!")
    plt.savefig('HillClimber-fit.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    

    

    
    
        
    
def hillClimbing():
    printStuff(24,20,0)
    printStuff(10,20,1)
    print("Time to find best solutions for 10 and 24 cities:")

    
    #Will take some time (around 15 sec.)
    #The function approximates the time it takes to find the best solutions for all cities.
    performance(1,1)


    
    
    

hillClimbing()
    
    

    
