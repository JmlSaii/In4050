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

#Function that finds the distance traveled, given the number of cities.
def optimaDistance(n):
    #Making a list of number of cities, ex. 4 cities: [0,1,2,3]
    listOfcities=[i for i in range(n)]
    #Making all possible perturbations based on number of cities.
    cityPermutations=list(itertools.permutations(listOfcities))
    #Listing city distance and the name of city.
    with open("european_cities.csv", "r") as f:
        data = list(csv.reader(f, delimiter=';'))
        cityNames = data[0]
        cityTravelLengths=data[1:]
    #Here all the posible distances will be listed.   
    listDistanceTraveled=[]
    for i in range(len(cityPermutations)):
        #This parameter sums all the distances.
        distanceTraveled=0
        #List of a certain(i) perturbation.
        cityPermutation=cityPermutations[i]
        for j in range(len(cityPermutation)-1):
            #Distance traveled for a certain combination of cities.
            distanceTraveled=distanceTraveled+float(cityTravelLengths[cityPermutation[j]][cityPermutation[j+1]])
        listDistanceTraveled.append(distanceTraveled)
    
    #Here we find the best distance, that is its place in the list.
    bestPermutation=listDistanceTraveled.index(min(listDistanceTraveled))
    #Shortest distance traveled.
    shortestRoute=round(listDistanceTraveled[bestPermutation],5)
    #The combination of cities that make the shortest possible route.
    shortestCombination=cityPermutations[bestPermutation]
    #The ordering of the cities.
    cityOrder=[]
    for i in range(len(shortestCombination)):
        cityOrder.append(cityNames[shortestCombination[i]])
    #Returing shotest route and the cities traveled.
    return shortestRoute,cityOrder



def exhaustiveSearch():
    totalDistance=optimaDistance(6)[0]
    plan=optimaDistance(6)[1]
    plot_plan(plan)
    plt.savefig('ExhaustiveSearchMap1.png', dpi=300, bbox_inches='tight')
    print("6 cities: ")
    print("Best sequence of cities: ",plan)
    print("Shortest distance: ",totalDistance)
    
    totalDistance=optimaDistance(10)[0]
    plan=optimaDistance(10)[1]
    plot_plan(plan)
    plt.savefig('ExhaustiveSearchMap2.png', dpi=300, bbox_inches='tight')
    print("10 cities: ")
    print("Best sequence of cities: ",plan)
    print("Shortest distance: ",totalDistance)
    
    cityCoord=[]
    bestRoute=[]
    timeCoord=[]
    for i in range(5,10):
        startTime = time.time()
        cityCoord.append(i)
        bestRoute.append(optimaDistance(i+1)[0])
        endTime=time.time()
        timeCoord.append(endTime-startTime)
        
    plt.figure()
    plt.xlabel('City (n)')
    plt.ylabel('Time (sec)')
    plt.title(('Estimated run time - Exhaustive Search'))
    plt.xlim([0, 24.3])
    print("Time for finding the best route for 6 cities: ", round(timeCoord[len(timeCoord)-5],3)," sec")
    print("Time for finding the best route for 7 cities: ", round(timeCoord[len(timeCoord)-4],3)," sec")
    print("Time for finding the best route for 8 cities: ", round(timeCoord[len(timeCoord)-3],3)," sec")
    print("Time for finding the best route for 9 cities: ", round(timeCoord[len(timeCoord)-2],3)," sec")
    print("Time for finding the best route for 10 cities: ", round(timeCoord[len(timeCoord)-1],3)," sec")
    timeCoord=np.log(timeCoord)
    C=np.polyfit(cityCoord, timeCoord, 1)
    cityCoord=np.linspace(0,24,1000)
    timeCoord=np.exp(C[0]*cityCoord+C[1])
    plt.plot(cityCoord,timeCoord,color="blue")
    print("Time to find optimal solution for 24 cities using exhaustive search algorithm: ",round((np.exp(C[0]*24+C[1]))/31536000000000,2)," million years!")
    plt.savefig('ExhaustiveSearch.png', dpi=300, bbox_inches='tight')
    plt.show() 

exhaustiveSearch()