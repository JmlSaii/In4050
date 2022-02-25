#Using all the functions from pure GA

from GAHLfunctions import*

#modified operator function.
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

#Selection population
def hybridSelection(popList,mode):  
    #Saving the parents
    parents=random.sample(popList,int(len(popList)/10))
    parents=np.array(parents)
    #Saving the fitted parents and their distances
    fitedParents=parents.tolist()
    distanceFP=[]
    for i in range(len(fitedParents)):
        sW=swapingNeighbor(fitedParents[i])
        fitedParents[i]=sW[1]
        distanceFP.append(sW[0])
    parents=parents.tolist()
    
    #2 of the best fitted parents will be choosen
    if mode=="Lamarck":
        gene=distanceFP.index(min(distanceFP))
        parent=fitedParents[gene]
        fitedParents.remove(parent)
        distanceFP.remove(min(distanceFP))
        

    #2 parents will be choosen that are potentially best fitted.
    if mode=="Baldwin":
        gene=distanceFP.index(min(distanceFP))
        parent=parents[gene]
        parents.remove(parent)
        distanceFP.remove(min(distanceFP))

    return parent

#Selection of new population.
def eliteSelection(nPop,popList): 
    selectPop=[-1]*2*nPop
    for k in range(nPop):
        if distanceTour(popList[nPop+k])<distanceTour(popList[k]):
            selectPop[k]=popList[nPop+k]
        else:
            selectPop[k]=popList[k]   
    return selectPop

#Choosing intial population
def initialPop(nPop,nCities,mode):
    popList=[-1]*nPop*2
    for i in range(nPop):
            popList[i]=randomSearch(nCities)
    return popList


#Main function. From here everything is controlled.
def hybridGA(nPop,nCities,generation,mode):
    stime=time.time()
    popList=initialPop(nPop,nCities,mode)
    Gen=[]
    Fitness=[]
    iteration=0
    while iteration<generation:
        Gen.append(iteration)
        Fitness.append(currentBestGen(nPop,popList[0:nPop]))
        popList=eliteSelection(nPop,operators(nPop,popList,mode))
        iteration=iteration+1
    etime=time.time()
    diffTime=etime-stime
    return Gen,Fitness,nCities,generation,diffTime,nPop,popList[0:nPop]

#Calculating times for the 2 diffrent algorithms
def hybridTours(mode):
    nPop=20
    nCities=[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
    bestDist=[3168,3637,4817,4828,5273,6326,6334,6820,7375,8132,
          8386,8635,8697,9262,9888,9964,10031,10055,10838]
    Time=[]
    val=1
    for j in range(val):
        for i in range(len(nCities)):
            popList=initialPop(nPop,nCities[i],mode)
            stime=time.time()
            fitness=currentBestGen(nPop,popList)
            while bestDist[i]<fitness:
                fitness=currentBestGen(nPop,popList[0:nPop])
                popList=eliteSelection(nPop,operators(nPop,popList,mode))
           
            etime=time.time()
            Time.append(etime-stime)

    averageTime=[]
    for i in range(len(nCities)):
        tot=0
        for j in range(val):
            tot=tot+ Time[j*len(nCities)+i]
        tot=tot/val
        averageTime.append(tot)
    
    for i in range(len(averageTime)):
        if averageTime[i]==0:
            averageTime[i]=0.000001

    x=nCities     
    y=np.log(averageTime)
    C=np.polyfit(x,y,1)
    x=np.linspace(0,25,1000)
    y=np.exp(C[0]*x+C[1])
    return x,y

    
#Plotting the 3 diffrent algorithms
def compareTimePlot():
    plt.figure()
    plt.title(('Hybrid GA'))
    plt.xlabel('City (n)')
    plt.ylabel('Time (sec)')
    
    modes=["Lamarck","Baldwin"]
    Color=["blue","red"]
    for i in range(3):
        coord=hybridTours(modes[i])
        plt.plot(coord[0],coord[1],color=Color[i])
        
    plt.legend(['Lamarckian', 'Baldwinian'])
    plt.savefig('hybridTimes.png', dpi=300, bbox_inches='tight')
    

   
#plotting and printing table statstics.
def compareTablePlot():
    plt.figure()
    modes=["Lamarck","Baldwin"]
    Color=["blue","red","green"]
    parm=[]
    for i in range(2):
        parm.append(hybridGA(100,10,100,modes[i]))
               
        plt.plot(parm[i][0],parm[i][1],color=Color[i])
        
    plt.legend(['Lamarckian', 'Baldwinian', 'Pure GA'])
    plt.xlabel("Generations")
    plt.ylabel("Distance")
    plt.title("Hybrid GA - 10 cities")
    plt.savefig('hybridGen10.png', dpi=300, bbox_inches='tight')
    
    lastGenB=[]
    for i in range(len(parm[1][6])):
        lastGenB.append(distanceTour(parm[1][6][i]))
    
    lastGenL=[]
    for i in range(len(parm[0][6])):
        lastGenL.append(distanceTour(parm[0][6][i]))
        

    tours=100*100*100+10*100
    #stats for 10 cities
    data = {'Name':['Number of cities', 'Generations', 'Populaion', 'Routes ',
                    'time (s)','Best distance', 'Worst distance', 'Average distance', 
                    'Standard deviation'],
        
        'Baldwinian':[parm[0][2], parm[1][3], parm[0][5], tours,round(parm[1][4],2), round(min(lastGenB),2), round(max(lastGenB),2), round(np.mean(lastGenB),2),round(np.std(lastGenB),2)],
        'Lamarckian':[parm[0][2], parm[0][3], parm[0][5], tours,round(parm[0][4],2), round(min(lastGenL),2), round(max(lastGenL),2), round(np.mean(lastGenL),2),round(np.std(lastGenL),2)]}
    df = pd.DataFrame(data)
    parm=[]
    print(df)
    
    plt.figure()
    for i in range(2):
        parm.append(hybridGA(100,24,100,modes[i]))
               
        plt.plot(parm[i][0],parm[i][1],color=Color[i])
        
    plt.legend(['Lamarckian', 'Baldwinian', 'Pure GA'])
    plt.xlabel("Generations")
    plt.ylabel("Distance")
    plt.title("Hybrid GA - 24 cities")
    plt.savefig('hybridGen20.png', dpi=300, bbox_inches='tight')

    lastGenB=[]
    for i in range(len(parm[1][6])):
        lastGenB.append(distanceTour(parm[1][6][i]))
    
    lastGenL=[]
    for i in range(len(parm[0][6])):
        lastGenL.append(distanceTour(parm[0][6][i]))
    #stats for 24 cities
    data = {'Name':['Number of cities', 'Generations', 'Populaion', 'Routes ',
                    'time (s)','Best distance', 'Worst distance', 'Average distance', 
                    'Standard deviation'],
        
        'Baldwinian':[parm[0][2], parm[1][3], parm[0][5], tours,round(parm[1][4],2), round(min(lastGenB),2), round(max(lastGenB),2), round(np.mean(lastGenB),2),round(np.std(lastGenB),2)],
        'Lamarckian':[parm[0][2], parm[0][3], parm[0][5], tours,round(parm[0][4],2), round(min(lastGenL),2), round(max(lastGenL),2), round(np.mean(lastGenL),2),round(np.std(lastGenL),2)]}
    df = pd.DataFrame(data)
    print("------------")
    print(df)


                                                                        
def mainHybrid():
    compareTimePlot()
    compareTablePlot()
    
#compareTimePlot()    
compareTablePlot()     
        
        
