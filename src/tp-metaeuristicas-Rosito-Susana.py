import sys
from math import sqrt
import collections
import copy

graphGuie = None
criterioParadaBusqLocal = 1000
criterioParadaTabuSearch = 1000
tamListaTabu = 40

class Node:
    def __init__(self, idty, frequency=None):
        self.id = idty
        self.freq = frequency
        self.indexfreq = None
        self.conflictsToNeigh = {}
    
    def assignFreq(self, frequency, index):
        self.freq = frequency         
        self.indexfreq = index

    def isFreqNone(self):
        return self.freq == None

    def addConfNeigh(self, neigh, conf):
        self.conflictsToNeigh[neigh].insert(1, conf)

    def getSumNewConflictNeighs(self, idxfreq):
        conf = 0 
        for keyN in self.conflictsToNeigh.keys():
            if self.getNeigh(keyN)[0].indexfreq == idxfreq:
                conf += graphGuie.getVtx(self.id).getConflictNeigh(self.getNeigh(keyN)[0].id)
        return conf

    def getSumConflictNeighs(self):
        conf = 0 
        for keyN in self.conflictsToNeigh.keys():
            conf += self.getConflictNeigh(keyN)
        return conf

    def getConflictNeigh(self, neigh):
        return self.getNeigh(neigh)[1]

    def getNeigh(self, neigh):
        return self.conflictsToNeigh[neigh]
    
    def getNeighs(self):
        return self.conflictsToNeigh
            
    def addNeigh(self, idtyN, neighb, conf):
         if not (idtyN in self.conflictsToNeigh.keys()):
            self.conflictsToNeigh[idtyN] = [neighb, conf]
            
    def __str__(self):
        res = "Node id: %s freq: %s  \n" % (self.id, self.indexfreq) 
        res += "Neighbours and Conflict: \n"
        for node in self.conflictsToNeigh.keys():
            res += str(self.getNeigh(node)[0].id) + " assign Freq: " + str(self.getNeigh(node)[0].indexfreq) + " " + str(self.conflictsToNeigh[node][1]) + " \n"
        return res 
        
class Graph:
    def __init__(self, freqs, countUsedFreqs):
        self.frequency = freqs
        self.countUsedFrequency = countUsedFreqs
        self.costTotalSol = 0
        self.cotaCostUsedFreq = 0
        self.costHastaAhoraUsedNewFreqTotal = 0
        self.idxFreq = 0
        self.listAdj = {}
    
    #------- Metodos del Grafo --------

    # Agrega un vertice al grado    
    def addVtx(self, vtx): 
        vid = vtx.id
        if vid not in self.getVtxs():
            self.listAdj[vid] = vtx

    # Devuelve vertice a partir de id
    def getVtx(self, idt): 
        return self.listAdj[idt]

    # Devuelve todos los vertices del grafo
    def getVtxs(self): 
        return list(self.listAdj.keys())

    # Agrega el vecino edge al vertice vertex con conflicto conf 
    def addEdge(self, vertex, edge, conf):
        self.getVtx(vertex).addNeigh(edge, self.getVtx(edge), conf)

    # devuelve los vecinos del vertice vertex
    def getNeighbours(self, vertex):
        return self.getVtx(vertex).getNeighs()

    def getCostTotalSol(self):
        return self.costTotalSol

    #------- Abajo - Metodos relacionados al problema --------
    
    # cambia la frecuencia de un vertice
    def changeFreqToVertex(self, vertex, idxfreq):
        changeVert = self.getVtx(vertex) 
        freq = self.frequency[idxfreq]
        changeVert.assignFreq(freq, idxfreq)
        
        for neigh in changeVert.getNeighs():
            vecConf = changeVert.getNeigh(neigh)
            idxVecFreq = vecConf[0].indexfreq
            if(idxVecFreq == idxfreq):
                if(vecConf[1] == 0):
                    print("misma frecuencia")
                    costVec = graphGuie.getVtx(vertex).getConflictNeigh(vecConf[0].id) 
                    print(vecConf[0].id)
                    print(costVec)
                    changeVert.addConfNeigh(vecConf[0].id, costVec)
                    vecConf[0].addConfNeigh(vertex, costVec)
                    print(vertex)
                    print(costVec)
            elif (vecConf[1] != 0):
                print("distinta frecuencia")
                costVec = 0
                changeVert.addConfNeigh(vecConf[0].id, costVec)
                vecConf[0].addConfNeigh(vertex, costVec)
                print(vertex)
                print(costVec) 

    # Asigna una frecuencia al vertice vertex
    def assignVertFreq(self, vertex, idxfreq):
        self.countUsedFreq(idxfreq)
        self.getVtx(vertex).assignFreq(self.frequency[idxfreq], idxfreq)
    
    # Actualiza el contador de uso de la frecuencia freq // Ademas suma el costo de uso costTotalSol    
    def countUsedFreq(self, idxfreq):
        if self.countUsedFrequency[idxfreq] == 0:
            self.addCostUsedFreqToTotal(idxfreq)
            #self.addCostoHastaAhoraFreq(idxfreq)
        self.incrementCounterUsedFreq(idxfreq)
    
    # estas funcion esta bien     
    def addCostUsedFreqToTotal(self, idxfreq):
        self.costTotalSol += self.frequency[idxfreq]
    
    # estas funcion esta bien     
    def addConflFreqToTotal(self, conflFreq):
        self.costTotalSol += conflFreq
    # estas funcion esta bien
    def addCostoHastaAhoraFreq(self, idxfreq):
        self.costHastaAhoraUsedNewFreqTotal += self.frequency[idxfreq]
    
    # estas funcion esta bien     
    def incrementCounterUsedFreq(self, idxfreq):
        self.countUsedFrequency[idxfreq] += 1

    # estas funcion esta bien     
    def decrementCounterUsedFreq(self, idxfreq):
        self.countUsedFrequency[idxfreq] -= 1

    def getNextFreqBarata(self, freqVecinos):
        nextFreq = None
        minc = 100000*100000
        idx = None
        for fqidx in range(0, len(self.frequency)):
            freq = self.frequency[fqidx]
            if fqidx not in freqVecinos:
                cost = 0 if self.countUsedFrequency[fqidx] >= 1 else self.frequency[fqidx]
                if cost < minc:
                    nextFreq = freq
                    minc = cost
                    idx = fqidx
                    break
        return (nextFreq, idx)

    # devuelve si la frecuencia ya fue usada
    def isUsedFreq(self, freqidx):
        return self.countUsedFrequency[freqidx] >= 1
    
    # devuelve la frecuencia mas barata
    def getFreqBarata(self):
        return 0   

    # devuelve el costo de la frecuencia
    def costFreq(self, freqidx):
        return 0 if self.countUsedFrequency[freqidx] >= 1 else self.frequency[freqidx]
        
    # calculando Cota de uso de Freq        
    def addConflitcTotal(self, conflict): 
        self.cotaCostUsedFreq += conflict
    
    def calcularCostChangeFreq(self, vertexidx, freqidx):
        vertex = self.getVtx(vertexidx)
        freqidxVer = vertex.indexfreq
        # Dist
        costNewFreq = self.costFreq(freqidx)
        confNewFreq = vertex.getSumNewConflictNeighs(freqidx) 
        # Old
        confOldFreq = vertex.getSumConflictNeighs() 
        costUsedOldFreq = 0 if (self.countUsedFrequency[freqidxVer]-1) >= 1 else self.frequency[freqidxVer]

        print("costo Total")
        print(self.costTotalSol)
        
        print("cost New Freq")
        print(costNewFreq)
        
        print("conf Old Freq")
        print(confOldFreq)
        
        print("conf New Freq")
        print(confNewFreq)
        
        print("cost Used Old Freq")
        print(costUsedOldFreq)
        
        costNewTotal = self.costTotalSol + costNewFreq - confOldFreq + confNewFreq - costUsedOldFreq
        print("Costo New Total")
        print(costNewTotal)

        return costNewTotal

    def changeVertFreq(self, vertexidx, freqidx):
        vertex = self.getVtx(vertexidx)
        freqidxVer = vertex.indexfreq
        self.decrementCounterUsedFreq(freqidxVer) # actualiza el contador de anterior frecuencia
        # Dist
        costNewFreq = self.costFreq(freqidx)
        confNewFreq = vertex.getSumNewConflictNeighs(freqidx) 
        # Old
        confOldFreq = vertex.getSumConflictNeighs() 
        costUsedOldFreq = self.costFreq(freqidxVer)
        
        self.changeFreqToVertex(vertexidx, freqidx) # actualiza el contador de la nueva frecuencia  suma el nuevo costo : costNewFreq
        self.incrementCounterUsedFreq(freqidx)

        costNewTotal = costNewFreq + confNewFreq - confOldFreq - costUsedOldFreq
        self.addConflFreqToTotal(costNewTotal)

               
    def printToReprep(self):
        resgraph = ""
        resgraph += str(len(self.getVtxs())) + ":"
        for n in self.getVtxs():
            for vec in self.getVtx(n).getNeighs():
                resgraph += str(n) + "-" +  str(vec) + ","
        print(resgraph)

    def __str__(self):
        res = ""
        for n in self.getVtxs():
            res += "vertice \n" + str(self.getVtx(n)) + " \n"        
        return res

def read_line():
    ''' Este es un metodo que lee cada linea del archivo que toma como entrada y
        el mimo sabe terminar cuando no hay mas lineas que recorrer '''
    line = next(sys.stdin).strip()
    while(len(line) == 0):
        line = next(sys.stdin).strip()
    return line
    
def NofSol(sol):
    None

def funcObjetivo(sol):
    None

def tabu_search(solFactibles, sol):
    '''Dado s pertenece a S
        mientras s no es optimo local:
            reemplazar s con s' pertenece a N(s) tq f(s') < f(s)
       retornar s. 
    '''

########### BUSQUEDA LOCAL - PREPARANDO TODO PARA META HEURISTICA TABU SEARCH #################

def searchBestNeighbour(sol):
    print("buscar neighbour")
    solParcial = copy.deepcopy(sol)
    print(range(1, len(sol.getVtxs())+1))
    print(range(0, len(sol.frequency)))

    for ant in range(1, len(sol.getVtxs())+1):
        for freq in range(0, len(sol.frequency)):
            if not (solParcial.getVtx(ant).indexfreq == freq):
                costoSolParc = solParcial.calcularCostChangeFreq(ant, freq)
                print("costo sol_parc")
                print(costoSolParc)
                print("cost de Sol ")
                print(sol.getCostTotalSol())
                if(costoSolParc < sol.getCostTotalSol()):
                    print("es menor")
                    solParcial.changeVertFreq(ant, freq)
                    return solParcial
    return solParcial
                    
def busquedaLocal(solGlobal):
    print("busqueda local")
    solGlobal = searchBestNeighbour(solGlobal) # despues vere que hago con i 
    return solGlobal

##################################################################

########### HEURISTICA CONSTRUCTIVA GOLOSA #################
def calcularAssignVecinosConFreq(vecCurrentIdx, neighbours, newgraph):
    indexFreq = None
    costAssign = None
    vecFreqVecinos = []
    FreqToVecinos = {}

    # separa los vecinos por frecuencia
    for neigh in neighbours:
        vec = neigh[0]
        freqIndexVec = vec.indexfreq
        vecFreqVecinos.append(freqIndexVec)
        if freqIndexVec not in FreqToVecinos.keys():
            FreqToVecinos[freqIndexVec] = [vec]
        else:
            FreqToVecinos[freqIndexVec].append(vec)
   
    # compara y se fija cual es lo que mas conviene considerando todos los conflictos si coinciden la frecuencia 
    for freqIndx in FreqToVecinos.keys():
        conflVec = 0
        for vect in FreqToVecinos[freqIndx]:
            conflVec += graphGuie.getVtx(vecCurrentIdx).getConflictNeigh(vect.id)
        
        if (indexFreq == None and costAssign == None):
            indexFreq = freqIndx
            costAssign = conflVec
        elif conflVec < costAssign:
            indexFreq = freqIndx
            costAssign = conflVec
    
    # compara con la proxima frecuencia a ver con cual se queda
    newFreq = newgraph.getNextFreqBarata(vecFreqVecinos)
    if(newFreq[0] != None):
        freqNewVec = newFreq[0]
        costNewVec = newgraph.costFreq(newFreq[1]) # paso indice indice
        if(costNewVec < costAssign):
            indexFreq = newFreq[1]
            costAssign = costNewVec
     
    # asigna la frecuencia al vertice actual        
    newgraph.assignVertFreq(vecCurrentIdx, indexFreq) # asignamos la frecuencia
    
    for neigh in neighbours: # actualizamos los conflictos si los hay.
        vec = neigh[0]
        idxVecFreq = vec.indexfreq
        if(idxVecFreq == indexFreq):
            costVec = graphGuie.getVtx(vecCurrentIdx).getConflictNeigh(vec.id)
            newgraph.addConflFreqToTotal(costVec)
            newgraph.getVtx(vecCurrentIdx).addConfNeigh(vec.id, costVec)
            newgraph.getVtx(vec.id).addConfNeigh(vecCurrentIdx, costVec)

def vecinosSinConFreq(vecinos):
    vecSinfreq = []
    vecConfreq = []
    for vec in vecinos:
        if (vecinos[vec][0].isFreqNone()):
            vecSinfreq.append(vecinos[vec])
        else:
            vecConfreq.append(vecinos[vec])
    return (vecSinfreq, vecConfreq)

def heuristicaConstructivaGolosa(newgraph):
    for vertex in newgraph.getVtxs():
        vecSin, vecCon = vecinosSinConFreq(newgraph.getNeighbours(vertex))
        if len(vecCon) == 0:
            newgraph.assignVertFreq(vertex, newgraph.getFreqBarata())
        else:
            calcularAssignVecinosConFreq(vertex, vecCon, newgraph)
    return newgraph

###########################################################

def exec_main():
    
    global graphGuie 
    global criterioParadaBusqLocal 
    global criterioParadaTabuSearch 
    global tamListaTabu 

    criterioParadaBusqLocal = 1000
    criterioParadaTabuSearch = 1000
    tamListaTabu = 40

    antFreqConf = read_line().split(" ")
    n = int(antFreqConf[0])  # numero de antenas
    t = int(antFreqConf[1])  # cantidad de frecuencias para usar
    E = int(antFreqConf[2])  # aristas con conflictos si se asignan la misma frecuencia

    frequency = []
    contUsedFrequency = []

    # no estan ordenadas (Tambien quiero reducir la cantidad de frecuencias usadas)
    for freq in range(0, t):
        costFreq = int(read_line().split(" ")[0])
        frequency.append(costFreq)
    
    frequency.sort()
    
    for freq in range(0, t):
        contUsedFrequency.append(0)
        
    graphGuie = Graph(frequency, contUsedFrequency)
    graphSinAssignFreq = Graph(frequency, contUsedFrequency)
    
    for nod in range(1, n+1):
        graphGuie.addVtx(Node(nod))
        graphSinAssignFreq.addVtx(Node(nod))

    for e in range(0, E):
        edgeConf = read_line().split(" ")
        
        vi   = int(edgeConf[0])
        vj   = int(edgeConf[1])
        conf = int(edgeConf[2])
        graphGuie.addConflitcTotal(conf)

        graphGuie.addEdge(vi, vj, conf)
        graphGuie.addEdge(vj, vi, conf)
        
        graphSinAssignFreq.addEdge(vi, vj, 0)
        graphSinAssignFreq.addEdge(vj, vi, 0)

    #print(frequency)
    #print(contUsedFrequency)
    resultHeuristConstGolosa = heuristicaConstructivaGolosa(graphSinAssignFreq)    
    resultBusqLocal = busquedaLocal(resultHeuristConstGolosa)

    print(";;;;;;;;;;;;;;;")
    #print(graphGuie.cotaCostUsedFreq)    
    #print(resultHeuristConstGolosa)
    #print(frequency)
    #print("*************")
    #print(contUsedFrequency)
    #print("//////////////////")
    print(resultBusqLocal.getCostTotalSol())
    #print(resultGraph)

exec_main()
