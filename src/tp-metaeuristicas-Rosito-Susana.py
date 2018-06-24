import sys
from math import sqrt
import collections

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

    def getNeigh(self, neigh):
        return self.conflictsToNeigh[neigh]
    
    def addConfNeigh(self, neigh, conf):
        self.getNeigh(neigh)[1] = conf 

    def getConflictNeigh(self, neigh):
        return self.getNeigh(neigh)[1]

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

    #------- Abajo - Metodos relacionados al problema --------
    
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

    def getNextFreqBarata(self, freqVecinos):
        nextFreq = None
        min = 100000*100000
        idx = None
        for fqidx in range(0, len(self.frequency)):
            freq = self.frequency[fqidx]
            if fqidx not in freqVecinos:
                cost = 0 if self.countUsedFrequency[fqidx] >= 1 else self.frequency[fqidx]
                if cost < min:
                    nextFreq = freq
                    min = cost
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
def searchBestNeighbour(sol, iter):
    solParcial = sol
    for (a=1 ; a < N+1 ; a++) # por cada antena
        for (f = 0 ; f < T ; f++) # por cada frecuencia
            if (solParcial.get_freq(a) == f) continue
            solParcial.set_freq(a,f)
            if (solParcial.cost() < s.cost())
                return solParcial


def busquedaLocal(solGlobal, cantAntenas, cantFreq):
    iteration = 0 
    while ( no_parada ) # tengo que definir criterio de parada : Cant iteracciones
        solParcial = buscar_vecino_mejor(solGlobal, i, cantAntenas, cantFreq) # despues vere que hago con i 
        solGlobal = solParcial 
        i++
    return solGlobal

##################################################################

########### HEURISTICA CONSTRUCTIVA GOLOSA #################
def calcularAssignVecinosConFreq(vecCurrentIdx, neighbours, graphguie, newgraph):
    freqAssign = None
    costAssign = None
    indexFreq = None
    vecFreqVecinos = []
    for neigh in neighbours:
        vec = neigh[0]
        freqVec = vec.freq
        costVec = graphguie.getVtx(vecCurrentIdx).getConflictNeigh(vec.id)
        vecFreqVecinos.append(newgraph.getVtx(vec.id).indexfreq)
        if (freqAssign == None and costAssign == None):
            freqAssign = freqVec
            costAssign = costVec
            indexFreq = newgraph.getVtx(vec.id).indexfreq

        elif costVec < costAssign:
            freqAssign = freqVec
            costAssign = costVec
            indexFreq = newgraph.getVtx(vec.id).indexfreq
    
    print(vecFreqVecinos)

    newFreq = newgraph.getNextFreqBarata(vecFreqVecinos)
    if(newFreq[0] != None):
        freqNewVec = newFreq[0]
        costNewVec = newgraph.costFreq(newFreq[1]) # paso indice indice
        if(costNewVec < costAssign):
            freqAssign = freqNewVec
            costAssign = costNewVec
            indexFreq = newFreq[1]

    #print("freqAssign") 
    #print(freqAssign)       
    #print("costoAssig")
    #print(costAssign)        
    #print("index de frecuencia con el cual me quedo")
    #print(indexFreq)

    newgraph.assignVertFreq(vecCurrentIdx, indexFreq) # asignamos la frecuencia
    
    for neigh in neighbours: # actualizamos los conflictos si los hay.
        vec = neigh[0]
        idxVecFreq = vec.indexfreq
        #print("index vec freq")
        #print(idxVecFreq)
        if(idxVecFreq == indexFreq):
            #print("son iguales")
            costVec = graphguie.getVtx(vecCurrentIdx).getConflictNeigh(vec.id)
            #print("costo que se suma")
            #print(costVec)
            newgraph.addConflFreqToTotal(costVec)

            newgraph.getVtx(vecCurrentIdx).addConfNeigh(vec.id, costVec)
            newgraph.getVtx(vec.id).addConfNeigh(vecCurrentIdx, costVec)
    #print("costo total hasta ahora")        
    #print(newgraph.costTotalSol)        

def vecinosSinConFreq(vecinos):
    vecSinfreq = []
    vecConfreq = []
    for vec in vecinos:
        if (vecinos[vec][0].isFreqNone()):
            vecSinfreq.append(vecinos[vec])
        else:
            vecConfreq.append(vecinos[vec])
    return (vecSinfreq, vecConfreq)

def heuristicaConstructivaGolosa(graphguie, newgraph):
    for vertex in newgraph.getVtxs():
        print(newgraph.countUsedFrequency)
        vecSin, vecCon = vecinosSinConFreq(newgraph.getNeighbours(vertex))
        if len(vecCon) == 0:
            newgraph.assignVertFreq(vertex, newgraph.getFreqBarata())
        else:
            calcularAssignVecinosConFreq(vertex, vecCon, graphguie, newgraph)
    return newgraph

###########################################################

def exec_main():

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
    resultGraph = heuristicaConstructivaGolosa(graphGuie, graphSinAssignFreq)    

    print(";;;;;;;;;;;;;;;")
    #print(graphGuie.cotaCostUsedFreq)    
    #print(resultGraph)
    #print(frequency)
    #print("*************")
    #print(contUsedFrequency)
    #print("//////////////////")
    #print(resultGraph.costTotalSol)
    print(resultGraph)

exec_main()
