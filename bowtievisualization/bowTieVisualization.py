
#%%
"""
    bowtievisualization is an OpenSource python package for the analysis of
    networkx networks under the bow-tie structure analysis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    Open Points:
    - TODO test graph without SCC
    - TODO setup.py + dependencies
    - TODO
    - TODO
    - TODO
    - TODO

"""
import matplotlib.pyplot as plt
plt.style.use('ggplot')

import matplotlib.mlab as mlab
import numpy as np
import networkx as nx
import math

import time
from datetime import datetime

def _printDict(dic):
    for keys,values in dic.items():
        print(keys , ": ", values)

        
def _calcPYgiven2PointsAndPx(Ax,Ay,Bx,By,Px):
    Py = ((Px - Ax) * ((By-Ay)/(Bx-Ax))) + Ay
    return Py


def __plotLabel(x, y, text, alineation= "left", fontSize = 12):
    y = y - 0.025  # shift y-value for label so that it's below the artist
    plt.text(x, y, text, ha=alineation, family='sans-serif', size=fontSize, bbox={'facecolor':'lightgray', 'alpha':0.5, 'pad':2})

class createBowTieNetworkValues:
    """"""
    #TODO


class BowTieNetworkValues:
    """
    Base class for analyzing a directed graphs' bowtie structure values.

    A BowTieNetworkValues stores the number of nodes of each bowtie component 
    and the percent of each component from the whole graph.

    Parameters
    ----------
    nrNodesAllGraph : input integer (default: 0)
        Total number of nodes in the graph
    nrNodesWeaklyLCC : 
        TODO
    nrNodesOCC = 0
    nrNodesIn = 0
    nrNodesSCC = 0
    nrNodesOut = 0
    nrNodesTubes = 0 
    nrNodesTendrilsIn = 0 
    nrNodesTendrilsOut = 0
    connectedComponentsSizes = []

    """
    def __init__(self,  nrNodesTubes = 0,nrNodesTendrilsIn = 0, nrNodesIn = 0, nrNodesSCC = 0, nrNodesOut = 0, nrNodesTendrilsOut = 0, nrNodesOCC = 0, nrNodesUnidentified = 0, connectedComponentsSizes = []):
        self.nrNodesAllGraph = nrNodesTubes + nrNodesTendrilsIn + nrNodesIn + nrNodesSCC + nrNodesOut + nrNodesTendrilsOut + nrNodesOCC 
        self.nrNodesWeaklyLCC = nrNodesTubes + nrNodesTendrilsIn + nrNodesIn + nrNodesSCC + nrNodesOut + nrNodesTendrilsOut 
        self.nrNodesTubes = nrNodesTubes
        self.nrNodesTendrilsIn = nrNodesTendrilsIn
        self.nrNodesIn = nrNodesIn
        self.nrNodesSCC = nrNodesSCC
        self.nrNodesOut = nrNodesOut
        self.nrNodesTendrilsOut = nrNodesTendrilsOut
        self.nrNodesOCC = nrNodesOCC
        self.nrNodesUnidentified = nrNodesUnidentified
        self.pctWeaklyLCCNodes = 0 if (self.nrNodesAllGraph is None or self.nrNodesAllGraph == 0) else round((self.nrNodesWeaklyLCC / self.nrNodesAllGraph * 100),2)
        self.pctTubeNodes = 0 if (self.nrNodesAllGraph is None or self.nrNodesAllGraph == 0) else round((nrNodesTubes / self.nrNodesAllGraph * 100),2)
        self.pctTendrilsInNodes = 0 if (self.nrNodesAllGraph is None or self.nrNodesAllGraph == 0) else round((nrNodesTendrilsIn / self.nrNodesAllGraph * 100),2)
        self.pctInNodes = 0 if (self.nrNodesAllGraph is None or self.nrNodesAllGraph == 0) else round((nrNodesIn / self.nrNodesAllGraph * 100),2)
        self.pctSCCNodes = 0 if (self.nrNodesAllGraph is None or self.nrNodesAllGraph == 0) else round((nrNodesSCC / self.nrNodesAllGraph * 100),2)
        self.pctOutNodes = 0 if (self.nrNodesAllGraph is None or self.nrNodesAllGraph == 0) else round((nrNodesOut / self.nrNodesAllGraph * 100),2)
        self.pctTendrilsOutNodes = 0 if (self.nrNodesAllGraph is None or self.nrNodesAllGraph == 0) else round((nrNodesTendrilsOut / self.nrNodesAllGraph * 100),2)
        self.pctOCCNodes = 0 if (self.nrNodesAllGraph is None or self.nrNodesAllGraph == 0) else round((nrNodesOCC / self.nrNodesAllGraph * 100),2)
        self.connectedComponentsSizes = [] if (connectedComponentsSizes is None or len(connectedComponentsSizes)==0)else connectedComponentsSizes

    def __repr__(self):
        """
        Prints the BowTieNetworkValues keys and values as a list
        """
        return "<BowTieNetworkValues: %s>" % _printDict(self.__dict__)

class BowTieVisualizationValues:
    """
    Base class for that contains values for the bow-tie shape visualization.
    """
    def __init__(self, hIn = 0, bIn = 0, rSCC = 0, hOut = 0, bOut = 0, hOCC = 0, wOCC = 0,
        bTIn = 0, hTIn = 0, bTOut = 0, hTOut = 0, bTubes = 0,  hTubes = 0, areaIn = 0,
        areaOut = 0, areaOCC = 0, areaSCC = 0, areaTendrilsIn = 0, areaTendrilsOut = 0, areaTubes = 0):
        self.hTubes = hTubes
        self.bTubes = bTubes
        self.hTIn = hTIn
        self.bTIn = bTIn
        self.hIn = hIn
        self.bIn = bIn
        self.rSCC = rSCC
        self.hOut = hOut
        self.bOut = bOut
        self.hTOut = hTOut
        self.bTOut = bTOut
        self.hOCC = hOCC
        self.wOCC = wOCC
        self.areaTubes = areaTubes
        self.areaTendrilsIn = areaTendrilsIn
        self.areaIn = areaIn
        self.areaSCC = areaSCC
        self.areaOut = areaOut
        self.areaTendrilsOut = areaTendrilsOut
        self.areaOCC = areaOCC

    def __repr__(self):
        """
        Prints the BowTieVisualizationValues keys and values as a list
        """
        return "<BowTieVisualizationValues: %s>" % _printDict(self.__dict__)
        occupiedArea = self.areaTubes + self.areaTendrilsIn + self.areaIn + self.areaSCC + self.areaOut + self.areaTendrilsOut + self.areaOCC
        print("relative area Tendrils In: ", self.areaTendrilsIn/occupiedArea)


def getBowTieNetworkValues(gx, debug = False):
    # TODO validate for fully connected network (0 values...)
    """
    From a given Networkx.DiGraph, it calculates values for the Bow Tie based
    on the network topology, like number of nodes in each Bow Tie component
    and the percentage value of each component towards the whole network

    Parameters
    ----------
    gx : Networkx.DiGraph
        path to edgelist file
    separator : str
        character separating the nodes
    weighted : bool
        is a weight given? if ``True`` it is the last element in the edge
        (i.e. ``a,b,2``)
    directed : bool
        are the edges directed or undirected
    header: bool
        if true skip the first row, useful if header row in file

    Returns
    -------
    Network:
        a ``Network`` object obtained from the edgelist
    """
    def __printDebug(*val):
        if debug:
            print(list(val))

    nrNodesAllGraph = 0
    nrNodesWeaklyLCC = 0
    nrNodesOCC = 0
    nrNodesIn = 0
    nrNodesSCC = 0
    nrNodesOut = 0
    nrNodesTubes = 0 
    nrNodesTendrilsIn = 0 
    nrNodesTendrilsOut = 0
    connectedComponentsSizes = []

    # 0. compute nodes in graph
    nrNodesAllGraph = gx.order()

    # 1. compute connected components
    wc = nx.weakly_connected_components(gx)
    weaklyCC = list(wc)

    for i in range(0, len(weaklyCC)):
        connectedComponentsSizes.append(len(weaklyCC[i]))

    connectedComponentsSizes.remove(max(connectedComponentsSizes))

    # 2. keep the largest, the rest are OCCs  
    LCC = max(nx.weakly_connected_components(gx),key=len)
    nrNodesWeaklyLCC = len(LCC)
    LCCn = gx.subgraph(LCC)
    nrNodesOCC = nrNodesAllGraph - nrNodesWeaklyLCC

    # 3. for the largest then compute the strongly CC
    SCC = max(nx.strongly_connected_components(gx.subgraph(LCC)),key=len) 
    SCCn = gx.subgraph(SCC)
    nrNodesSCC = SCCn.order()

    # get diameter for radius of the ego_graph
    #TODO ASK maybe too optimistic
    graphDiameter = nx.diameter(LCCn.to_undirected())*2
    __printDebug("graphDiameter: ", graphDiameter)

    #4. get OUT nodes
    for v in SCCn.nodes():
        __printDebug(v)
        outEgoGraph = nx.ego_graph(LCCn, v, radius=graphDiameter)
        __printDebug( outEgoGraph.nodes())
        break  # done only once because from all nodes in SCC the OUT can be reached

    nodesOut = outEgoGraph.nodes() - SCCn.nodes()
    nrNodesOut = len(nodesOut)

    #5. get IN nodes
    LCCnInverse = LCCn.reverse()
    for v in SCCn.nodes():
        inEgoGraph = nx.ego_graph(LCCnInverse, v, radius=graphDiameter)
        break  # done only once because from all nodes in SCC the OUT can be reached

    nodesIn = inEgoGraph.nodes() - SCCn.nodes()
    nrNodesIn = len(nodesIn)

    # print results
    __printDebug("nodesIn", nodesIn)
    __printDebug("nodesOut", nodesOut)
    __printDebug("SCCn.nodes()", SCCn.nodes())

    #6. get TendrilsIn, TendrilsOut and Tubes
    restNodes = LCCn.nodes() - nodesIn - nodesOut - SCCn.nodes()

    restNodesInitialLenght = len(restNodes)
    nrIt = 0

    nodesTubes = set()
    nodesTendrilsIn = set()
    nodesTendrilsOut = set()

    while (nrNodesTendrilsIn + nrNodesTubes + nrNodesTendrilsOut < restNodesInitialLenght) and nrIt <100:
            nrIt += 1
            __printDebug("restNodes", restNodes)
            for v in restNodes:
                __printDebug("v: ", v)
                isTendrilsIn = False
                isTendrilsOut = False
                isTubes = False

                testEgoGraph = nx.ego_graph(LCCn, v, radius=graphDiameter)
                testEgoGraphInverse = nx.ego_graph(LCCnInverse, v, radius=graphDiameter)

                __printDebug("testEgoGraph: ", testEgoGraph.nodes())
                __printDebug("testEgoGraphInverse: ", testEgoGraphInverse.nodes())

                # get TendrilsIn
                intersectionWithOut = nodesOut.intersection(testEgoGraph.nodes())
                __printDebug("intersectionWithOut: ", intersectionWithOut, len(intersectionWithOut))
                if len(intersectionWithOut) <= 0 :
                    isTendrilsIn = True
                    __printDebug("isTendrilsIn")
                    nodesTendrilsIn.update(testEgoGraphInverse.nodes() - nodesIn)
                    nrNodesTendrilsIn = len(nodesTendrilsIn)
                    break

                if not isTendrilsIn:
                    intersectionInverseWithIN = nodesIn.intersection(testEgoGraphInverse.nodes())
                    __printDebug("intersectionInverseWithIN: ", intersectionInverseWithIN, len(intersectionInverseWithIN))
                    # get Tubes
                    if len(intersectionInverseWithIN) > 0:
                        isTubes = True
                        __printDebug("isTubes")
                        nodesTubes.update(testEgoGraphInverse.nodes() - intersectionInverseWithIN)
                        nrNodesTubes = len(nodesTubes)
                        break
                    # get TendrilsOut
                    else:
                        isTendrilsOut = True
                        __printDebug("isTendrilsOut")
                        nodesTendrilsOut.update(testEgoGraphInverse.nodes() - intersectionInverseWithIN)
                        nrNodesTendrilsOut = len(nodesTendrilsOut)
                        break
            
            # goes out of break to modify the restNodes set
            if isTendrilsIn:
                for dv in nodesTendrilsIn:
                    restNodes.discard(dv)
            
            if isTubes:
                for dv in nodesTubes:
                    restNodes.discard(dv)

            if isTendrilsOut:
                for dv in nodesTendrilsOut:
                    restNodes.discard(dv)

            __printDebug("nrNodesTendrilsIn: ", nrNodesTendrilsIn,  ", nrNodesTubes: ", nrNodesTubes, ", nrNodesTendrilsOut: ", nrNodesTendrilsOut, ", Sum: " )
            __printDebug("nodesTubes: ", nodesTubes)
            __printDebug("nodesTendrilsIn: ", nodesTendrilsIn)
            __printDebug("nodesTendrilsOut: ", nodesTendrilsOut)

    bowTieNetworkValues = BowTieNetworkValues(
        nrNodesOCC = nrNodesOCC,
        nrNodesIn = nrNodesIn,
        nrNodesSCC = nrNodesSCC,
        nrNodesOut = nrNodesOut,
        nrNodesTubes = nrNodesTubes,
        nrNodesTendrilsIn = nrNodesTendrilsIn,
        nrNodesTendrilsOut = nrNodesTendrilsOut,
        connectedComponentsSizes = connectedComponentsSizes
    )

    return bowTieNetworkValues

def getBowTieVisualizationValues(gx):
    bowTieNetworkValues = getBowTieNetworkValues(gx)
    _getBowTieVisualizationValues(bowTieNetworkValues)


def _getBowTieVisualizationValues(bowTieNetworkValues, debug = False):

    def __printDebug(*val):
        if debug:
            print(list(val))

    largestInSide = False
    largestOut = False
    largestSCC = False
    largestTendrilsIn = False
    largestTendrilsOut = False
    equalInOut = False

    bowTieVisualizationValues = BowTieVisualizationValues()

    # horizontal layouts
    # get area for Tubes and max height of Tubes
    if (bowTieNetworkValues.pctTubeNodes <= 5):
        bowTieVisualizationValues.hTubes = 0.05
    else:
        bowTieVisualizationValues.hTubes = bowTieNetworkValues.pctTubeNodes / 100

    __printDebug("hTubes: ", bowTieVisualizationValues.hTubes)


    # get area for OCC and max height of OCC
    if (bowTieNetworkValues.pctOCCNodes<= 5):
        bowTieVisualizationValues.hOCC = 0.05
    else:
        bowTieVisualizationValues.hOCC = bowTieNetworkValues.pctOCCNodes / 100

    __printDebug("hOCC: ", bowTieVisualizationValues.hOCC)

    # largest component in central area (In, SCC, Out, TendrilsIn, TendrilsOut)
    if (bowTieNetworkValues.nrNodesIn == bowTieNetworkValues.nrNodesOut):
        equalInOut = True

    if (bowTieNetworkValues.pctSCCNodes > max(bowTieNetworkValues.pctInNodes, bowTieNetworkValues.pctTendrilsInNodes, bowTieNetworkValues.pctOutNodes, bowTieNetworkValues.pctTendrilsOutNodes)):
        largestSCC = True
    else:
        if (bowTieNetworkValues.pctInNodes + bowTieNetworkValues.pctTendrilsInNodes >= bowTieNetworkValues.pctOutNodes + bowTieNetworkValues.pctTendrilsOutNodes):
            largestInSide = True
            if (bowTieNetworkValues.pctInNodes < bowTieNetworkValues.pctTendrilsInNodes):
                largestTendrilsIn = True
        else:
            largestOut = True
            if (bowTieNetworkValues.pctOutNodes < bowTieNetworkValues.pctTendrilsOutNodes):
                largestTendrilsOut = True

    __printDebug("equalInOut: ", equalInOut)
    __printDebug("largestSCC: ", largestSCC)
    __printDebug("largestInSide: ", largestInSide)
    __printDebug("largestTendrilsIn: ", largestTendrilsIn)
    __printDebug("largestOut: ", largestOut)
    __printDebug("largestTendrilsOut: ", largestTendrilsOut)


    # calculate lengths x-axe of central components
    areaOccupied = bowTieNetworkValues.pctTendrilsInNodes + bowTieNetworkValues.pctInNodes + bowTieNetworkValues.pctSCCNodes + bowTieNetworkValues.pctOutNodes +  bowTieNetworkValues.pctTendrilsOutNodes
    pctAreaTIn = bowTieNetworkValues.pctTendrilsInNodes/areaOccupied
    pctAreaIn = bowTieNetworkValues.pctInNodes/areaOccupied
    pctAreaOut = bowTieNetworkValues.pctOutNodes/areaOccupied
    pctTOut = bowTieNetworkValues.pctTendrilsOutNodes/areaOccupied
    pctAreaSCC = bowTieNetworkValues.pctSCCNodes/areaOccupied

    bowTieVisualizationValues.rSCC = ( pctAreaSCC * 6/8 ) / 2
    __printDebug("pctAreaSCC: ", pctAreaSCC)
    __printDebug("rSCC: ", bowTieVisualizationValues.rSCC)

    bowTieVisualizationValues.bTIn = (1 - 0.05 -bowTieVisualizationValues.rSCC) * pctAreaTIn
    bowTieVisualizationValues.hIn = (1 - 0.05 - bowTieVisualizationValues.rSCC) * pctAreaIn
    bowTieVisualizationValues.hOut = (1 - 0.05 - bowTieVisualizationValues.rSCC) * pctAreaOut
    bowTieVisualizationValues.bTOut = (1 - 0.05 - bowTieVisualizationValues.rSCC) * pctTOut

    __printDebug("bTIn: ", bowTieVisualizationValues.bTIn)
    __printDebug("hIn: ", bowTieVisualizationValues.hIn)
    __printDebug("hOut: ", bowTieVisualizationValues.hOut)
    __printDebug("bTOut: ", bowTieVisualizationValues.bTOut)

    __printDebug("horizontal middle component: ",  bowTieVisualizationValues.bTIn +  bowTieVisualizationValues.hIn + bowTieVisualizationValues.rSCC*2 +  bowTieVisualizationValues.hOut +  bowTieVisualizationValues.bTOut)

    bowTieVisualizationValues.areaSCC = math.pi * pow(bowTieVisualizationValues.rSCC, 2)

    __printDebug("areaSCC: ", bowTieVisualizationValues.areaSCC)

    maxPctArea = max(bowTieNetworkValues.pctTendrilsInNodes, bowTieNetworkValues.pctInNodes, bowTieNetworkValues.pctSCCNodes , bowTieNetworkValues.pctOutNodes, bowTieNetworkValues.pctTendrilsOutNodes)
    __printDebug("maxPctArea: ", maxPctArea)


    #all height are now calculated based on the are of SCC
    if bowTieVisualizationValues.areaSCC > 0:
    # largestTendrilsIn - Tendrils In have the largest height in the central area
        if bowTieNetworkValues.nrNodesTendrilsIn > 0:
            bowTieVisualizationValues.areaTendrilsIn = bowTieNetworkValues.pctTendrilsInNodes * (bowTieVisualizationValues.areaSCC/bowTieNetworkValues.pctSCCNodes)
            bowTieVisualizationValues.hTIn = (bowTieVisualizationValues.areaTendrilsIn/bowTieVisualizationValues.bTIn)
            __printDebug("calc areaTendrilsIn: ", bowTieVisualizationValues.bTIn*bowTieVisualizationValues.hTIn)
        # calc In
        if (bowTieNetworkValues.nrNodesIn > 0):
            bowTieVisualizationValues.areaIn = bowTieNetworkValues.pctInNodes * (bowTieVisualizationValues.areaSCC/bowTieNetworkValues.pctSCCNodes)
            bowTieVisualizationValues.bIn = ((2 * bowTieVisualizationValues.areaIn) / bowTieVisualizationValues.hIn) 
        # calc out
        if(bowTieNetworkValues.nrNodesOut>0):
            bowTieVisualizationValues.areaOut = bowTieNetworkValues.pctOutNodes * (bowTieVisualizationValues.areaSCC/bowTieNetworkValues.pctSCCNodes)
            bowTieVisualizationValues.bOut = ((2 * bowTieVisualizationValues.areaOut) / bowTieVisualizationValues.hOut)
        # calc TendrilOut
        if(bowTieNetworkValues.nrNodesTendrilsOut>0):
            bowTieVisualizationValues.areaTendrilsOut = bowTieNetworkValues.pctTendrilsOutNodes * (bowTieVisualizationValues.areaSCC/bowTieNetworkValues.pctSCCNodes)
            bowTieVisualizationValues.hTOut = (bowTieVisualizationValues.areaTendrilsOut/bowTieVisualizationValues.bTOut)

    else:
        #TODO
        print("Undefined action for graphs without SCC")


    """
    # get heights
    if largestSCC:
        #TODO
        maxPctArea = 0
    else:
        if (largestInSide):
            if (largestTendrilsIn):
                # largestTendrilsIn - Tendrils In have the largest height in the central area
                bowTieVisualizationValues.hTIn = 1 - bowTieVisualizationValues.hOCC - bowTieVisualizationValues.hTubes
                bowTieVisualizationValues.areaTendrilsIn = bowTieVisualizationValues.bTIn * bowTieVisualizationValues.hTIn
                # calc TendrilOut
                if(bowTieNetworkValues.nrNodesTendrilsOut>0):
                    bowTieVisualizationValues.areaTendrilsOut = bowTieNetworkValues.pctTendrilsOutNodes * (bowTieVisualizationValues.areaTendrilsIn/bowTieNetworkValues.pctTendrilsInNodes)
                    bowTieVisualizationValues.hTOut = (bowTieVisualizationValues.areaTendrilsOut/bowTieVisualizationValues.bTOut) /100
                # calc In
                if (bowTieNetworkValues.nrNodesIn > 0):
                    bowTieVisualizationValues.areaIn = bowTieNetworkValues.pctInNodes * (bowTieVisualizationValues.areaTendrilsIn/bowTieNetworkValues.pctTendrilsInNodes)
                    bowTieVisualizationValues.bIn = ((2 * bowTieVisualizationValues.areaIn) / bowTieVisualizationValues.hIn) /100
                # calc out
                if(bowTieNetworkValues.nrNodesOut>0):
                    bowTieVisualizationValues.areaOut = bowTieNetworkValues.pctOutNodes * (bowTieVisualizationValues.areaTendrilsIn/bowTieNetworkValues.pctTendrilsInNodes)
                    bowTieVisualizationValues.bOut = ((2 * bowTieVisualizationValues.areaOut) / bowTieVisualizationValues.hOut)/100
            else:
                bowTieVisualizationValues.bIn = 1 - bowTieVisualizationValues.hOCC - bowTieVisualizationValues.hTubes
                bowTieVisualizationValues.areaIn = bowTieVisualizationValues.bIn * bowTieVisualizationValues.hIn / 2
                # calc out
                if(bowTieNetworkValues.nrNodesOut>0):
                    bowTieVisualizationValues.bOut = (2 * bowTieVisualizationValues.areaIn * bowTieNetworkValues.pctOutNodes) / (bowTieVisualizationValues.hOut * bowTieNetworkValues.pctInNodes)
                    bowTieVisualizationValues.areaOut = bowTieVisualizationValues.bOut * bowTieVisualizationValues.hOut / 2
                # calc tendrils in
                if(bowTieNetworkValues.nrNodesTendrilsIn>0):
                    bowTieVisualizationValues.areaTendrilsIn = bowTieVisualizationValues.areaIn * bowTieNetworkValues.pctTendrilsInNodes / bowTieNetworkValues.pctInNodes
                    bowTieVisualizationValues.hTIn = bowTieVisualizationValues.areaTendrilsIn / bowTieVisualizationValues.bTIn
                # calc tendrils out
                if(bowTieNetworkValues.nrNodesTendrilsOut>0):
                    bowTieVisualizationValues.areaTendrilsOut = bowTieVisualizationValues.areaIn * bowTieNetworkValues.pctTendrilsOutNodes / bowTieNetworkValues.pctInNodes
                    bowTieVisualizationValues.hTOut = bowTieVisualizationValues.areaTendrilsOut / bowTieVisualizationValues.bTOut

        else:
            if (largestTendrilsOut):
                bowTieVisualizationValues.hTOut = 1 - bowTieVisualizationValues.hOCC - bowTieVisualizationValues.hTubes  
                bowTieVisualizationValues.areaTendrilsOut = bowTieVisualizationValues.bTOut * bowTieVisualizationValues.hTOut
                # calc TendrilIn
                if(bowTieNetworkValues.nrNodesTendrilsIn>0):
                    bowTieVisualizationValues.areaTendrilsIn = bowTieVisualizationValues.areaTendrilsOut * bowTieNetworkValues.pctTendrilsInNodes/bowTieNetworkValues.pctTendrilsOutNodes
                    bowTieVisualizationValues.hTIn = bowTieVisualizationValues.areaTendrilsIn/bowTieVisualizationValues.bIn
                # calc In
                if(bowTieNetworkValues.nrNodesIn>0):
                    bowTieVisualizationValues.bIn = (2 * bowTieVisualizationValues.areaTendrilsIn * bowTieNetworkValues.pctInNodes) / (bowTieVisualizationValues.hIn * bowTieNetworkValues.pctTendrilsInNodes)
                    bowTieVisualizationValues.areaIn = bowTieVisualizationValues.bIn * bowTieVisualizationValues.hIn / 2
                # calc out
                bowTieVisualizationValues.bOut = (2 * bowTieVisualizationValues.areaOut * bowTieNetworkValues.pctOutNodes) / (bowTieVisualizationValues.hOut * bowTieNetworkValues.pctOutNodes)
                bowTieVisualizationValues.areaOut = bowTieVisualizationValues.bOut * bowTieVisualizationValues.hOut / 2

            else:
                bowTieVisualizationValues.bOut = 1 - bowTieVisualizationValues.hOCC - bowTieVisualizationValues.hTubes 
                bowTieVisualizationValues.areaOut  = bowTieVisualizationValues.bOut * bowTieVisualizationValues.hOut / 2
                # calc In
                if(bowTieNetworkValues.nrNodesIn>0):
                    bowTieVisualizationValues.bIn = (2 * bowTieVisualizationValues.areaOut * bowTieNetworkValues.pctInNodes) / (bowTieVisualizationValues.hIn * bowTieNetworkValues.pctOutNodes)
                    bowTieVisualizationValues.areaIn = bowTieVisualizationValues.bIn * bowTieVisualizationValues.hIn / 2
                # calc tendrils in
                if(bowTieNetworkValues.nrNodesTendrilsIn>0):
                    bowTieVisualizationValues.areaTendrilsIn = bowTieVisualizationValues.areaIn * bowTieNetworkValues.pctTendrilsInNodes / bowTieNetworkValues.pctInNodes
                    bowTieVisualizationValues.hTIn = bowTieVisualizationValues.areaTendrilsIn / bowTieVisualizationValues.bTIn
                # calc tendrils out
                if(bowTieNetworkValues.nrNodesTendrilsOut>0):
                    bowTieVisualizationValues.areaTendrilsOut = bowTieVisualizationValues.areaOut * bowTieNetworkValues.pctTendrilsOutNodes / bowTieNetworkValues.pctOutNodes
                    bowTieVisualizationValues.hTOut = bowTieVisualizationValues.areaTendrilsOut / bowTieVisualizationValues.bTOut
    """
    # get area of OCC and width (widht is horizontal)
    if(bowTieNetworkValues.nrNodesOCC>0):
        if(bowTieNetworkValues.nrNodesIn>0):
            bowTieVisualizationValues.areaOCC = bowTieVisualizationValues.areaIn * bowTieNetworkValues.pctOCCNodes/bowTieNetworkValues.pctInNodes
        elif(bowTieNetworkValues.nrNodesOut>0):
            bowTieVisualizationValues.areaOCC = bowTieVisualizationValues.areaOut * bowTieNetworkValues.pctOCCNodes/bowTieNetworkValues.pctOutNodes
        else:
            #TODO
            bowTieVisualizationValues.areaOCC = 0
        bowTieVisualizationValues.wOCC = (4 * bowTieVisualizationValues.areaOCC) / (math.pi * bowTieVisualizationValues.hOCC)


    # get area of Tubes
    if(bowTieNetworkValues.nrNodesTubes>0):
        bowTieVisualizationValues.areaTubes = bowTieVisualizationValues.areaIn * bowTieNetworkValues.pctTubeNodes/bowTieNetworkValues.pctInNodes
        bowTieVisualizationValues.bTubes = bowTieVisualizationValues.areaTubes/bowTieVisualizationValues.hTubes

    return bowTieVisualizationValues


def plotPieNetworkValues(bowTieNetworkValues):
    # get basic plot of node distribution
    labels = ('IN','SCC','OUT','Tubes','TendrilsIn','TendrilsOut','OCC')
    sizes = (bowTieNetworkValues.nrNodesIn, bowTieNetworkValues.nrNodesSCC, bowTieNetworkValues.nrNodesOut, bowTieNetworkValues.nrNodesTubes, bowTieNetworkValues.nrNodesTendrilsIn, bowTieNetworkValues.nrNodesTendrilsOut, bowTieNetworkValues.nrNodesOCC)

    fig, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.show()


def plotBowTie(bowTieNetworkValues, sizeFactor = 2, saveSVGfile = False, svgFileSufix = "0", printVisualizationValues = False):
    """
    Plots the bow tie shape visualization based on the network characteristics.
    Parameters
    ----------
    bowTieNetworkValues : object of class BowTieNetworkValues
    TODO
    dictBowTieVisualizationValues : dictionary
        Name of the 
    sizeFactor : integer, default = 2
        Size of the plot related to the standard ploting in the python client
    Returns
    -------
    dictBowTieVisualizationValues : plot
        dictBowTieVisualizationValues plot ??
    """
    from matplotlib.patches import Circle, Polygon, Ellipse, Arrow
    from matplotlib.collections import PatchCollection
    import matplotlib.pyplot as plt

    bowTieVisualizationValues = _getBowTieVisualizationValues(bowTieNetworkValues)

    # Fixing random state for reproducibility
    np.random.seed(19101101)
    fig, ax = plt.subplots()
    #plt.axis("off")
    patches = []

    # size
    params = plt.gcf()
    plotSize = params.get_size_inches()
    params.set_size_inches( (plotSize[0]*sizeFactor, plotSize[1]*sizeFactor) , forward=True)

    # maximal height of component shape 
    maxHeight = max(bowTieVisualizationValues.hTIn, bowTieVisualizationValues.bIn, (bowTieVisualizationValues.rSCC*2), bowTieVisualizationValues.bOut, bowTieVisualizationValues.hTOut)
    YAxe = 1 - bowTieVisualizationValues.hTubes - (maxHeight / 2) - 0.01

    # Rectangle Tendrils In
    RTinx = (1 - (bowTieVisualizationValues.bTIn + bowTieVisualizationValues.hIn + 2*bowTieVisualizationValues.rSCC + bowTieVisualizationValues.hOut + bowTieVisualizationValues.bTOut + 0.1))/2
    if(bowTieNetworkValues.nrNodesTendrilsIn>0):
        coordsTInStr = str(RTinx) + " " + str(YAxe - bowTieVisualizationValues.hTIn/2) + ";" + str(RTinx) + " "+ str(YAxe + bowTieVisualizationValues.hTIn/2) + ";" + str(RTinx + bowTieVisualizationValues.bTIn) + " " + str(YAxe + bowTieVisualizationValues.hTIn/2) + ";" + str(RTinx + bowTieVisualizationValues.bTIn) + " " + str(YAxe - bowTieVisualizationValues.hTIn/2)
        coordsTIn = np.matrix(coordsTInStr)
        polygon = Polygon(coordsTIn, True)
        patches.append(polygon)
        labelTIn = str("Tendrils In\n") + str(bowTieNetworkValues.pctTendrilsInNodes) + str(" %\n n=") + str(bowTieNetworkValues.nrNodesTendrilsIn)
        __plotLabel(RTinx , YAxe - bowTieVisualizationValues.hTIn/2 - 0.02, labelTIn , alineation = "center")


    # Triangle IN
    Ax = RTinx + bowTieVisualizationValues.bTIn + 0.05
    if(bowTieNetworkValues.nrNodesIn>0):
        coordsInStr = str(Ax) + " " + str(YAxe + bowTieVisualizationValues.bIn/2) + ";" + str(Ax) + " "+ str(YAxe - bowTieVisualizationValues.bIn/2) + ";" + str(bowTieVisualizationValues.hIn + Ax) + " " + str(YAxe)
        coordsIn = np.matrix(coordsInStr)
        polygon = Polygon(coordsIn, True)
        patches.append(polygon)
        labelIn = str("IN\n") + str(bowTieNetworkValues.pctInNodes) + str(" %\n n=") + str(bowTieNetworkValues.nrNodesIn)
        __plotLabel(Ax + bowTieVisualizationValues.hIn/2, YAxe ,  labelIn, alineation="center") #TODO optimization label according to shape size...

    # Circle SCC
    Vx = Ax + bowTieVisualizationValues.hIn + bowTieVisualizationValues.rSCC
    if(bowTieNetworkValues.nrNodesSCC>0):
        circle = Circle((Vx, str(YAxe)), bowTieVisualizationValues.rSCC)
        patches.append(circle)
        labelSCC = str("SCC\n") + str(bowTieNetworkValues.pctSCCNodes) + str(" %\n n=") + str(bowTieNetworkValues.nrNodesSCC)
        __plotLabel(Vx, YAxe ,  labelSCC, alineation="center")

    # Triangle OUT
    AOx = Vx + bowTieVisualizationValues.rSCC
    if(bowTieNetworkValues.nrNodesOut>0):
        coordsOutStr = str(AOx) + " " + str(YAxe) + ";" + str(AOx + bowTieVisualizationValues.hOut) +" " + str(YAxe + bowTieVisualizationValues.bOut/2) + ";" + str(AOx + bowTieVisualizationValues.hOut) + " " + str(YAxe - bowTieVisualizationValues.bOut/2)
        coordsOut = np.matrix(coordsOutStr)
        polygon = Polygon(coordsOut, True)
        patches.append(polygon)
        labelOut = str("OUT\n") + str(bowTieNetworkValues.pctOutNodes) + str(" %\n n=") + str(bowTieNetworkValues.nrNodesOut)
        __plotLabel(AOx + bowTieVisualizationValues.hOut/2, YAxe ,  labelOut, alineation="center")

    # Rectangle Tendrils Out
    RTOutx = AOx + bowTieVisualizationValues.hOut + 0.05
    if(bowTieNetworkValues.nrNodesTendrilsOut>0):
        coordsTOutStr = str(RTOutx) + " " + str(YAxe - bowTieVisualizationValues.hTOut/2) + ";" + str(RTOutx) + " "+ str(YAxe + bowTieVisualizationValues.hTOut/2) + ";" + str(RTOutx + bowTieVisualizationValues.bTOut) + " " + str(YAxe + bowTieVisualizationValues.hTOut/2) + ";" + str(RTOutx + bowTieVisualizationValues.bTOut) + " " + str(YAxe - bowTieVisualizationValues.hTOut/2)
        coordsTOut = np.matrix(coordsTOutStr)
        polygon = Polygon(coordsTOut, True)
        patches.append(polygon)
        labelTOut = str("Tendrils Out\n") + str(bowTieNetworkValues.pctTendrilsOutNodes) + str(" %\n n=") + str(bowTieNetworkValues.nrNodesTendrilsOut)
        __plotLabel(RTOutx, YAxe - bowTieVisualizationValues.hTOut/2 - 0.02,  labelTOut)

    # Rectangle Tubes
    RTubesx = Vx
    if(bowTieNetworkValues.nrNodesTubes>0):
        coordsTInStr = str(RTubesx + bowTieVisualizationValues.bTubes/2) + " " + str(1) + ";" + str(RTubesx + bowTieVisualizationValues.bTubes/2) + " "+ str(1 - bowTieVisualizationValues.hTubes) + ";" + str(RTubesx - bowTieVisualizationValues.bTubes/2) + " " + str(1 - bowTieVisualizationValues.hTubes) + ";" + str(RTubesx - bowTieVisualizationValues.bTubes/2) + " " + str(1)
        coordsTIn = np.matrix(coordsTInStr)
        polygon = Polygon(coordsTIn, True)
        patches.append(polygon)
        labelTubes = str("Tubes\n") + str(bowTieNetworkValues.pctTubeNodes) + str(" %\n n=") + str(bowTieNetworkValues.nrNodesTubes)
        __plotLabel(RTubesx, 1 - bowTieVisualizationValues.hTubes/2,  labelTubes, alineation="center")

    #TODO
    labelOCC = str("OCC\n") + str(bowTieNetworkValues.pctOCCNodes) + str(" %\n n=") + str(bowTieNetworkValues.nrNodesOCC)

    if(bowTieNetworkValues.nrNodesOCC>0):
        if (bowTieNetworkValues.pctOCCNodes < 20):
            # elipse
            ellipse = Ellipse((Vx, 0.2), 0.2, 0.1)
            patches.append(ellipse)
            __plotLabel(Vx, 0.2,  labelOCC, alineation="center")
        else:
            RTubesx = Vx
            coordsTInStr = str(RTubesx + bowTieVisualizationValues.hOCC/2) + " " + str(1) + ";" + str(RTubesx + bowTieVisualizationValues.wOCC/2) + " "+ str(1 - bowTieVisualizationValues.hTubes) + ";" + str(RTubesx - bowTieVisualizationValues.bTubes/2) + " " + str(1 - bowTieVisualizationValues.hTubes) + ";" + str(RTubesx - bowTieVisualizationValues.bTubes/2) + " " + str(1)
            coordsTIn = np.matrix(coordsTInStr)
            polygon = Polygon(coordsTIn, True)
            patches.append(polygon)


    # Arrows for Tubes rectangle
    if(bowTieNetworkValues.nrNodesTubes>0):
        Sx = Ax + (bowTieVisualizationValues.hIn/2)
        Sy = _calcPYgiven2PointsAndPx(Ax, YAxe+bowTieVisualizationValues.bIn/2, Ax + bowTieVisualizationValues.hIn, YAxe, Ax + bowTieVisualizationValues.hIn/2)
        Ex = RTubesx - bowTieVisualizationValues.bTubes/2
        Ey = 1 - (bowTieVisualizationValues.hTubes/2)
        arrowIn = Arrow( Sx, Sy , Ex-Sx, Ey-Sy, width=0.03)
        patches.append(arrowIn)

        SAOutx = RTubesx + bowTieVisualizationValues.bTubes/2
        SAOuty = _calcPYgiven2PointsAndPx(AOx, YAxe, AOx + bowTieVisualizationValues.hOut, YAxe + bowTieVisualizationValues.bOut/2, AOx + bowTieVisualizationValues.hOut/2)
        arrowOut = Arrow( SAOutx, Ey ,  AOx + bowTieVisualizationValues.hOut/2- SAOutx, SAOuty - Ey, width=0.03)
        patches.append(arrowOut)


    #TODO colors
    colors = 100*np.random.rand(len(patches))
    p = PatchCollection(patches, alpha=0.4)
    p.set_array(np.array(colors))
    ax.add_collection(p)

    if saveSVGfile:
        imgName = str("bowTie" + str(svgFileSufix) + str(datetime.today().strftime("%d-%m-%y %H %M %S")) + ".svg")
        plt.savefig(imgName)
    plt.show()

    if printVisualizationValues:
        print(bowTieVisualizationValues)





def showBowTieVisualization(gx, plotPieValues = True, printBowTieNetworkValues = True, saveSVGfile = False, plotSizeFactor = 2, svgFileSufix = " ", printVisualizationValues = False):
    """
    Displays the Bow Tie visualization of the given network with other options.
    Parameters
    ----------
    gx : Networkx DiGraph
        Networkx DiGraph
    plotPieValues : boolean, default True
        Name of the 
    saveSVGfile : boolean, default False
        Saves a svg file of the Bow Tie plot visualization on the disc.
    plotSizeFactor = number, default 2
        Based on the default size of the plot it multiplies it size
    svgFileSufix = string, default " "
        Sufix for the svg file of the plot visualization
    Returns
    -------
    ?? : ??
        plot??
    """

    bowTieNetworkValues = getBowTieNetworkValues(gx)
    if plotPieValues:
        plotPieNetworkValues(bowTieNetworkValues)
    if printBowTieNetworkValues:
        print(bowTieNetworkValues)
    plotBowTie(bowTieNetworkValues, sizeFactor=plotSizeFactor, saveSVGfile=saveSVGfile, svgFileSufix=svgFileSufix, printVisualizationValues = printVisualizationValues)



def getEnsambleBowTieNetworkValues(gx, samples = 5, model = "configuration"):
    """
    Returns the averaged Bow Tie Network Values for graphs created based on the
    network topology of the given graph under a certain model.
    Parameters
    ----------
    gx : Networkx DiGraph
        Networkx DiGraph
    samples : number of graph in the ensamble
        Name of the 
    model : "configuration", "ER", ""
        Model with wich the graphs will be created, valid modesl are:
        "configuration" : directed_configuration_model
        "ER" : Erdős-Rényi Gnp graph
        default is "configuration", this model will be used if no other valid name is given.
    Returns
    -------
    R : bowTieNetworkValues
        Averaged Network Values for all the generated graphs
    """
    din = list(d for n, d in gx.in_degree())
    dout = list(d for n, d in gx.out_degree())
    n = gx.order()
    m = nx.number_of_edges(gx)
    p = m/(n*(n-1))

    nrNodesAllGraph = 0
    nrNodesWeaklyLCC = 0
    nrNodesOCC = 0
    nrNodesIn = 0
    nrNodesSCC = 0
    nrNodesOut = 0
    nrNodesTubes = 0
    nrNodesTendrilsIn = 0
    nrNodesTendrilsOut = 0

    if (model == "ER"):
        print("Gnp values. n:",n, ", p:", p)


    # add the values for each realization
    for i in range(samples):
        if (model == "ER"):
            g = nx.fast_gnp_random_graph(n=n, p=p, directed = True)
        else :
            g = nx.directed_configuration_model(din,dout)
        
        btV = getBowTieNetworkValues(g)
        nrNodesAllGraph += btV.nrNodesAllGraph
        nrNodesWeaklyLCC += btV.nrNodesWeaklyLCC
        nrNodesOCC += btV.nrNodesOCC
        nrNodesIn += btV.nrNodesIn
        nrNodesSCC += btV.nrNodesSCC
        nrNodesOut += btV.nrNodesOut
        nrNodesTubes += btV.nrNodesTubes
        nrNodesTendrilsIn += btV.nrNodesTendrilsIn
        nrNodesTendrilsOut += btV.nrNodesTendrilsOut


    # average the values
    nrNodesAllGraph = nrNodesAllGraph/samples
    nrNodesWeaklyLCC = nrNodesWeaklyLCC/samples
    nrNodesOCC = nrNodesOCC/samples
    nrNodesIn = nrNodesIn/samples
    nrNodesSCC = nrNodesSCC/samples
    nrNodesOut = nrNodesOut/samples
    nrNodesTubes = nrNodesTubes/samples
    nrNodesTendrilsIn = nrNodesTendrilsIn/samples
    nrNodesTendrilsOut = nrNodesTendrilsOut/samples


    # create dictionary
    bowTieNetworkValues = BowTieNetworkValues(
        nrNodesOCC = nrNodesOCC,
        nrNodesIn = nrNodesIn,
        nrNodesSCC = nrNodesSCC,
        nrNodesOut = nrNodesOut,
        nrNodesTubes = nrNodesTubes,
        nrNodesTendrilsIn = nrNodesTendrilsIn,
        nrNodesTendrilsOut = nrNodesTendrilsOut
    )

    return bowTieNetworkValues



def showBowTieVisualizationWithEnsamble(gx, samples = 5, model = "configuration", plotPieValues = True, printBowTieNetworkValues = True, saveSVGfile = False, plotSizeFactor = 2, printVisualizationValues=False):
    print("---------------------------- Original graph ----------------------------")
    showBowTieVisualization(gx = gx, plotPieValues=plotPieValues, saveSVGfile=saveSVGfile, plotSizeFactor=plotSizeFactor, svgFileSufix="Original", printVisualizationValues=printVisualizationValues)

    print("---------------------------- Ensamble values graph ----------------------------")
    ensambleBowTieNewtorkValues = getEnsambleBowTieNetworkValues(gx, samples = samples, model = model)

    if plotPieValues:
        plotPieNetworkValues(ensambleBowTieNewtorkValues)
    if printBowTieNetworkValues:
        print(ensambleBowTieNewtorkValues)
    
    plotBowTie(ensambleBowTieNewtorkValues, sizeFactor=plotSizeFactor, saveSVGfile=saveSVGfile, svgFileSufix="Ensamble")
    if printVisualizationValues:
        print(_getBowTieVisualizationValues(ensambleBowTieNewtorkValues))

