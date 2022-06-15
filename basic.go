// Basic modules for compact routing. Generates random coordinates and levels, computes clusters

// Cristina Basescu


package main

import (
	"flag"
	"gopkg.in/dedis/onet.v2/log"
	"math"
	"math/rand"
	"strconv"
	"time"
	"fmt"
)

var RndSrc *rand.Rand

type Node struct {
	Name string
	X int
	Y int
	Level int
	ADist []float64
	pDist []string
	cluster map[string]bool
	Bunch map[string]bool
}

// generate random coordinates for nodes
func genNodes(N int, SpaceMax int, K int) []Node{
	nodes := make([]Node, N)
	for i := 0 ; i < N ; i++ {
		nodes[i].X = RndSrc.Intn(SpaceMax)
		//nodes[i].Y = RndSrc.Intn(SpaceMax)
		nodes[i].Y = RndSrc.Intn(SpaceMax)
		nodes[i].Name = strconv.Itoa(i)
		nodes[i].ADist = make([]float64, K)
		for j := 0 ; j < K ; j++ {
			nodes[i].ADist[j] = math.MaxInt64
		}
		nodes[i].pDist = make([]string, K)
		nodes[i].cluster = make(map[string]bool, 0)
		nodes[i].Bunch = make(map[string]bool, 0)
	}



	return nodes
}

func genLevels(N int, K int, nodes []Node) {

	prob := 1.0 / math.Pow(float64(N), 1.0/float64(K))

	for lvl := 0 ; lvl < K ; lvl++ {
		for i := 0; i < N; i++ {
			if nodes[i].Level == lvl - 1 {
				rnd := RndSrc.Float64()
				if rnd < prob {
					nodes[i].Level = lvl
				}
			}
		}
	}

	for i := 0; i < N; i++ {
		log.Lvl1(nodes[i].X, nodes[i].Y, nodes[i].Level)
	}

}

func euclidianDist(a Node, b Node) float64{
	return math.Sqrt(math.Pow(float64(a.X - b.X), 2.0) + math.Pow(float64(a.Y - b.Y), 2.0))
}

func computeADist(nodes []Node, K int) {
	N := len(nodes)
	for i := 0; i < N; i++ {
		crtNode := &nodes[i]
		for lvl := K-1 ; lvl >= 1 ; lvl-- {
			for j := 0; j < N; j++ {
				dist := euclidianDist(nodes[i], nodes[j])
				if nodes[j].Level >= lvl && dist < crtNode.ADist[lvl] {
					crtNode.ADist[lvl] = dist
					crtNode.pDist[lvl] = nodes[j].Name
				}
			}
		}

		crtNode.ADist[0] = 0
		crtNode.pDist[0] = crtNode.Name

		for lvl := K-1 ; lvl >= 1 ; lvl-- {
			if crtNode.ADist[lvl] == crtNode.ADist[lvl-1] {
				crtNode.pDist[lvl-1] = crtNode.pDist[lvl]
			}
		}
	}
}

// This code assumes that each node can maintain an unlimited nr of links to other nodes:
// we basically start with a full mesh as a graph. Because we start with a full mesh, when building a node's cluster with
// Dijkstra, the node ends up maintaining direct links to all nodes in its cluster. Note that, in LEOs, each node has at most 5 direct links,
// thus we probably need something smarter here.

// Each node also maintains direct links to the nodes in its bunch. We most probably cannot get around this without relaxing the resulting
// latency guarantees, so we need to see what impact it has on LEOs.

func computeCluster(nodeIdx int, nodes []Node, K int) {
	crtNode := &nodes[nodeIdx]

	for i := 0 ; i < len(nodes) ; i++ {
		if i != nodeIdx {
			targetLvl := crtNode.Level + 1

			if targetLvl > K {
				log.Error("lvl too large")
			}

			dist := euclidianDist(*crtNode, nodes[i])
			if targetLvl == K || dist < crtNode.ADist[targetLvl] {
				crtNode.cluster[nodes[i].Name] = true
				nodes[i].Bunch[crtNode.Name] = true
			}
		}
	}

	// add everyone who's closer to me than to anyone else at the next level

}


func alg(K int, nodes []Node) {
	for i := K-1 ; i >=0 ; i-- {
		for j := 0 ; j < len(nodes) ; j++ {
			if nodes[j].Level >= i {
				computeCluster(j, nodes, K)
			}
		}
	}
}

// This function only computes the distance between nodeU and nodeV but not the actual path.
// The path should be easy to extract though.
func compactDistance (nodeU Node, nodeV Node, K int, nodes []Node) float64 {
	// same node check
	if nodeU.Name == nodeV.Name {
		return 0.0
	}

	w := nodeU
	i := 0

	for {
		found := false

		//log.LLvl1("w=", w.Name, "i=", i)

		if nodeV.Bunch[w.Name] {
			found = true
			break
		}

		i++
		if found {
			break
		}

		if i == K {
			log.Error("shouldn't get here!", nodeU.Name, nodeV.Name, nodeU.Bunch, nodeV.Bunch)
		}
		aux := nodeU
		nodeU = nodeV
		nodeV = aux

		for _, n := range nodes {
			if n.Name == nodeU.pDist[i] {
				w = n
				break
			}
		}
	}

	return euclidianDist(w,nodeU) + euclidianDist(w,nodeV)
}


func main() {

	K := flag.Int("K", 3, "Number of levels.")
	N := flag.Int("N", 6, "Number of validators.")
	SpaceMax := flag.Int("SpaceMax", 150, "Coordinate space size.")

	flag.Parse()
	RndSrc = rand.New(rand.NewSource(time.Now().UnixNano()))
	nodes := genNodes(*N, *SpaceMax, *K)
	genLevels(*N, *K, nodes)

	computeADist(nodes, *K)

	alg(*K, nodes)

	for i := 0 ; i < *N; i++ {
		for j := i+1 ; j < *N ; j++ {
			fmt.Println("triangle_graph.add_edge(", i, ",", j, ", weight=", euclidianDist(nodes[i], nodes[j]), ")")
		}
	}

	lvl0 := make([]int,0)
	lvl1 := make([]int,0)
	lvl2 := make([]int,0)
	for i := 0 ; i < *N ; i++ {
		if nodes[i].Level ==  0  {
			lvl0 = append(lvl0,i)
		}
		if nodes[i].Level ==  1  {
			lvl1 = append(lvl1,i)
		}
		if nodes[i].Level ==  2  {
			lvl2 = append(lvl2,i)
		}
	}

	fmt.Println("a_list=[", lvl0, ",", lvl1, ",", lvl2, ",[]]")


	for i := 0 ; i < *N ; i++ {
		//log.Lvl1(nodes[i].X, nodes[i].Y, nodes[i].Level, nodes[i].ADist, nodes[i].pDist)
		log.Lvl1(nodes[i].Name, "<",nodes[i].X,",", nodes[i].Y, ">", "lvl=", nodes[i].Level, "cluster = ", nodes[i].cluster, "bunch = ", nodes[i].Bunch)
	}
/*
	for i := 0 ; i < *N ; i++ {
		for j := 0 ; j < *N ; j++ {
			//log.Lvl1(nodes[i].X, nodes[i].Y, nodes[i].Level, nodes[i].ADist, nodes[i].pDist)
			sanityCheckPasses := false
			if euclidianDist(nodes[i], nodes[j]) * (2.0 * float64(*K) - 1.0) >= compactDistance(nodes[i], nodes[j], *K, nodes) {
				sanityCheckPasses = true
			}



			log.Lvl1(nodes[i].Name, nodes[j].Name, "direct distance", euclidianDist(nodes[i], nodes[j]), "compact distance", compactDistance(nodes[i], nodes[j], *K, nodes), "sanity check passes", sanityCheckPasses)

			if sanityCheckPasses == false {
				// shouldn't happen!
				break
			}
		}
	}
*/
}


