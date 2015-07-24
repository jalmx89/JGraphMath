/**
 * JGraphMath
 * Copyright (C) 2015 Jeremiah N. Hankins
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

package com.jnhankins.graph;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * This class's main method is an experimental method for quickly generating
 * the set of all non-isomorphic graphs of a given order.
 * 
 * @author Jeremiah N. Hankins
 */
public class Test {
    public static void main(String[] args) {
        final int minOrder = 1;
        final int maxOrder = 9;
        final int numThreads = 16;

        // For each order of graph...
        for (int order = minOrder; order <= maxOrder; order++) {

            // Keep track of how long it takes to generate all graphs of this 
            // order
            long startTime = System.nanoTime();

            // Create an array which will hold each isomphically unique graph.
            List<Graph> iso = new ArrayList();

            // Create an map from BFS-SLR numbers keys to values which are lists 
            // of graphs which have the key's BFS-SLR number
            Map<long[], List<Graph>> map = new TreeMap(new Comparator<long[]>() {
                @Override
                public int compare(long[] a, long[] b) {
                    if (a.length < b.length) return -1;
                    if (b.length < a.length) return 1;
                    for (int i=a.length-1; i>=0; i--) {
                        if (a[i] < b[i]) return -1;
                        if (b[i] < a[i]) return 1;
                    }
                    return 0;
                }
            });

            // Create an array to contain all the threads that will be employed
            Thread[] threads = new Thread[numThreads];

            // Create and start all the threads...
            for (int t=0; t<threads.length; t++) {
                final int o = order; // order
                final int s = t;     // seed number

                // Create a new thread
                threads[t] =  new Thread() {
                    @Override
                    public void run() {
                        // Each thread must keep track of the isomorphically
                        // disticnt graphs that it encounters locally 
                        List<Graph> localIso = new ArrayList();

                        // Each thread must also keep a map from BFS-SLR numbers
                        // to lists containing isomorphically distinct graphs it
                        // encounters wich have that number
                        Map<long[], List<Graph>> localMap = new TreeMap(new Comparator<long[]>() {
                            @Override
                            public int compare(long[] a, long[] b) {
                                if (a.length < b.length) return -1;
                                if (b.length < a.length) return 1;
                                for (int i=a.length-1; i>=0; i--) {
                                    if (a[i] < b[i]) return -1;
                                    if (b[i] < a[i]) return 1;
                                }
                                return 0;
                            }
                        });

                        // Dont use the thread if the seed is too big
                        if (s >= Graph.getNumGraphs(o))
                            return;

                        // Create a graph using the order and seed
                        Graph g = Graph.makeFromBigInt(o, BigInteger.valueOf(s));

                        // Iterate through all of the graphs of the specified 
                        // order, incrementing the graph by the total number of 
                        // threads. Stop once the highest bit in the graph's 
                        // ajacency matrix has been set, since any graph with 
                        // the highest bit set must be the complement of some 
                        // graph which has previously been encounterd.
                        mainLoop : do {

                            // Get the BFS-SLR number for the current graph
                            long[] gkey = g.getBfsSlrTreeNumber();

                            // Get the list of graphs for the specified number
                            List<Graph> list = localMap.get(gkey);

                            // If no graphs have been encoutned with the number,
                            // then create a new list and store it in the map
                            if (list == null) {
                                list = new ArrayList();
                                localMap.put(gkey, list);
                            }

                            // Ensure that the graph is not isomorphic with any
                            // graph that has the same BFS-SLR number
                            for (int i=0; i<list.size(); i++)
                                if (g.isIsomorphic(list.get(i)))
                                    continue mainLoop;

                            // Assert: The graph is not isomorphic with any 
                            // graph previously encountered by this thread

                            // Create a copy of the graph, and store it in the 
                            // list and the map
                            Graph c = Graph.makeCopy(g);
                            list.add(c);
                            localIso.add(c);

                            // Create a complement of the graph
                            Graph d = Graph.makeCopy(g);
                            d.complement();

                            // If the graph's complement is not isomorphic with 
                            // itself, then it must not be isomorphic with any 
                            // other graph encountered by this thread so far
                            if (!d.isIsomorphic(g)) {
                                // Get the complement's number
                                long[] dkey = d.getBfsSlrTreeNumber();
                                // Get the number's list
                                list = localMap.get(dkey);
                                // Create the list if it doesnt exist
                                if (list == null) {
                                    list = new ArrayList();
                                    localMap.put(dkey, list);
                                }
                                // Store the complement
                                list.add(d);
                                localIso.add(d);
                            }

                            // Increment the graph by the number of threads (so 
                            // that each thread tests a unique set of labled 
                            // graphs), and stop once all graphs have been 
                            // tested
                        } while (g.increment(numThreads) && !g.getEdge(g.getMaximumSize()-1));

                        // Synchronize with the other threads
                        synchronized (iso) {
                            // For every unuque graph that this thread has 
                            // encounterd
                            outerLoop : for (Graph h : localIso) {
                                // Get the graph's number
                                long[] hkey = h.getBfsSlrTreeNumber();
                                // Get the number's list from the global map
                                List<Graph> list = map.get(hkey); 
                                // Create the list if it does't exist
                                if (list == null) {
                                    list = new ArrayList();
                                    map.put(hkey, list);
                                }
                                // Ensure the graph is not in the list 
                                // (i.e. found by annother thread)
                                for (int i=0; i<list.size(); i++)
                                    if (h.isIsomorphic(list.get(i)))
                                        continue outerLoop;
                                // Store the graph globally
                                list.add(h);
                                iso.add(h);
                            }
                        }
                    }
                };
                // Start the thread
                threads[t].start();
            }

            // Rejoin all the threads
            for (int t=0; t<threads.length; t++)
                try { threads[t].join(); } catch (Exception ignore) {}
            
            // Display results
            System.out.println();
            System.out.println("Finished order: "+order+" size: "+iso.size());
            System.out.println("Elapsed Time: "+((System.nanoTime()-startTime)*1e-9)+" sec");
        }
    }
}
