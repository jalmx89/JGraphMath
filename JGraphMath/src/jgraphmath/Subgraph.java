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

package jgraphmath;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Random;

/**
 * {@code Subgraph} is a
 * <a href="https://en.wikipedia.org/wiki/Glossary_of_graph_theory#Subgraphs">
 * subgraph</a> of a {@link Graph}. Subgraphs are backed by their supergraph so
 * that changes to the subgraph will be reflected in the supergraph and vice
 * versa.
 * <p>
 * Though {@code Subgraph} extends {@code Graph} and provides all of the same
 * methods and functionality, {@code Subgraph} accesses the data contained in a
 * {@code Graph} using an vertex index mapping function and generally will not
 * perform as well as a {@code Graph}. However the performance of a subgraph of
 * a subgraph will not be further reduced. If {@code A} is a subgraph of 
 * {@code Graph} {@code G}, and {@code B} is a subgraph of {@code A}, then
 * the supergraph of {@code B} will be {@code G}, not {@code A}.
 * 
 * @author Jeremiah N. Hankins
 */
public class Subgraph extends Graph {
    /**
     * The supergraph that contains this subgraph.
     */
    private final Graph supergraph;
    
    /**
     * The mapping from this subgraph's vertex indices to the supergraph's
     * vertex indices.
     */
    private final int[] supergraphMap;
    
    /**
     * Constructs a new {@code Subgraph} backed by the specified
     * supergraph that uses the specified mapping from this graph's vertex
     * indices to the supergraph's vertex indices.
     * <p>
     * If a vertex has index {@code a} in this subgraph, then the index of that
     * vertex in the supergraph is {@code supergraphMap[a]}.
     * 
     * @param supergraph the graph containing the new graph
     * @param supergraphMap the vertex index map
     */
    @SuppressWarnings("LeakingThisInConstructor")
    public Subgraph(Graph supergraph, int[] supergraphMap) {
        super(supergraphMap.length, null);
        this.supergraph = supergraph;
        this.supergraphMap = Arrays.copyOf(supergraphMap, supergraphMap.length);
        supergraph.addSubgraph(this);
    }
    
    /**
     * Returns the graph that contains this graph.
     * 
     * @return the graph that contains this graph
     */
    public Graph getSupergraph() {
        return supergraph;
    }

    /**
     * Called by the supergraph to inform this subgraph that the supergraph has
     * been modified.
     */
    protected void supergraphModified() {
        super.modified();
    }
    
    @Override
    public void modified() {
        supergraph.modified();
    }
    
    @Override
    public boolean getEdge(long index) {
        int a = (int)trinalgeInvFloor(index);
        int b = (int)(index-triangular(a));
        a = supergraphMap[a];
        b = supergraphMap[b];
        return supergraph.getEdge(a, b);
    }
    
    @Override
    public void setEdge(long index, boolean value) {
        int a = (int)trinalgeInvFloor(index);
        int b = (int)(index-triangular(a));
        a = supergraphMap[a];
        b = supergraphMap[b];
        supergraph.setEdge(a, b, value);
        modified();
    }
    
    @Override
    public long getSize() {
        if (size == null) {
            long z = 0;
            long s = getMaximumSize();
            for (long i=0; i<s; i++) {
                if (getEdge(i))
                    z++;
            }
            size = z;
        }
        return size;
    }
    
    @Override
    public void randomize() {
        Random rand = new Random();
        long s = getMaximumSize();
        for (long i=0; i<s; i++) {
            setEdge(i, rand.nextBoolean());
        }
    }
    
    @Override
    public void complement() {
        long s = getMaximumSize();
        for (int i=0; i<s; i++) {
            setEdge(i, !getEdge(i));
        }
    }
    
    @Override
    public boolean increment() {
        long s = getMaximumSize();
        for(int i=0; i<s; i++) {
            if(getEdge(i)) {
                setEdge(i, false);
            } else {
                setEdge(i, true);
                return true;
            }
        }
        return false;
    }
    
    @Override
    public BitSet toBitSet() {
        BitSet bitset = new BitSet();
        int index = 0;
        for (Boolean b : this)  {
            if (b)
                bitset.set(index);
            index++;
        }
        return bitset;
    }

    @Override
    public BigInteger toBigInteger() {
        BigInteger bigint = BigInteger.ZERO;
        int index = 0;
        for (Boolean b : this)  {
            if (b)
                bigint.setBit(index);
            index++;
        }
        return bigint;
    }
    
    @Override
    public Subgraph[] getComponentGraphs() {
        if (componentGraphs == null) {
            int[][] comps = getComponents();
            Subgraph[] subs = new Subgraph[comps.length];
            for (int i=0; i<subs.length; i++) {
                int[] c = comps[i];
                int[] q = new int[c.length];
                for (int j=0; j<c.length; j++)
                    q[j] = supergraphMap[q[j]];
                subs[i] = new Subgraph(supergraph, comps[i]);
            }
            componentGraphs = subs;
        }
        return componentGraphs;
    }
}
