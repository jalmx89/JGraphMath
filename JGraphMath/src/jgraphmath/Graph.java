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

import java.lang.ref.WeakReference;
import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.LongBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

/**
 * {@code Graph} represents a non-weighted simple
 * <a href="https://en.wikipedia.org/wiki/Graph_(mathematics)">graphs</a> with
 * at least two vertices.
 * <p>
 * At its core {@code Graph} stores the graph it represents using an boolean
 * <a href="https://en.wikipedia.org/wiki/Adjacency_matrix">adjacency matrix
 * </a>. However since the an undirected graph's adjacency matrix is a
 * <a href="https://en.wikipedia.org/wiki/Symmetric_matrix">symmetric matrix
 * </a>, it is stored as a
 * <a href="https://en.wikipedia.org/wiki/Triangular_matrix">triangular matrix
 * </a>. The boolean values of the lower triangle of the adjacency matrix are 
 * bit-packed into a {@code long} array. Though bit-packing reduces the
 * performance of random access to the matrix when compared to a {@code boolean}
 * array, it has advantages: 1) Bit-packing uses one eighth as much memory as a
 * boolean array. 2) Some operations such as finding the graphs complement and
 * testing for equality can be performed extremely quickly using bitwise
 * operations. 3) In theory, given sufficient heap space, much larger graphs can
 * be represented.
 * <p>
 * Example usage:
 * <pre>{@code 
 * // The following code fills a list with every connected graph that has fewer
 * // than 8 verticies. Note that the resultant list will contain many 
 * // isomorphically equivalent graphs.
 * 
 * List<Graph> list = new ArrayList();
 * for (int v = 1; v < 8; v++) {
 *    Graph g = new Graph(v);
 *    do {
 *       if (g.isConnected()) {
 *          list.add(new Graph(h));
 *       }
 *    } while (g.increment());
 * }
 * }</pre>
 * 
 * @author Jeremiah N. Hankins
 */
public class Graph implements Iterable<Boolean> {
    
    /**
     * The number of vertices.
     */
    protected final int order;
    
    /**
     * The number of edges that can be stored in the matrix.
     * @see #getMaximumSize(long)
     */
    protected final long matrixCapacity;
    
    /**
     * The precomputed value of {@code matrixCapacity%64}.
     */
    protected final int capMod64;
    
    /**
     * The bit-packed representation of the lower triangle of the graph's 
     * boolean adjacency matrix.
     */
    protected final long[] matrix;
    
    /** 
    * The number of edges in this graph.
     * <p>
     * Uses lazy instantiation.
     * @see #getSize()
     */
    protected Long size = null;
    
    /**
     * The adjacency list representation of the graph.
     * <p>
     * Uses lazy instantiation.
     * @see #getAdjacencyList()
     */
    protected int[][] adjacency = null;
    
    /**
     * A jagged array representation of the disjoint components of this graph
     * whose union form the graph.
     * <p>
     * Uses lazy instantiation.
     * @see #getComponents()
     */
    protected int[][] components = null;
    
    /**
     * An array of subgraphs representing of the disjoint components of this 
     * graph whose union form the graph. 
     * <p>
     * Uses lazy instantiation.
     * @see #getComponentGraphs()
     */
    protected Subgraph[] componentGraphs = null;
    
    /**
     * A list of weak references to subgraphs of this graph. This list is used
     * to notify subgraphs that have not been collected by the GC about
     * modifications to their supergraph (this graph).
     */
    protected List<WeakReference<Subgraph>> subgraphs;
    
    /** 
     * Constructs a new {@code Graph} with the specified number of vertices and
     * the specified backing array.
     * 
     * @param order the number of vertices
     * @param matrix the bit-packed lower triangle of the adjacency matrix to
     * use
     */
    protected Graph(int order, long[] matrix) {
        this.order = order;
        matrixCapacity = triangular(order-1);
        capMod64 = (int)matrixCapacity & 0x3F;
        this.matrix = matrix;
    }
    
    /**
     * Constructs a new {@code Graph} with the specified number of vertices
     * and no edges.
     * 
     * @param order the number of vertices
     */
    public Graph(int order) {
        this.order = order;
        matrixCapacity = triangular(order-1);
        capMod64 = (int)matrixCapacity & 0x3F;
        matrix = new long[(int)((matrixCapacity+63L)/64L)];
    }
    
    /**
     * Constructs a new {@code Graph} by copying the specified graph.
     *
     * @param g the graph to copy
     */
    public Graph(Graph g) {
        order = g.order;
        matrixCapacity = g.matrixCapacity;
        capMod64 = g.capMod64;
        matrix = Arrays.copyOf(g.matrix, g.matrix.length);
    }
    
    /**
     * Registers a subgraph of this graph so that it can receive modification
     * notifications when this graph is modified.
     * 
     * @param subgraph the subgraph to register
     */
    protected void addSubgraph(Subgraph subgraph) {
        if (subgraphs == null)
            subgraphs = new ArrayList();
        subgraphs.add(new WeakReference(subgraph));
    }
    
    /**
     * Signals that the graphs adjacency matrix has been modified.
     * <p>
     * If this method is overridden, the overriding method should invoke
     * {@code super.modified()} so that the subgraphs are still notified about
     * the modification.
     */
    protected void modified() {
        // Clear cached data
        size = null;
        adjacency = null;
        components = null;
        componentGraphs = null; 
        
        // If there are subgraphs attached to this graph
        if (subgraphs != null) {
            // For every supergraph reference
            for (int i=subgraphs.size()-1; i>=0; i--) {
                // Try to get the supergraph
                Subgraph subgraph = subgraphs.get(i).get();
                // If the supergraph has been garbage collected
                if (subgraph == null) {
                    // Remove it from the list
                    subgraphs.remove(i);
                } else {
                    // Otherise, signal that the supergraph has been modified
                    subgraph.modified();
                }
            }
        }
    }
    
    /**
     * Returns the number of vertices in the graph.
     * 
     * @return the number of vertices
     */
    public int getOrder() {
        return order;
    }
    
    /**
     * Returns the maximum number of edges this graph of this order can contain.
     * 
     * @return the maximum size of a graph of this order
     */
    public long getMaximumSize() {
        return matrixCapacity;
    }
    
    /**
     * Returns {@code true} if the graph has an edge at the specified index.
     * Edge indexes are taken from the lower triangle of the graph's adjacency
     * matrix ordered from left to right and top to bottom.
     * <p>
     * For example, the indexes of the edges in the adjacency matrix of a graph
     * with five vertices are: <pre>
     *   - 0 1 3 6
     *   0 - 2 4 7
     *   1 2 - 5 8
     *   3 4 5 - 9
     *   6 7 8 9 -
     * </pre>
     * 
     * @param index the index of the edge being tested
     * @return {@code true} if the graph has an edge at the specified index
     */
    public boolean getEdge(long index) {
        return ((matrix[(int)(index>>>6)] & (1L << index))) != 0;
    }
    
    /**
     * Returns {@code true} if the graph has an edge from the vertex at index 
     * {@code a} to the vertex at index {@code b}.
     * 
     * @param a the index of the first vertex
     * @param b the index of the second vertex
     * @return {@code true} if there is an edge connecting the two vertices
     */
    public boolean getEdge(int a, int b) {
        if (a > b) // lower triangle
            return getEdge((long)a*(a-1)/2+b);
        if (a < b) // upper triangle
            return getEdge((long)b*(b-1)/2+a);
        return false;
    }
    
    /**
     * Sets the presence of the edge at the specified index in the graph. If the
     * {@code value} argument is {@code true} the edge will contain the
     * specified edge after this method call, otherwise it will not.
     * <p>
     * Edge indexes are taken from the lower triangle of the graph's adjacency
     * matrix ordered from left to right and top to bottom.
     * <p>
     * For example, the indexes of the edges in the adjacency matrix of a graph
     * with five vertices are: <pre>
     *   - 0 1 3 6
     *   0 - 2 4 7
     *   1 2 - 5 8
     *   3 4 5 - 9
     *   6 7 8 9 -
     * </pre>
     * 
     * @param index the index of the edge being set
     * @param value determines whether or not the edge will be contained in
     * the graph after the call
     */
    public void setEdge(long index, boolean value) {
        if (value)
            matrix[(int)(index>>>6)] |= (1L << index);
        else
            matrix[(int)(index>>>6)] &= ~(1L << index);
        modified();
    }
    
    /**
     * Sets the presence of the edge from the vertex at index {@code a} to the
     * vertex at index {@code b}. If the {@code value} argument is {@code true}
     * the edge will contain the specified edge after this method call,
     * otherwise it will not.
     *
     * @param a the index of the first vertex
     * @param b the index of the second vertex
     * @param value if {@code true} the graph will contain the edge after the
     * call, otherwise it will not
     */
    public void setEdge(int a, int b, boolean value) {
        if (a > b)      // lower triangle
            setEdge((long)a*(a-1)/2+b, value);
        else if (a < b) // upper triangle
            setEdge((long)b*(b-1)/2+a, value);
    }
    
    /**
     * Randomizes the edges contained in the graph. This method randomizes the
     * entries in the graph's adjacency matrix.
     */
    public void randomize() {
        Random random = new Random();
        for (int i=0; i<matrix.length-1; i++)
            matrix[i] = random.nextLong();
        matrix[matrix.length-1] = (random.nextLong() >>> (64-capMod64));
        modified();
    }
    
    /**
     * Modifies this graph so that it represents its complement after the call.
     * If two vertices share an edge before the call, they will not share an
     * edge after the call. If two vertices did not share an edge before the
     * call, they will share an edge after the call.
     */
    public void complement() {
        for (int i=0; i<matrix.length; i++)
            matrix[i] = ~matrix[i];
        matrix[matrix.length-1] &= (0xFFFFFFFFFFFFFFFFL >>> (64-capMod64));
        modified();
    }
    
    /**
     * Modifies this graph in such a way that the value returned by 
     * {@link #toBigInteger()} will be incremented by one after the call.
     * By repeatedly invoking this method it is possible to generate every graph
     * (i.e. every possible edge configuration) for a graph of this
     * {@link #getOrder() order}. Returns {@code false} if the graph 
     * "rolls over" and transitions from a complete graph to an empty graph,
     * otherwise {@code true}.
     * 
     * @return {@code false} if the graph transitions from a complete graph to
     * an empty graph, otherwise {@code true}
     */
    public boolean increment() {
        modified();
        int n = matrix.length-1;
        for (int i=0; i<n; i++) {
            matrix[i]++;
            if (matrix[i] == 0L) {
                return true;
            }
        }
        matrix[n]++;
        if ((matrix[n] >>> capMod64) != 0) {
            matrix[n] = 0;
            return false;
        }
        return true;
    }
    
    /**
     * Returns an iterator that traverses the lower triangle of the graph's
     * adjacency matrix by reading bits from the lower triangle from left to
     * right and top to bottom, using the first bit read as the least.
     * 
     * @return an iterator that traverses the graph's adjacency matrix
     */
    @Override
    public Iterator<Boolean> iterator() {
        return new Iterator<Boolean>() {
            long index = 0;
            @Override
            public boolean hasNext() {
                return index < matrixCapacity;
            }
            @Override
            public Boolean next() {
                return getEdge(index++);
            }
        };
    }
    
    /**
     * Converts the lower triangle of the graph's adjacency matrix into a binary
     * string by reading bits from the lower triangle from left to right and top
     * to bottom so that the first character in the string represents the
     * first bit read.
     * 
     * @return a binary string representation of the graph
     */
    public String toBinString() {
        char[] arr = new char[(int)matrixCapacity];
        int idx = 0;
        for (Boolean b : this)
            arr[idx++] = b ? '1' : '0';
        return new String(arr);
    }
    
    /**
     * Converts the lower triangle of the graph's adjacency matrix into a 
     * decimal string by reading bits from the lower triangle from left to right
     * and top to bottom so that the least significant bit in the binary 
     * representation of the returned decimal number represents the
     * first bit read.
     * <p>
     * Equivalent to {@code toBigInteger().toString()}.
     * 
     * @return a decimal string representation of the graph
     */
    public String toDecString() {
        return toBigInteger().toString();
    }
    
    /**
     * Converts the lower triangle of the graph's adjacency matrix into a
     * {@code boolean} array by reading bits from the lower triangle from left
     * to right and top to bottom so that the entry at index 0 in the array is
     * the first value read.
     *
     * @return a {@code boolean} array representation of the graph
     */
    public boolean[] toBooleanArray() {
        boolean[] arr = new boolean[(int)matrixCapacity];
        int idx = 0;
        for (Boolean b : this) {
            arr[idx++] = b;
        }
        return arr;
    }

    /**
     * Converts the lower triangle of the graph's adjacency matrix into a
     * {@code BitSet} by reading bits from the lower triangle from left to right
     * and top to bottom, using the first bit read as the least significant bit.
     *
     * @return a {@code BitSet} representation of the graph
     */
    public BitSet toBitSet() {
        return BitSet.valueOf(matrix);
    }

    /**
     * Converts the lower triangle of the graph's adjacency matrix into a
     * {@code long} array by reading bits from the lower triangle from left to
     * right and top to bottom. Values are bit packed into the {@code long}
     * values such that the first bit read is stored in the least significant
     * bit of {@code long} at index 0 in the returned array.
     *
     * @return a {@code long} array representation of the graph
     */
    public long[] toLongArray() {
        return Arrays.copyOf(matrix, matrix.length);
    }

    /**
     * Converts the lower triangle of the graph's adjacency matrix into a
     * {@code BigInteger} by reading bits from the lower triangle from left to
     * right and top to bottom, using the first bit read as the least
     * significant bit.
     *
     * @return a {@code BigInteger} representation of the graph
     */
    public BigInteger toBigInteger() {
        ByteBuffer bb = ByteBuffer.allocate(matrix.length*8);
        LongBuffer lb = bb.asLongBuffer();
        for (int i = matrix.length-1; i >= 0; i--)
            lb.put(matrix[i]);
        return new BigInteger(1, bb.array());
    }
    
    /**
     * Prints the lower triangle of the graph's adjacency matrix to the standard
     * output stream.
     */
    public void printLowerTriangle() {
        StringBuilder str = new StringBuilder();
        for (int row = 1; row < order; row++) {
            for (int col = 0; col < row; col++)
                str.append(getEdge(row, col)?"1 ":"0 ");
            str.append("\n");
        }
        System.out.println(str);
    }
    
    /**
     * Prints the upper triangle of the graph's adjacency matrix to the standard
     * output stream.
     */
    public void printUpperTriangle() {
        StringBuilder str = new StringBuilder();
        for (int row = 0; row < order-1; row++) {
            int col = 0;
            for (; col <= row; col++)
                str.append("  ");
            for (; col < order; col++)
                str.append(getEdge(row, col)?"1 ":"0 ");
            str.append("\n");
        }
        System.out.println(str);
    }
    
    /**
     * Prints the graph's adjacency matrix to the standard output stream.
     */
    public void printMatrix() {
        StringBuilder str = new StringBuilder();
        for (int row = 0; row < order; row++) {
            for (int col = 0; col < order; col++)
                str.append(getEdge(row, col)?"1 ":"0 ");
            str.append("\n");
        }
        System.out.println(str);
    }
    
    /**
     * Returns the number of edges in this graph.
     * 
     * @return the number of edges
     */
    public long getSize() {
        if (size == null) {
            long[] m = matrix;
            long total = 0;
            for (int i=0; i<m.length; i++)
                total += countSetBits(m[i]);
            size = total;
        }
        return size;
    }
    
    /**
     * Returns {@code true} if the graph has no edges.
     * 
     * @return {@code true} if the graph has no edges.
     */
    public boolean isEmpty() {
        return getSize() == 0;
    }
    
    /**
     * Returns {@code true} if the graph is a complete graph. A graph is
     * complete if every vertex is connected to every other vertex by an edge.
     * 
     * @return {@code true} if the graph is a complete graph
     */
    public boolean isComplete() {
        return getSize() == getMaximumSize();
    }
    
    /**
     * Returns an adjacency list representation of the graph.
     * <p>
     * The returned array is a jagged array of integers. The length of the array
     * is equal to the number of vertices, such that each vertex has its own
     * sub-array. The sub-array contains the indexes of the vertices that share
     * an edge with the vertex that owns the sub-array.
     * <p>
     * For example, if the graph has five vertices, then the returned array,
     * {@code edges}, then {@code edges.length} will be {@code 5}. If the vertex
     * at index {@code 2} within the graph has three edges, then {@code edges[2].length}
     * will be {@code 3}. If vertex {@code 4} is connected to vertices 
     * {@code 0}, {@code 1}, and {@code 3}, then {@code edges[4]} will be the
     * list {@code {0, 1, 3}}.
     * 
     * @return an edge list representation of the graph
     */
    public int[][] getAdjacencyList() {
        if (adjacency == null) {
            int[][] e = new int[order][];
            int row[] = new int[order];
            for (int r=0; r<order; r++) {
                int d = 0;
                for (int c=0; c<order; c++)
                    if (getEdge(r, c))
                        row[d++] = c;
                e[r] = Arrays.copyOf(row, d);
            }
            adjacency = e;
        }
        return adjacency;
    }
    
    /**
     * Returns a jagged array representation of the disjoint components of this
     * graph.
     * <p>
     * If the returned array has length {@code n}, then this graph has {@code n}
     * disjoint components whose union is equal to this graph. If
     * {@code getComponents()[0].length == m}, then the first component of this
     * graph contains {@code m} vertices and their indexes (in {@code this}
     * graph) are stored in {@code getComponents()[0]}, and so forth.
     *
     * @return a jagged array representation of the disjoint components of this
     * graph
     */
    public int[][] getComponents() {
        if (components == null) {
            int[][] edges = getAdjacencyList();
            // The total number of components
            int compCount = 0;
            // The array of components
            int[][] comps = new int[order][];
            // The current component
            int[] currComp = new int[order];
            // The vertex queue
            int[] queue = new int[order];
            int queueFront = 0;
            int queueBack = 0;
            // Keep track of which vertexes have been visted
            boolean[] visited = new boolean[order];
            // Keep track of the next vertex that has not been visted
            int visitPos = 0;
            // Continue untill all verticies have been visited
            while (queueFront < order) {
                // The size of the current compoennt
                int currCompSize = 0;
                // Find the next unvisited verted
                while (visitPos < order && visited[visitPos])
                    visitPos++;
                // Mark it as visted
                visited[visitPos] = true;
                // Add it to the queue
                queue[queueBack++] = visitPos;
                // Until the queue is empty
                while (queueFront < queueBack) {
                    // Get the next vertex from the front of the queue
                    int v = queue[queueFront++];
                    // Add the vertex to the componenet
                    currComp[currCompSize++] = v;
                    // Get the array of vertices adjacent to this vertex
                    int[] edgesV = edges[v];
                    // Get the number of adjacent vertices
                    int len = edgesV.length;
                    // For each adjacent vertex
                    for (int i=0; i<len; i++) {
                        int w = edgesV[i];
                        // If the adjacent vertex has not been visted...
                        if (!visited[w]) {
                            // Enqueue the connected vertex
                            queue[queueBack++] = w;
                            // And mark it as visted
                            visited[w] = true;
                        }
                    }
                }
                // Store the completed component
                comps[compCount++] = Arrays.copyOf(currComp, currCompSize);
            }
            // Truncate the components arrays
            components = Arrays.copyOf(comps, compCount);
        }
        return components;
    }
    
    /**
     * Returns an array of disjoint subgraphs whose union is equal to this 
     * graph.
     * 
     * @return an array of disjoint subgraphs
     */
    public Subgraph[] getComponentGraphs() {
        if (componentGraphs == null) {
            int[][] comps = getComponents();
            Subgraph[] subs = new Subgraph[comps.length];
            for (int i=0; i<subs.length; i++)
                subs[i] = new Subgraph(this, comps[i]);
            componentGraphs = subs;
        }
        return componentGraphs;
    }
    
    /**
     * Returns {@code true} if the graph is connected. A graph is connected if
     * there is there is a path from every vertex to every other vertex.
     * 
     * @return {@code true} if the graph is connected
     */
    public boolean isConnected() {
        return getComponents().length == 1;
    }
    
    /**
     * Returns {@code true} if the graph is a cycle graph.
     * 
     * @return {@code true} if the graph is a cycle graph
     */
    public boolean isCycleGraph() {
        return isConnected() && getSize() == order;
    }
    
    /**
     * Returns {@code true} if the specified object is a {@code Graph} with the
     * same number of vertices as this graph and its adjacency matrix is
     * identical to the adjacency matrix of this graph. This method does not
     * check for isomorphic equivalence.
     * 
     * @param obj the reference object with which to compare
     * @return {@code true} if {@code obj} is a graph whose adjacency matrix is
     * identical to this graphs adjacency matrix
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == null)
            return false;
        if (!(obj instanceof Graph))
            return false;
        Graph g = (Graph)obj;
        if (g.order != order)
            return false;
        for (int i=0; i<matrix.length; i++)
            if (g.matrix[i] != matrix[i])
                return false;
        return true;
    }
    
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 79 * hash + this.order;
        hash = 79 * hash + Arrays.hashCode(this.matrix);
        return hash;
    }
    
    /**
     * Returns {@code true} if the this graph is equal to the specified graph
     * under the provided vertex mapping. 
     * <p>
     * In more precise terms, if for every pair of vertex indices {@code a} and 
     * {@code b} it holds that 
     * {@code (this.getEdge(a, b) == g.getEdge(map[a], map[b])}, then this
     * method will return {@code true}, otherwise it will return {@code false}.
     * This method also returns {@code false} if the order of the provided graph
     * is not equal to the order of this graph.
     * 
     * @param g the graph to test equality against
     * @param map the vertex map
     * @return {@code true} if the this graph is equal to the specified graph
     * under the provided vertex mapping
     */
    public boolean equals(Graph g, int[] map) {
        if (g.order != order)
            return false;
        int i=0;
        for (int r=1; r<order; ++r)
            for (int c=0; c<r; ++c)
                if (getEdge(i++) != g.getEdge(map[r], map[c]))
                    return false;
        return true;
    }
    
    /**
     * Returns {@code true} if this graph is isomorphic to the specified graph.
     * <p>
     * This method is implemented using the Johnson Trotter permutation
     * algorithm.
     * 
     * @param g the graph with which to test for isomorphism
     * @return {@code true} if this graph is isomorphic to the specified graph
     * @see #getIsomorphicMapJohnsonTrotter(jgraphmath.Graph)
     */
    public boolean isIsomorphicMapJohnsonTrotter(Graph g) {
        return getIsomorphicMapJohnsonTrotter(g) != null;
    }
    
    /**
     * Returns an isomorphic mapping using the Johnson Trotter permutation
     * algorithm, or {@code null} if no isomorphic mapping exists.
     * <p>
     * If an isomorphism exists between this graph and the specified graph, then
     * the returned map will contain an isomorphic mapping from the vertex 
     * indexes in the specified map to the vertex indices in this map.
     * 
     * @param g the graph with which to find the isomorphic mapping
     * @return an isomorphic mapping, or {@code null}
     */
    public int[] getIsomorphicMapJohnsonTrotter(Graph g) {
        if (g.order != order)
            return null;
        int[] prm = new int[order];     // permutation
        int[] inv = new int[order];     // inverse permutation
        int[] dir = new int[order];     // direction = +1 or -1
        for(int i=0; i<order; i++) {
            dir[i] = -1;
            prm[i] = i;
            inv[i] = i;
        }
        return getIsomorphicMapJohnsonTrotter(g, 0, prm, inv, dir);
    }
    
    /**
     * Recursive helper method used by
     * {@link #getIsomorphicMapJohnsonTrotter(Graph.Graph)}.
     *
     * @param g the graph with which to find the isomorphic mapping
     * @param n the "position" of the algorithm in the map
     * @param prm the current permutation
     * @param inv the inverse permutation
     * @param dir the current "direction" of the algorithm
     * @return an isomorphic mapping, or {@code null}
     * @see #getIsomorphicMapJohnsonTrotter(jgraphmath.Graph, int, int[], int[], int[]) 
     */
    private int[] getIsomorphicMapJohnsonTrotter(Graph g, int n, int[] prm, int[] inv, int[] dir) { 
        if (prm.length <= n) {
            if (equals(g, prm)) {
                return prm;
            } else {
                return null;
            }
        }
        int[] map = getIsomorphicMapJohnsonTrotter(g, n+1, prm, inv, dir);
        if (map != null) {
            return map;
        }
        for(int i=0; i<=n-1; i++) {
            int z = prm[inv[n] + dir[n]];
            prm[inv[n]] = z;
            prm[inv[n] + dir[n]] = n;
            inv[z] = inv[n];
            inv[n] = inv[n] + dir[n];  
            map = getIsomorphicMapJohnsonTrotter(g, n+1, prm, inv, dir);
            if (map != null) {
                return map;
            }
        }
        dir[n] = -dir[n];
        return null;
    }
    
    /**
     * Returns the n-th term in the sequence of triangular numbers. This is
     * sequence <a href="https://oeis.org/A000217">A000217 in OEIS</a>. The n-th
     * term in the sequence is given by the formula {@code n*(n+1)/1}. Beginning
     * with the 0th triangle number, the sequence is 0, 1, 3, 6, 10, 15, 21, 28,
     * 36, 45, 55, 66, 78, 91, 105, 120, 136, 153, 171, 190, 210, 231, 253, 276,
     * 300, 325, 351, 378, 406...
     *
     * @param n the index of the triangular number
     * @return the n-th triangular number
     */
    public static long triangular(long n) {
        return n*(n+1)/2;
    }
    
    /**
     * Returns the index of the specified nonnegative triangular number in
     * the sequence of triangular numbers, or {@code -1} if the specified number
     * is negative or not triangular.
     * 
     * @param n the triangular number
     * @return the position of the triangular number in the sequence, or {@code -1}
     */
    public static long triangularInv(long n) {
        long r = (long)Math.sqrt(8*n+1);
        return (r*r == n)? (r-1)/2: -1;
    }
    
    /**
     * Returns the index of the largest triangular number that is less than the
     * specified number.
     * 
     * @param n the number
     * @return the index of the largest triangular number that is less than the
     * specified number
     */
    public static long trinalgeInvFloor(long n) {
        return ((long)Math.sqrt(8*n+1)-1)/2;
    }
    
    /**
     * Returns the largest triangular number that is less than the specified
     * number.
     * 
     * @param n the number
     * @return the largest triangular number that is less than the specified
     * number
     */
    public static long triangularFloor(long n) {
        return triangular(trinalgeInvFloor(n));
    }
    
    /**
     * Returns the maximum number of edges a graph with the specified number of
     * vertices could have. In other words, returns the number of edges in a
     * complete graph of the specified order.
     *
     * @param order the number of vertices
     * @return the number of edges in a complete graph of the specified order
     */
    public static long getMaximumSize(long order) {
        return triangular(order-1);
    }
    
    /**
     * Returns the approximate number of possible edge configurations for a
     * graph with the specified number of vertices. This is given by the
     * exponential equation {@code 2^(v*(v-1)/2)}, where {@code v} is the number
     * of vertices.
     *
     * @param order the order of the graph
     * @return the approximate number of possible edge configurations for a
     * graph with the specified number of vertices
     */
    public static double getNumGraphs(long order) {
        return Math.pow(2, getMaximumSize(order));
    }
    
    /**
     * Returns the number of bits that are set in the specified 64 bit integer.
     * In other words, this method counts the number of {@code 1}s in the binary
     * representation of the specified number.
     * <p>
     * For example, {@code 151} can be written in binary as {@code 1001011}. So,
     * if {@code 151} were supplied as the argument, this method would return
     * {@code 5}.
     * 
     * @param i the number whose set bits will be counted
     * @return the number of bits set in the specified {@code long}
     */
    public static long countSetBits(long i) {
        i = i - ((i >>> 1) & 0x5555555555555555L);
        i = (i & 0x3333333333333333L) + ((i >>> 2) & 0x3333333333333333L);
        return (((i + (i >>> 4)) & 0xF0F0F0F0F0F0F0FL) * 0x101010101010101L) >>> 56;
    }
    
    /**
     * Reverses the order of the bits in the specified 64 bit integer.
     * 
     * @param x the original number
     * @return the bit order reversed number
     */
    public static long reverseBits(long x) {
        x = (((x & 0xAAAAAAAAAAAAAAAAL) >>> 1) | ((x & 0x5555555555555555L) << 1));
        x = (((x & 0xCCCCCCCCCCCCCCCCL) >>> 2) | ((x & 0x3333333333333333L) << 2));
        x = (((x & 0xF0F0F0F0F0F0F0F0L) >>> 4) | ((x & 0x0F0F0F0F0F0F0F0FL) << 4));
        x = (((x & 0xFF00FF00FF00FF00L) >>> 8) | ((x & 0x00FF00FF00FF00FFL) << 8));
        return((x >> 16) | (x << 16));
    }
}
