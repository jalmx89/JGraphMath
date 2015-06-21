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
import java.util.Iterator;
import java.util.List;
import java.util.Random;

/**
 * {@code Graph} represents a simple, unweighted labeled
 * <a href="https://en.wikipedia.org/wiki/Graph_(mathematics)">graph</a>.
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
 * // The following code fills a list with every connected labeled graph that 
 * // has fewer than 8 verticies.
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
     * A mask with the lowest {@code (matrixCapacity % 64)} bits set. 
     */
    protected final long highMask;
    
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
     * An array containing the degree sequence.
     * <p>
     * Uses lazy instantiation.
     * @see #getDegreeSequence()
     */
    protected int[] degrees = null;
    
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
    protected Graph[] componentGraphs = null;
    
    
    
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
        matrixCapacity = getMaximumSize(order);
        highMask = ~(~0L << (matrixCapacity & 0x3F));
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
        matrixCapacity = getMaximumSize(order);
        highMask = ~(-1L << (matrixCapacity & 0x3F));
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
        highMask = g.highMask;
        matrix = g.matrix.clone();
        if (g.size != null)
            size = g.size;
        if (g.degrees != null)
            degrees = Arrays.copyOf(g.degrees, g.degrees.length);
        if (g.adjacency != null) {
            adjacency = new int[g.adjacency.length][];
            for (int i=0; i<adjacency.length; i++) {
                adjacency[i] = g.adjacency[i].clone();
            }
        }
        if (g.components != null) {
            components = new int[g.components.length][];
            for (int i=0; i<components.length; i++) {
                components[i] = g.components[i].clone();
            }
        }
    }
    
    /**
     * Constructs a new {@code Graph} that is a subgraph of the specified graph
     * using the specified vertex index map.
     * <p>
     * The vertex at index {@code n} in the newly constructed graph will be 
     * mapped from the vertex at index {@code map[n} in the specified graph.
     */
    public Graph(Graph g, int[] map) {
        this(map.length);
        
        long[] m = matrix;
        long index = 0;
        for (int r=1; r<order; r++) {
            int R = map[r];
            for (int c=0; c<r; c++) {
                if (g.getEdge(R, map[c])) {
                    m[(int)(index>>>6)] |= (1L << index);
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
     * otherwise it will not. This method has no effect if {@code a==b}.
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
     * Signals that the graphs adjacency matrix has been modified.
     * <p>
     * If this method is overridden, the overriding method should invoke
     * {@code super.modified()} so that the subgraphs are still notified about
     * the modification.
     */
    protected void modified() {
        // Clear cached data
        size = null;
        degrees = null;
        adjacency = null;
        components = null;
        componentGraphs = null; 
    }
    
    
    
    /**
     * Randomizes the edges contained in the graph. This method randomizes the
     * entries in the graph's adjacency matrix.
     */
    public void randomize() {
        Random random = new Random();
        for (int i=0; i<matrix.length; i++)
            matrix[i] = random.nextLong();
        matrix[matrix.length-1] &= highMask;
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
        matrix[matrix.length-1] &= highMask;
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
        long[] m = matrix;
        int n = m.length-1;
        for (int i=0; i<n; i++)
            if (++m[i] != 0L)
                return true;
        long v = m[n];
        boolean b = (v != highMask);
        m[n] = (v + 1) & highMask;
        return b;
    }
    
    
    
    /**
     * Represents the order in which bits are read from the lower triangle of
     * a graph's adjacency matrix.
     */
    public static enum ReadOrder {
        /**
         * Left-to-right (top-to-bottom).
         * <pre>{@code
         * 0
         * 1 2
         * 3 4 5}</pre>
         */
        LRTB(true, true, true),
        
        /**
         * Left-to-right (bottom-to-top).
         * <pre>{@code
         * 5
         * 3 4
         * 0 1 2}</pre>
         */
        LRBT(true, true, false),
        
        /**
         * Right-to-left (top-to-bottom).
         * <pre>{@code
         * 0
         * 2 1
         * 5 4 3}</pre>
         */
        RLTB(true, false, true),
        
        /**
         * Right-to-left (bottom-to-top).
         * <pre>{@code
         * 5
         * 4 3
         * 2 1 0}</pre>
         */
        RLBT(true, false, false),
        
        /**
         * Top-to-bottom (left-to-right).
         * <pre>{@code
         * 0
         * 1 3
         * 2 4 5}</pre>
         */
        TBLR(false, true, true),
        
        /**
         * Bottom-to-top (left-to-right).
         * <pre>{@code
         * 2
         * 1 4
         * 0 3 5}</pre>
         */
        BTLR(false, true, false),
        
        /**
         * Top-to-bottom (right-to-left).
         * <pre>{@code
         * 3
         * 4 1
         * 5 2 0}</pre>
         */
        TBRL(false, false, true),
        
        /**
         * Bottom-to-top (right-to-left).
         * <pre>{@code
         * 5
         * 4 2
         * 3 1 0}</pre>
         */
        BTRL(false, false, false);
        
        boolean major;
        boolean rowdir;
        boolean coldir;
        
        ReadOrder(boolean major, boolean rowdir, boolean coldir) {
            this.rowdir = rowdir;
            this.coldir = coldir;
            this.major = major;
        }
        
        /**
         * Returns {@code true} if row-major order and {@code false} if
         * column-major order.
         * 
         * @return {@code true} if row-major order
         */
        public boolean getMajorOrder() {
            return major;
        }
        
        /**
         * Returns {@code true} if rows are read from left to right and
         * {@code false} if rows are read from right to left.
         * 
         * @return {@code true} if rows are read from left to right
         */
        public boolean getRowDirection() {
            return rowdir;
        }
        
        /**
         * Returns {@code true} if columns are read from top to bottom and
         * {@code false} if columns are read from bottom to top.
         * 
         * @return {@code true} if columns are read from top to bottom
         */
        public boolean getColDirection() {
            return rowdir;
        }
    }
    
    /**
     * Returns an iterator that traverses the lower triangle of the graph's
     * adjacency matrix in the order specified by {@code readOrder}.
     * 
     * @param readOrder the order in which the lower triangle of the graph's
     * adjacency matrix is traversed
     * @return an iterator that traverses the lower triangle of the graph's
     * adjacency matrix
     */
    public Iterator<Boolean> iterator(ReadOrder readOrder) {
        if (readOrder.major)
            if (readOrder.rowdir)
                if (readOrder.coldir)
//                    for (int r=1; r<order; r++)
//                        for (int c=0; c<r; c++)
                    return new Iterator<Boolean>() {
                        int r = 1;
                        int c = 0;
                        @Override public boolean hasNext() {
                            return r < order;
                        }
                        @Override public Boolean next() {
                            boolean b = getEdge(r, c);
                            if (++c >= r) {
                                r++;
                                c = 0;
                            }
                            return b;
                        }
                    };
                else 
//                    for (int r=1; r<order; r++)
//                        for (int c=r-1; c>=0; c--)
                    return new Iterator<Boolean>() {
                        int r = 1;
                        int c = r-1;
                        @Override public boolean hasNext() {
                            return r < order;
                        }
                        @Override public Boolean next() {
                            boolean b = getEdge(r, c);
                            if (--c < 0) {
                                r++;
                                c = r-1;
                            }
                            return b;
                        }
                    };
            else 
                if (readOrder.coldir)
//                    for (int r=order; r>=1; r--)
//                        for (int c=0; c<r; c++)
                    return new Iterator<Boolean>() {
                        int r = order;
                        int c = 0;
                        @Override public boolean hasNext() {
                            return r >= 1;
                        }
                        @Override
                        public Boolean next() {
                            boolean b = getEdge(r, c);
                            if (++c >= r) {
                                r--;
                                c = 0;
                            }
                            return b;
                        }
                    };
                else 
//                    for (int r=order; r>=1; r--)
//                        for (int c=r-1; c>=0; c--)
                    return new Iterator<Boolean>() {
                        int r = order;
                        int c = r-1;
                        @Override public boolean hasNext() {
                            return r >= 1;
                        }
                        @Override
                        public Boolean next() {
                            boolean b = getEdge(r, c);
                            if (--c < 0) {
                                r--;
                                c = r-1;
                            }
                            return b;
                        }
                    };
        else 
            if (readOrder.rowdir)
                if (readOrder.coldir)
//                    for (int c=0; c<order-1; c++)
//                        for (int r=c+1; r<order; r++)
                    return new Iterator<Boolean>() {
                        int c = 0;
                        int r = c+1;
                        @Override public boolean hasNext() {
                            return c < order-1;
                        }
                        @Override public Boolean next() {
                            boolean b = getEdge(r, c);
                            if (++r >= order) {
                                c++;
                                r = c+1;
                            }
                            return b;
                        }
                    };
                else 
//                    for (int c=order-2; c>=0; c--)
//                        for (int r=c+1; r<order; r++)
                    return new Iterator<Boolean>() {
                        int c = order-2;
                        int r = c+1;
                        @Override public boolean hasNext() {
                            return c >= 0;
                        }
                        @Override public Boolean next() {
                            boolean b = getEdge(r, c);
                            if (++r >= order) {
                                c--;
                                r = c+1;
                            }
                            return b;
                        }
                    };
            else 
                if (readOrder.coldir)
//                    for (int c=0; c<order-1; c++)
//                        for (int r=order-1; r>c; r--)
                    return new Iterator<Boolean>() {
                        int c = 0;
                        int r = order-1;
                        @Override public boolean hasNext() {
                            return c < order-1;
                        }
                        @Override public Boolean next() {
                            boolean b = getEdge(r, c);
                            if (--r <= c) {
                                c++;
                                r = order-1;
                            }
                            return b;
                        }
                    };
                else 
//                    for (int c=order-2; c>=0; c--)
//                        for (int r=order-1; r>c; r--)
                    return new Iterator<Boolean>() {
                        int c = order-2;
                        int r = order-1;
                        @Override public boolean hasNext() {
                            return c >= 0;
                        }
                        @Override public Boolean next() {
                            boolean b = getEdge(r, c);
                            if (--r <= c) {
                                c--;
                                r = order-1;
                            }
                            return b;
                        }
                    };
    }
    
    /**
     * Returns an iterator that traverses the lower triangle of the graph's
     * adjacency matrix by reading bits from left left-to-right and 
     * top-to-bottom.
     * <p>
     * Equivalent to (though more efficient than) calling
     * {@code iterator(ReadOrder.LRTB)}.
     *
     * @return an iterator that traverses the lower triangle of the graph's
     * adjacency matrix
     */
    @Override
    public Iterator<Boolean> iterator() {
        return new Iterator<Boolean>() {
            long index = 0;
            @Override public boolean hasNext() {
                return index < matrixCapacity;
            }
            @Override public Boolean next() {
                return getEdge(index++);
            }
        };
    }
    
    /**
     * Returns a binary (characters {@code '0'} and {@code '1'}) string
     * representation of the lower triangle of the graphs adjacency matrix using
     * the read order specified.
     * <p>
     * The first bit read, which is determined by the read order, will be stored
     * in index {@code 0} of the returned string.
     * 
     * @param readOrder the order that bits will be read
     * @return a binary string representation of the graph
     */
    public String toString(ReadOrder readOrder) {
        char[] arr = new char[(int)matrixCapacity];
        int idx = 0;
        Iterator<Boolean> it = iterator(readOrder);
        while (it.hasNext())
            arr[idx++] = it.next()? '1' : '0';
        return new String(arr);
    }
    
    /**
     * Returns the bit-packed lower triangle of the graph's adjacency matrix as
     * a little-endian binary string.
     * <p>
     * Equivalent to (thought more efficient than) calling 
     * {@code toString(BitOrder.RLTB).
     *
     * @return a binary string representation of the graph
     */
    @Override
    public String toString() {
        char[] arr = new char[(int)matrixCapacity];
        int idx = arr.length;
        for (Boolean b : this)
            arr[--idx] = b ? '1' : '0';
        return new String(arr);
    }
    
    /**
     * Returns a {@code BigInteger} representation of the lower triangle of the
     * graphs adjacency matrix using the read order specified.
     * <p>
     * Equivalent to {@code new BigInteger(toString(readOrder), 2)}.
     * 
     * @return a {@code BigInteger} representation of the graph
     */
    public BigInteger toBigInteger(ReadOrder readOrder) {
        return new BigInteger(toString(readOrder), 2);
    }
    
    /**
     * Returns a {@code BigInteger} representation of the bit-packed lower
     * triangle of the graphs adjacency matrix.
     * <p>
     * Equivalent to (though more efficient than) 
     * {@code toBigInteger(BitOrder.RLTB)}.
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
     * Prints the upper triangle of the graph's adjacency matrix to the standard
     * output stream.
     */
    public void printUpperTriangle() {
        StringBuilder str = new StringBuilder();
        for (int r=0; r<order-1; r++) {
            int c = 0;
            for (; c<=r; c++)
                str.append("  ");
            for (; c<order; c++)
                str.append(getEdge(r, c) ? "1 ": "0 ");
            str.append("\n");
        }
        System.out.println(str);
    }
    
    /**
     * Prints the lower triangle of the graph's adjacency matrix to the standard
     * output stream.
     */
    public void printLowerTriangle() {
        StringBuilder str = new StringBuilder();
        for (int r=1; r<order; r++) {
            for (int c=0; c<r; c++)
                str.append(getEdge(r, c) ? "1 " : "0 ");
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
     * Returns the
     * <a hfref="https://en.wikipedia.org/wiki/Degree_(graph_theory)#Degree_sequence">
     * degree sequence</a> of the graph.
     * 
     * @return the degree sequence
     */
    public int[] getDegreeSequence() {
        if (degrees == null) {
            int[][] e = getAdjacencyList();
            int len = e.length;
            int[] d = new int[len];
            for (int i=0; i<len; i++)
                d[i] = e[i].length;
            // This is a rubbish way to get decending order, but for the 
            // sorting methods for primative do not accept comparators and 
            // I'm too lazy to write a sorting algorithm just for this.
            Arrays.sort(d);
            for (int i=0; i<len; i++)
                d[i] = d[len-1-i];
            degrees = d;
        }
        return degrees;
    }
    
    /**
     * Returns a jagged array representation of the disjoint components of this
     * graph.
     * <p>
     * If the returned array has length {@code n}, then this graph has {@code n}
     * disjoint components whose union is equal to this graph. If
     * {@code getComponents()[0].length == m}, then the first component of this
     * graph contains {@code m} vertices and their indices (in {@code this}
     * graph) are stored in {@code getComponents()[0]}.
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
            // Continue until all verticies have been visited
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
            comps = Arrays.copyOf(comps, compCount);
            // Sort the components in ascending order by vertex count
            Arrays.sort(comps, (int[] a, int[] b)->(a.length-b.length));
            // Store the components
            components = comps;
        }
        return components;
    }
    
    /**
     * Returns an array of disjoint subgraphs whose union is equal to this 
     * graph.
     * 
     * @return an array of disjoint subgraphs
     */
    public Graph[] getComponentGraphs() {
        if (componentGraphs == null) {
            int[][] comps = getComponents();
            Graph[] subs = new Graph[comps.length];
            for (int i=0; i<subs.length; i++)
                subs[i] = new Graph(this, comps[i]);
            componentGraphs = subs;
        }
        return componentGraphs;
    }
    
    
    
    /**
     * Returns the number of disjoint components in this graph.
     * 
     * @return the number of disjoint components in this graph.
     */
    public int getNumComponents() {
        return getComponents().length;
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
     * Returns {@code true} if the graph is regular. A graph is regular when
     * every vertex has the same degree.
     * 
     * @return {@code true} if the graph is regular
     */
    public boolean isRegular() {
        int[] d = getDegreeSequence();
        int v = d[0];
        for (int i=1; i<d.length; i++)
            if (d[i] != v)
                return false;
        return true;
    }
    
    /**
     * Returns {@code true} if the graph is connected. A graph is connected if
     * there is there is a path from every vertex to every other vertex.
     * 
     * @return {@code true} if the graph is connected
     */
    public boolean isConnected() {
        long s = getSize();
        if (s < order-1)
            return false;
        if (s >= matrixCapacity-order+1)
            return true;
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
     * Returns {@code true} if the graph is connected and contains no cycles.
     * 
     * @return {@code true} if the graph is connected and contains no cycles
     */
    public boolean isTree() {
        return isConnected() && getSize() == order-1;
    }
    
    /**
     * Returns {@code true} if each of the graph's disjoint components is a tree
     * (contains no cycles).
     * 
     * @return {@code true} if each of the graph's disjoint components is a tree
     */
    public boolean isForest() {
        if (isConnected()) {
            return getSize() == order-1;
        } else {
            Graph[] sub = getComponentGraphs();
            for (int i=0; i<sub.length; i++)
                if (!sub[i].isTree())
                    return false;
            return true;
        }
    }
    
    
    
    /**
     * Returns the hash code for this {@code Graph}.
     * 
     * @return the hash code for this {@code Graph}
     */
    @Override
    public int hashCode() {
        int hash = 1;
        long[] m = matrix;
        for (int i=0; i<m.length; i++) {
            long e = m[i];
            hash = 31 * hash + (int)(e ^ (e >>> 32));
        }
        return hash;
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
        long i=0;
        for (int r=1; r<order; r++)
            for (int c=0; c<r; c++)
                if (getEdge(i++) != g.getEdge(map[r], map[c]))
                    return false;
        return true;
    }
    
    /**
     * Returns {@code true} if this graph is isomorphic to the specified graph.
     * 
     * @param g the graph with which to test for isomorphism
     * @return {@code true} if this graph is isomorphic to the specified graph
     */
    public boolean isIsomorphic(Graph g) {
        // Order
        if (g.order != order)
            return false;
        
        // Size
        if (g.getSize() != getSize())
            return false;
        
        // Degree list
        int[] gd = g.getDegreeSequence();
        int[] d = g.getDegreeSequence();
        for (int i=0; i<order; i++)
            if (gd[i] != d[i])
                return false;
            
        // Subcomponents
        int[][] gc = g.getComponents();
        int[][] c = g.getComponents();
        if (gc.length != c.length)
            return false;
        for (int i=0; i<c.length; i++)
            if (gc[i].length != c[i].length)
                return false;
        
        // Last resort, brute force it    
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
     * Returns a new empty graph on the specified number of vertices.
     * 
     * @param order the number of vertices
     * @return new empty graph
     */
    public static Graph makeEmptyGraph(int order) {
        return new Graph(order);
    }
    
    /**
     * Returns a new complete graph on the specified number of vertices.
     * 
     * @param order the number of vertices
     * @return a new complete graph
     */
    public static Graph makeCompleteGraph(int order) {
        Graph g = new Graph(order);
        g.complement();
        return g;
    }
    
    /**
     * Returns a new cycle graph on the specified number of vertices.
     * 
     * @param order the number of vertices
     * @return a new cycle graph
     */
    public static Graph makeCycleGraph(int order) {
        Graph g = new Graph(order);
        for (int i=0; i<order; i++)
            g.setEdge(i, (i+1)%order, true);
        return g;
    }
    
    /**
     * Constructs a new star graph on the specified number of vertices using
     * the vertex with the specified index as the inner vertex of the star.
     *  
     * @param order the number of vertexes in the graph
     * @param index the index of the vertex to use as the inner vertex
     * @return a new star graph
     */
    public static Graph makeStarGraph(int order, int index) {
        Graph g = new Graph(order);
        for (int i=0; i<order; i++)
            g.setEdge(i, index, true);
        return g;
    }
    
    
    
    /**
     * Returns an iterator that generates all labeled graphs with the specified
     * number of vertices.
     * 
     * @param order the number of vertices
     * @return an iterator that generates all labeled graphs with the specified
     * number of vertices
     * @see #increment()
     */
    public static Iterator<Graph> itLabeledGraphs(int order) {
        return new Iterator<Graph>() {
            final Graph g = new Graph(order);
            boolean b = true;
            
            @Override
            public boolean hasNext() {
                return b;
            }

            @Override
            public Graph next() {
                Graph h = new Graph(g);
                b = g.increment();
                return h;
            }
        };
    }
    
    /**
     * Returns a list of all labeled graphs with the specified number of
     * vertices.
     *
     * @param order the number of vertices
     * @return a list of all labeled graphs with the specified number of
     * vertices
     * @see #itLabeledGraphs(int) 
     * @see #increment() 
     */
    public static ArrayList<Graph> getLabledGraphs(int order) {
        ArrayList<Graph> list = new ArrayList();
        Iterator<Graph> it = itLabeledGraphs(order);
        while (it.hasNext())
            list.add(it.next());
        return list;
    }
    
    /**
     * Returns an iterator that generates all connected labeled graphs with the
     * specified number of vertices.
     * <p>
     * Calling the returned iterators {@code next()} method always returns the
     * same {@code Graph} object instance, however the graph is modified between
     * successive calls to {@code next()} by
     * {@link Graph#increment() incrementing} the graph until the next connected
     * graph is found.
     * 
     * @param order the number of vertices
     * @return an iterator that generates all connected labeled graphs with the
     * specified number of vertices
     * @see #isConnected()
     */
    public static Iterator<Graph> itConnectedGraphs(int order) {
        return new Iterator<Graph>() {
            final Graph g = makeStarGraph(order, 1);
            boolean b = true;
            
            @Override
            public boolean hasNext() {
                return b;
            }

            @Override
            public Graph next() {
                Graph h = new Graph(g);
                do {
                    b = g.increment();
                } while (b && !g.isConnected());
                return h;
            }
        };
    }
    
    /**
     * Returns a list of all connected labeled graphs with the specified number
     * of vertices.
     *
     * @param order the number of vertices
     * @return a list of all connected labeled graphs with the specified number
     * of vertices
     * @see #itConnectedGraphs(int) 
     * @see #isConnected() 
     */
    public static ArrayList<Graph> getConnectedGraphs(int order) {
        ArrayList<Graph> list = new ArrayList();
        Iterator<Graph> it = itConnectedGraphs(order);
        while (it.hasNext())
            list.add(it.next());
        return list;
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
     * Returns the number of unlabeled graphs with the specified number of
     * vertices. This is given by the exponential equation
     * {@code 2^(v*(v-1)/2)}, where {@code v} is the number of vertices.
     * 
     *
     * @param order the order of the graph
     * @return the approximate number of possible edge configurations for a
     * graph with the specified number of vertices
     * @throws ArithmeticException if the number of unlabeled graphs with the
     * specified number of vertices is too large to be represented by a
     * {@code long}
     */
    public static long getNumGraphs(long order) {
        long s = getMaximumSize(order);
        if (s > 63)
            throw new ArithmeticException();
        return 1L << s;
    }
    
    /**
     * Returns the approximate number of unlabeled graphs with the specified
     * number of vertices.. This is given by the exponential equation
     * {@code 2^(v*(v-1)/2)}, where {@code v} is the number of vertices.
     *
     * @param order the order of the graph
     * @return the approximate number of possible edge configurations for a
     * graph with the specified number of vertices
     */
    public static double getNumGraphsAprox(long order) {
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
