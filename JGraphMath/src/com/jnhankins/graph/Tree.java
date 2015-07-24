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

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Stack;

/**
 * A rooted tree on a graph.
 * <pre>
 * Example usage:
 * {@code 
 * 
 * @author Jeremiah N. Hankins
 */
public class Tree {
    
    /**
     * The graph that this tree is on. (The graph that was used to construct
     * this tree.)
     */
    protected final Graph graph;
    
    /**
     * The root node.
     */
    protected final TreeNode root;
    
    /**
     * A list of unique vertex indices contained in this tree sorted in
     * ascending order.
     */
    protected final int[] coverage;
    
    /**
     * The number of nodes in the tree.
     */
    protected final int size;
    
    /**
     * Constructs a new {@code Tree}. The {@code coverage} array is stored by
     * direct reference, i.e. not copied, and is sorted by this constructor.
     * 
     * @param graph the graph that this tree is on
     * @param root the root node
     * @param coverage a list of unique vertex indices contained in this tree
     * @param size the number of nodes in the tree
     */
    protected Tree(Graph graph, TreeNode root, int[] coverage, int size) {
        this.graph = graph;
        this.root = root;
        this.coverage = coverage;
        this.size = size;
        Arrays.sort(coverage);
    }
    
    /**
     * Returns the graph that this tree is on (i.e. the graph that was used to
     * construct the tree). Note that the graph or this may have been modified
     * since this tree was constructed.
     * 
     * @return the graph that this tree is on
     */
    public Graph getGraph() {
        return graph;
    }
    
    /**
     * Returns the root node of the tree.
     * 
     * @return the root node of the tree
     */
    public TreeNode getRootNode() {
        return root;
    }
    
    /**
     * Returns the number of nodes in the tree.
     * 
     * @return return size
     */
    public int getSize() {
        return size;
    }
    
    /**
     * Returns the number of unique vertex indices in the tree.
     * 
     * @return the number of unique vertex indices in the tree
     */
    public int getCoverage() {
        return coverage.length;
    }
    
    /**
     * Returns {@code true} if this tree is connected to the specified tree. In
     * other words, returns {@code true} if this tree and the specified tree are
     * both on the same graph and both contain at least one vertex in common.
     * 
     * @param t the tree to test for connectivity with
     * @return {@code true} if this tree is connected to the specified tree
     */
    public boolean isConnected(Tree t) {
        if (this.getGraph() !=  t.getGraph())
            return false;
        int[] coverA = coverage;
        int[] coverB = t.coverage;
        int a = coverA[0];
        int b = coverB[0];
        for (int i=1, j=1; i<coverA.length && j<coverB.length; ) {
            if      (a < b) a = coverA[i++];
            else if (b < a) b = coverB[j++];
            else return true;
        }
        return false;
    }
    
    /**
     * Returns canonical number of this tree as an unsigned little-endian
     * integer packed into a {@code long} array. If the tree was not in
     * canonical form before this method was called, it will be after. Use
     * {@link #getCanonicalBitLength()} to determine how many bits in the
     * returned {@code long} array are being used to represent the canonical
     * number.
     *
     * @return the canonical tree number
     * @see #getCanonicalBitLength() 
     */
    public long[] getCanonicalNumber() {
        return root.getCanonicalNumber();
    }

    /**
     * Returns the number of bits being used to represent this tree's
     * canonical number, or -1 if this node's subtree has not yet been
     * transformed into its canonical form.
     * 
     * @return the number of bits need to represent the the canonical
     * number, or {@code -1} if unknown
     * @see #getCanonicalNumber() 
     */
    public int getCanonicalBitLength() {
        return root.getCanonicalBitLength();
    }
        
    /**
     * Returns an iterator that traverses the nodes in the tree using pre-order
     * traversal. The behavior of the returned iterator is undefined if the tree
     * is modified before the iteration has completed (i.e. the iterator is not
     * fast-fail and does not check for concurrent modifications).
     * 
     * @return a pre-order traversal tree node iterator
     */
    public Iterator<TreeNode> iteratorPreOrder() {
        return root.iteratorPreOrder();
    }

    /**
     * Returns an iterator that traverses the nodes in the tree using pre-order
     * traversal. The behavior of the returned iterator is undefined if the tree
     * is modified before the iteration has completed (i.e. the iterator is not
     * fast-fail and does not check for concurrent modifications).
     *
     * @return a post-order traversal tree node iterator
     */
    public Iterator<TreeNode> iteratorPostOrder() {
        return root.iteratorPostOrder();
    }

    /**
     * Returns an iterator that traverses the tree using level-order (aka
     * breadth-first search). The behavior of the returned iterator is undefined
     * if the tree is modified before the iteration has completed (i.e. the
     * iterator is not fast-fail and does not check for concurrent
     * modifications).
     *
     * @return a level-order traversal tree node iterator
     */
    public Iterator<TreeNode> iteratorLevelOrder() {
        return root.iteratorLevelOrder();
    }
    
    @Override
    public String toString() {
        return root.toString();
    }
    
    /**
     * A node wrapping a vertex in a {@code Tree}.
     * 
     * @author Jeremiah N. Hankins
     */
    public static class TreeNode {
        /**
         * The index of the vertex wrapped by this tree node.
         */
        protected final int vertex;
        
        /**
         * Child nodes.
         */
        protected TreeNode[] children;
        
        /**
         * The canonical number of this node's subtree stored as an unsigned
         * little-endian integer packed into a {@code long} array, or
         * {@code null} if the subtree has not yet been transformed into its
         * canonical form.
         */
        protected long[] canonicalNumber = null;
        
        /**
         * The number of bits being used to represent this subtree's canonical
         * number, or -1 if this node's subtree has not yet been transformed
         * into its canonical form.
         */
        protected int canonicalBits = -1;
        
        /**
         * Constructs a new {@code TreeNode}.
         * 
         * @param vertex index of the vertex
         */
        protected TreeNode(int vertex) {
            this.vertex = vertex;
        }
        
        /**
         * Returns the index of the vertex that this tree node is wrapping. The
         * index is the index assigned to the vertex in the graph that the tree
         * is on.
         * 
         * @return returns the index of the vertex wrapped by this node
         */
        public int getVertexIndex() {
            return vertex;
        }
        
        /**
         * Returns the number of nodes that are children of this node.
         * 
         * @return the number of child nodes
         */
        public int getChildCount() {
            return children == null? 0 : children.length;
        }
        
        /**
         * Returns the child node with the specified index.
         * 
         * @param childIndex the index of the child to get
         * @return the child node with the specified index
         * @throws IndexOutOfBoundsException if {@code childIndex} is not in range
         */
        public TreeNode getChild(int childIndex) {
            if (children == null || childIndex < 0 || children.length <= childIndex)
                throw new IndexOutOfBoundsException("childIndex is not in range [0, getChildCount()): childIndex="+childIndex+", getChildCount()="+getChildCount());
            return children[childIndex];
        }
        
        /**
         * Returns canonical number of the subtree rooted at this node as an
         * unsigned little-endian integer packed into a {@code long} array. If
         * this node's subtree was not in canonical form before this method is
         * called, it will be after. Use {@link #getCanonicalBitLength()} to
         * determine how many bits in the returned {@code long} array are being
         * used to represent the canonical number.
         *
         * @return the canonical number
         * @see #getCanonicalBitLength() 
         */
        public long[] getCanonicalNumber() {
            if (canonicalNumber == null) {
                // Basecase: 10
                if (children == null) {
                    canonicalNumber = new long[] { 0b10 };
                    canonicalBits = 2;
                } else {
                    // Recursivly generate the canonical number for every child
                    // node and total up the number of bits that will be need
                    // for this node's number
                    int totalBits = 2;
                    for (int c=0; c<children.length; c++) {
                        TreeNode child = children[c];
                        child.getCanonicalNumber();
                        totalBits += child.canonicalBits;
                    }
                    // Sort the children by canonical number in descending order
                    Arrays.sort(children, (TreeNode a, TreeNode b) -> {
                        if (a.canonicalBits < b.canonicalBits) return 1;
                        if (b.canonicalBits < a.canonicalBits) return -1;
                        long[] canonA = a.canonicalNumber;
                        long[] canonB = b.canonicalNumber;
                        for (int i=canonA.length-1; i>=0; i--) {
                            if (canonA[i] < canonB[i]) return 1;
                            if (canonB[i] < canonA[i]) return -1;
                        }
                        return 0;
                    });
                    // Create an array with enough room to hold this nodes's
                    // canonical number
                    long[] number = new long[(totalBits+63)/64];
                    // Keep track of the bit index. Note: The least significant 
                    // bit is always 0, and we start copying in bits at index 1
                    int bitIndex = 1; 
                    // For each of the children...
                    for (int c=children.length-1; c >= 0; c--) {
                        TreeNode child  = children[c];
                        // Get the child's canonical number
                        long[] childNumber = child.canonicalNumber;
                        int childBits = child.canonicalBits;
                        // Copy the child number into the current node's number
                        bitwiseCopy(childNumber, 0, number, bitIndex, childBits);
                        // Increment the bit index
                        bitIndex += childBits;
                    }
                    // Store a 1 at the most significant bit
                    number[bitIndex/64] |= (1L << (bitIndex%64));
                    // Store the canonical number
                    canonicalBits = bitIndex+1;
                    canonicalNumber = number;
                }
            }
            return canonicalNumber;
        }
        
        /**
         * Returns the number of bits being used to represent this subtree's
         * canonical number, or -1 if this node's subtree has not yet been
         * transformed into its canonical form.
         * 
         * @return the number of bits need to represent the the canonical
         * number, or {@code -1} if unknown
         * @see #getCanonicalNumber() 
         */
        public int getCanonicalBitLength() {
            return canonicalBits;
        }
        
        /**
         * Returns {@code true} if this node is a leaf node (i.e. it has no
         * children).
         * 
         * @return {@code true} if this node is a leaf node
         */
        public boolean isLeaf() {
            return children == null;
        }
        
        /**
         * Returns an iterator that traverses the nodes in the subtree rooted at
         * this node using pre-order traversal. The behavior of the returned
         * iterator is undefined if the tree is modified before the iteration
         * has completed (i.e. the iterator is not fast-fail and does not check
         * for concurrent modifications).
         * 
         * @return a pre-order traversal tree node iterator
         */
        public Iterator<TreeNode> iteratorPreOrder() {
            class IteratorPreOrder implements Iterator<TreeNode> {
                final Stack<TreeNode> stack = new Stack();

                IteratorPreOrder() {
                    stack.push(TreeNode.this);
                }

                @Override
                public boolean hasNext() {
                    return !stack.isEmpty();
                }

                @Override
                public TreeNode next() {
                    if (!hasNext())
                        throw new NoSuchElementException();
                    TreeNode node = stack.pop();
                    if (!node.isLeaf())
                        for (int i=node.getChildCount()-1; i>=0; i--)
                            stack.push(node.children[i]);
                    return node;
                }
            }
            return new IteratorPreOrder();
        }
        
        /**
         * Returns an iterator that traverses the nodes in the subtree rooted at
         * this node using pre-order traversal. The behavior of the returned
         * iterator is undefined if the tree is modified before the iteration
         * has completed (i.e. the iterator is not fast-fail and does not check
         * for concurrent modifications).
         * 
         * @return a post-order traversal tree node iterator
         */
        public Iterator<TreeNode> iteratorPostOrder() {
            class StackFrame {
                final TreeNode treeNode;
                int childIndex;

                StackFrame(TreeNode treeNode, int childIndex) {
                    this.treeNode = treeNode;
                    this.childIndex = childIndex;
                }
            }
            class IteratorPostOrder implements Iterator<TreeNode> {
                final Stack<StackFrame> stack = new Stack();

                IteratorPostOrder() {
                    stack.push(new StackFrame(TreeNode.this, 0));
                }

                @Override
                public boolean hasNext() {
                    return !stack.isEmpty();
                }

                @Override
                public TreeNode next() {
                    if (!hasNext())
                        throw new NoSuchElementException();
                    StackFrame frame = stack.peek();
                    TreeNode node = frame.treeNode;
                    while (frame.childIndex < node.getChildCount()) {
                        node = node.children[frame.childIndex++];
                        frame = new StackFrame(node, 0);
                        stack.push(frame);
                    }
                    stack.pop();
                    while (!stack.isEmpty()) {
                        frame = stack.peek();
                        if (frame.childIndex <= frame.treeNode.getChildCount())
                            break;
                        stack.pop();
                    }
                    return node;
                }
            }
            return new IteratorPostOrder();
        }
        
        /**
         * Returns an iterator that traverses the nodes in the subtree rooted at
         * this node using level-order (aka breadth-first search). The behavior
         * of the returned iterator is undefined if the tree is modified before
         * the iteration has completed (i.e. the iterator is not fast-fail and
         * does not check for concurrent modifications).
         *
         * @return a level-order traversal tree node iterator
         */
        public Iterator<TreeNode> iteratorLevelOrder() {
            class IteratorPreOrder implements Iterator<TreeNode> {
                final Queue<TreeNode> queue = new LinkedList();

                IteratorPreOrder() {
                    queue.add(TreeNode.this);
                }

                @Override
                public boolean hasNext() {
                    return !queue.isEmpty();
                }

                @Override
                public TreeNode next() {
                    if (!hasNext())
                        throw new NoSuchElementException();
                    TreeNode node = queue.remove();
                    TreeNode[] children = node.children;
                    if (children != null)
                        queue.addAll(Arrays.asList(children));
                    return node;
                }
            }
            return new IteratorPreOrder();
        }

        @Override
        public String toString() {
            StringBuilder str = new StringBuilder();
            str.append("(");
            str.append(vertex);
            if (children != null) {
                str.append(" ");
                for (TreeNode child: children) {
                    str.append(child);
                }
            }
            str.append(")");
            return str.toString();
        }
    }
    
    /**
     * Constructs and returns a new {@code Tree} generated by performing a
     * breadth-first search on the specified graph, beginning with the specified
     * node.
     * 
     * @param g graph that the tree will be on
     * @param rootIndex the index of the vertex from which the
     * breadth-first search will begin
     * @return a breadth-first search tree on the specified graph with the
     * specified root index
     * @throws NullPointerException if {@code g} is {@code null}
     * @throws IllegalArgumentException if {@code rootIndex} is not in range 
     * [0, {@code order})
     */
    public static Tree makeBfsTree(Graph g, int rootIndex) {
        int order = g.order;
        if (rootIndex < 0 || order <= rootIndex)
            throw new IllegalArgumentException("rootIndex is not in range [0,order): rootIndex="+rootIndex+", order="+order);
        
        // Get the ajacency list
        int[][] edges = g.getAdjacencyList();
        
        TreeNode[] queue = new TreeNode[order];
        int qidx = 0; // The index of the front of the queue
        int qend = 1; // The index beyond the last element in the queue
        int qchi = 1; // Elements between qchi (inclusive) and qend (exclusive)
                      // are children of the current node
        
        // Keepp track of which vertices have already been visited
        boolean[] marked = new boolean[order];
        
        // Add the root vertex to the queue and mark it as visted
        marked[rootIndex] = true;
        queue[0] = new TreeNode(rootIndex);
        
        // While the queue is not empty
        while (qidx < qend) {
            // Dequeue the next node
            TreeNode node = queue[qidx++];
            // Get the vertices adjacent to this node's vertex
            int[] adj = edges[node.vertex];
            // For each of the adjacenty vertices...
            for (int i=0; i<adj.length; i++) {
                int u = adj[i];
                // If the vertex has not been marked as visited....
                if (!marked[u]) {
                    // Mark the as visited
                    marked[u] = true;
                    // And enqueue a new node for the verted
                    queue[qend++] = new TreeNode(u);
                }
            }
            // If this node has children
            if (qchi < qend) {
                // Add the children to this node
                node.children = Arrays.copyOfRange(queue, qchi, qend);
                // Reposition reposition the child start index
                qchi = qend;
            }
        }
        
        // Create a list of all unique vertex indicies nodes in this tree
        int[] coverage = new int[qend];
        for (int i=0; i<coverage.length; i++)
            coverage[i] = queue[i].vertex;
        
        // Construct and return the tree
        return new Tree(g, queue[0], coverage, qend);
    }
    
    /**
     * Constructs and returns a new {@code Tree} generated by performing a
     * breadth-first search same-level repeats (BFS-SLR) search on the specified
     * graph, beginning with the specified root vertex.
     * <p>
     * A BFS-SLR tree is very similar to a BFS tree; it only differs from a BFS
     * tree in that a BFS-SLR tree allows the same graph vertex to appear more
     * than once in the tree, if depth at which the vertex the vertex is
     * encountered is equal, i.e. the repeated vertices appear on the same
     * level. Despite the fact that it allows vertices to be repeated, it is
     * still a tree since three no cycles between nodes in the tree (nodes wrap
     * vertices).
     * <p>
     * A BFS-SLR tree on a graph {@code G} rooted in vertex {@code v} contains
     * <i>every</i> shortest path from {@code v} to every other vertex in
     * {@code G} reachable from {@code v}.
     * <p>
     * A BFS-SLR, unlike a BFS tree, will always produce the same tree
     * regardless of the order that the vertices and edges are stored in memory.
     * <pre>{@code
     *    Graph G:         BFS-Tree(G, a):      BFS-SLR Tree(G, a):
     *                                                            
     *  a -- b -- c        a    (or)    a               a         
     *  |    |    |       / \          / \             / \        
     *  d -- e -- f      b   d        b   d           b   d       
     *                  / \           |   |          / \   \      
     *                 c   e          c   e         c   e   e     
     *                 |                  |         |   |   |     
     *                 f                  f         f   f   f     
     * }</pre>
     * 
     * @param g graph that the tree will be on
     * @param rootIndex the index of the vertex from which the
     * breadth-first search will begin
     * @return a breadth-first search tree on the specified graph with the
     * specified root index
     * @throws NullPointerException if {@code g} is {@code null}
     * @throws IllegalArgumentException if {@code rootIndex} is not in range 
     * [0, {@code order})
     */
    public static Tree makeBfsSlrTree(Graph g, int rootIndex) {
        int order = g.order;
        if (rootIndex < 0 || order <= rootIndex)
            throw new IllegalArgumentException("rootIndex is not in range [0,order): rootIndex="+rootIndex+", order="+order);
        
        // Get the ajacency list
        int[][] edges = g.getAdjacencyList();
        
        TreeNode[] queue = new TreeNode[order*order];
        int qidx = 0; // The index of the front of the queue
        int qend = 1; // The index beyond the last element in the queue
        int qlvl = 1; // The index beyond the last element of the current level
        int qchi = 1; // Elements between qchi (inclusive) and qend (exclusive)
                      // are children of the current node
        
        // Keepp track of which vertices have already been visited
        boolean[] marked = new boolean[order];
        
        // Keep track of the unique vertices
        int[] coverage = new int[order];
        int n = 0;
        
        // Add the root vertex to the queue and mark it as visted
        marked[rootIndex] = true;
        queue[0] = new TreeNode(rootIndex);
        
        // While the queue is not empty
        while (qidx < qend) {
            // For every node on the current level
            while (qidx < qlvl) {
                // Dequeue the next node
                TreeNode node = queue[qidx++];
                // Get the vertices adjacent to this node's vertex
                int[] adj = edges[node.vertex];
                // For each of the adjacenty vertices...
                for (int i=0; i<adj.length; i++) {
                    int u = adj[i];
                    // If the vertex has not been marked as visited....
                    if (!marked[u]) {
                        // And enqueue a new node for the vertex
                        queue[qend++] = new TreeNode(u);
                    }
                }
                // If this node has children
                if (qchi < qend) {
                    // Add the children to this node
                    node.children = Arrays.copyOfRange(queue, qchi, qend);
                    // Reposition reposition the child start index
                    qchi = qend;
                }
            }
            // Mark all of the verticies on the next level as visited
            while (qidx < qend) {
                int u = queue[qidx++].vertex;
                if (!marked[u]) {
                    marked[u] = true;
                    coverage[n++] = u;
                }
            }
            // Prepare the the next level
            qidx = qlvl;
            qlvl = qend;
        }
        
        // Construct and return the tree
        return new Tree(g, queue[0], coverage, qend);
    }
    
    /**
     * Copies bits from the source array to the destination array.
     * 
     * @param src the source array
     * @param srcPos the starting bit position in the source array
     * @param dst the destination array
     * @param dstPos the starting bit position in the destination array
     * @param length the number of bits to be copied
     * @throws NullPointerException if either {@code src} or {@code dst} is {@code null}
     */
    public static void bitwiseCopy(long[] src, int srcPos, long[] dst, int dstPos, int length) {
        if (length < 0)
            return;
        
        int srcI = srcPos / 64; // Source array index
        int srcB = srcPos % 64; // Source bit alignment
        int dstI = dstPos / 64; // Destination array index
        int dstB = dstPos % 64; // Destination bit aligntment

        int copyLength;
        long mask;
        
        if (srcB == dstB) {
            
            if (dstB != 0) {
                copyLength = Math.min(length, 64 - dstB);
                mask = ((1L << copyLength)-1) << dstB;
                dst[dstI] = (src[srcI++] & mask) | (dst[dstI] & ~mask);
                dstI++;
                length -= copyLength;
            }
                
            copyLength = length / 64;
            if (copyLength > 0) {
                System.arraycopy(src, srcI, dst, dstI, copyLength);
                srcI += copyLength;
                dstI += copyLength;
            }

            copyLength = length % 64;
            if (copyLength > 0) {
                mask = (1L << copyLength)-1;
                dst[dstI] = (src[srcI] & mask) | (dst[dstI] & ~mask);
            }
            
        } else {
            int off, ioff;
            
            if (srcB < dstB) {
                off = dstB - srcB;
                ioff = 64 - off;
                
                copyLength = Math.min(length, 64 - dstB);
                mask = ((1L << copyLength)-1) << dstB;
                dst[dstI] = ((src[srcI] << off) & mask) | (dst[dstI] & ~mask);
                dstI++;
                length -= copyLength;
            } else {
                off = srcB - dstB;
                ioff = 64 - off;
                
                copyLength = Math.min(length, 64 - srcB);
                mask = ((1L << copyLength)-1) << dstB;
                dst[dstI] = ((src[srcI++] >> off) & mask) | (dst[dstI] & ~mask);
                length -= copyLength;
                
                if (length > 0) {
                    copyLength = Math.min(length, off);
                    mask = ((1L << copyLength)-1) << ioff;
                    dst[dstI] = ((src[srcI] >> ioff) & mask) | (dst[dstI] & ~mask);
                    dstI++;
                    length -= copyLength;
                }
            }
            
            copyLength = length / 64;
            while (copyLength-- > 0) {
                dst[dstI] = (src[srcI++] >>> ioff) & (src[srcI] << off);
                dstI++;
            }
            length = length % 64;
            
            if (length > 0) {
                copyLength = Math.min(length, off);
                mask = ((1L << copyLength)-1);
                dst[dstI] = ((src[srcI++] >> ioff) & mask) | (dst[dstI] & ~mask);
                length -= copyLength;
            }
            
            if (length > 0) {
                mask = 1L << length;
                dst[dstI] = (src[srcI] & mask) | (dst[dstI] & ~mask);
            }
        }
    }
}
