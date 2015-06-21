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

import java.util.Iterator;
import java.util.List;

/**
 *
 * @author Jeremiah N. Hankins
 */
public class Test {
    public static void main(String[] args) {
        final int minOrder = 2;
        final int maxOrder = 4;
        
        for (int order = minOrder; order <= maxOrder; order++) {
            Iterator<Graph> it = Graph.itLabeledGraphs(order);
            while (it.hasNext()) {
                Graph g = it.next();
                g.printLowerTriangle();
                System.out.println(g.toString());
                System.out.println();
            }
        }
        
//        final int minOrder = 2;
//        final int maxOrder = 6;
//        
//        for (int order = minOrder; order <= maxOrder; order++) {
//            Graph g = Graph.newStar(order, 0);
//            
//            List<Graph> iso = new ArrayList();
//            iso.add(new Graph(g));
//            g.printLowerTriangle();
//            
//            search: while (g.increment()) {
//                if (!g.isConnected())
//                    continue;
//                for (Graph h : iso)
//                    if (g.isIsomorphic(h))
//                        continue search;
////                if (g.getSize() < s-1)  {
////                    System.out.println(g.getSize()+" "+s);
////                    g.printLowerTriangle();
////                }
//                
//                iso.add(new Graph(g));
//            }
//            System.out.println("size: "+iso.size()+"\n");
//        }
    }
    
    public static List<Integer> makeList(int order) {
        return null;
    }
}
