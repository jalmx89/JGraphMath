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

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Jeremiah N. Hankins
 */
public class Test {
    public static void main(String[] args) {
        
        final int maxOrder = 8;
        
        for (int order = 2; order <= maxOrder; order++) {
            Graph g = new Graph(order);
            List<Graph> iso = new ArrayList();
            iso.add(new Graph(g));
            
            incr: while (g.increment()) {
                for (int i=0; i<iso.size(); i++) {
                    Graph h = iso.get(i);
                    if (h.isIsomorphicMapJohnsonTrotter(g)) {
                        continue incr;
                    }
                }
                iso.add(new Graph(g));
            }
            
            System.out.println(iso.size()+" classes of isomorphic graphs with "+order+" vertexes");
        }
    }
}
