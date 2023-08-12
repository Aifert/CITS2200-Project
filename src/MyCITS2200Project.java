
import java.util.*;

public class MyCITS2200Project implements CITS2200Project {
    

    /**
     * CODE FOR :
     * Aifert Yet 23455873
     * Daniel Le 23625105
     * 
     * In this project, our vertices are represented by strings, the URL of the Wikipedia page.
     * Indexing by these strings can be messy and inefficient, so we instead assign each vertex a
     * unique integer ID between in the range [0, V ) . This ID will serve as an index into the
     * adjacency list, allowing us to find a vertex's list of neighbours in constant time. 
     * 
     * To allow us to
     * convert back and forth between the string and integer representations of our vertices,
     * we introduce a list of strings that can be indexed efficiently by vertex ID, and a map from the
     * vertex URL to its ID.
     */
    
    /**
     * SCC = Strongly Connected Component
     * To solve strongly connected components in Q3, we can use Kosaraju's Algorithm
     * it uses the original list, a transpose of the original list and a stack to find SCCs
     * Highlight : "a transpose of the original list" means the original list with all the edges
     *             reversed
     * https://www.youtube.com/watch?v=5wFyZJ8yH9Q reference video
     */

    /**
     * ArrayList to allow us to find the url(string) with an pre-assigned id(int)
     * be able to access them in constant time
     */
     private final ArrayList<String> idURL = new ArrayList<String>();

    /**
     * allows us to find the id(int) with url (string)
     * using linkedhashmap instead of hashmap because since linkedhashmap
     * is ordered it may make it easier to store and access different 
     * information needed
    */
    private final LinkedHashMap<String,Integer> URLid = new LinkedHashMap<>();

    /**
     * Original list for the graph
     */

    private final ArrayList<List<Integer>> originalList = new ArrayList<>();

    /**
     * Transposed list for the graph
     */

     private final ArrayList<List<Integer>> transposedList = new ArrayList<>();
    
    /**
     * Field for Hamiltonian Path
     */

    private int N;
    private int start;
    private int end;

    private double[][] adjacencyMatrix;

    private List<Integer> path = new ArrayList<>();
    private boolean ranSolver = false;

    /**
     * Now we have to fill the list with url id (int)
     * 
     * @param urlStart start of the edge
     * @param urlEnd end of the edge
     */
    @Override

    public void addEdge(String urlStart, String urlEnd){
        //make the url into vertices for our graph
        addVertex(urlStart);
        addVertex(urlEnd);

        //add edges
        int start = URLid.get(urlStart);
        int end = URLid.get(urlEnd);

        originalList.get(start).add(end);

        //transposed is just the start of original == end of transposed
        //vice versa
        transposedList.get(end).add(start);
    }

    public void addVertex(String url){
        if(!URLid.containsKey(url)){
            //add the url into our list
            idURL.add(url);
            //as this is a new entry the index of it will just be the size, since 
            //indexing in hashmaps start from 0 and size == index + 1
            //which means an id for the new entry
            URLid.put(url, URLid.size());

            //add into list
            originalList.add(new ArrayList<>());
            transposedList.add(new ArrayList<>());
        }
    }

    //Question 1 : Shortest Path
    //Write a method that, given a pair of pages, returns the minimum number of links you must follow to get from the first page to the second

    /**
     * BFS = Breath First Search
     * Assuming our graph is an unweighted directed graph, we can perform BFS to find the minimum number of edges between
     * start vertex and end vertex
     * that will give us the shortest path between the two given vertices
     */

    /**
     * After performing BFS on our graph, we will return an array which stores the length of shortest path 
     * from our source to each vertex it reaches
     * if there is no path it will take the original value of array
     * Time Complexity : O(V+E) due to BFS, but assuming BFS is not taken into account it is O(1)
     * 
     * @param urlStart
     * @param urlEnd
     * @return distance
     */
    @Override

    public int getShortestPath(String urlStart, String urlEnd){
        int [] result = BreathFirstSearch(URLid.get(urlStart));

        //return shortest path of the end index from start
        return result[URLid.get(urlEnd)];
    }

    /**
     * Finding the lengts of shortest paths by performing a BFS
     * which records vertices according to the number of edges from our start
     * we create an empty list to fill in with the distances of each vertex from our starting vertex
     * Time Complexity : O(V+E)
     * @param start
     * @return distances stored in array
     */
    private int[] BreathFirstSearch(int start){
        //Distance from start to each vertex in graph
        //using a linkedlist structure to keep every entry ordered
        Queue<Integer> queue = new LinkedList<>();
        int[] visited = new int[idURL.size()];

        //to show every vertex is unvisited at the start,
        //mark all entries in visited as -1
        Arrays.fill(visited, -1);

        //since we visit the start, we will mark it as 0
        //visited and add it into our queue
        visited[start] = 0;
        queue.add(start);

        //performing BFS
        while(!queue.isEmpty()){
            //remove start of this queue and store it 
            //for comparisons
            int curr = queue.remove();
            // retrieves the list of integers associated with curr integer key in map
            for(int next : originalList.get(curr)){
                if(visited[next] == -1){
                    //add distance to array
                    // +1 to mark as visited 
                    visited[next] = visited[curr] + 1;
                    //add vertex in our queue
                    queue.add(next);
                }
            }
        }  
        return visited;
    }   

    //Question 2 : Hamiltonian Path / Hamiltonian Cycle
    //Write a method that finds a Hamiltonian path in a Wikipedia page graph. 
    //A Hamiltonian path is any path in some graph that visits every vertex exactly once. 
    //This method will never be called for graphs with more than 20 pages. See here for more information about Hamiltonian paths.

    /*
     * Using Math.pow(2, n) instead of (1 << n): Java's Math.pow() function is not constant
     * time, but is rather logarithmic in the power. This introduced an O(log n) factor that is not
     * present when using bit-shifts.
     * 
     * Best Solution to use is dynamic programming known as Held-Karp or Bellman-Held-Karp algorithm.
     * We will operate directly on bitmasks. Which masks the value of an original 
     * bit with a mask, revealing only the bit we want
     * How bitmasks works is, if we apply the mask and it is still 1, we have visited, or else we have not visited that vertex
     * and can explore it 
     * https://www.youtube.com/watch?v=0j0TF_u4FYw 
     * https://www.youtube.com/watch?v=cY4HiiFHO1o
     */

     public String[] getHamiltonianPath(){
        //Setting up environment and adjacencyMatrix for the algorithm
        adjacencyMatrix = setUp(originalList);

        //Retrieve path
        List<Integer> path = getPath();

        //Construct path with id -> URL
        String[] hampath = new String[path.size()];
        for (int i = 0; i < hampath.length; i++){
            hampath[i] = idURL.get(path.get(i));
        }
        return hampath;
    }
    /**
     * Creating an adjacency Matrix with 1 being having an edge and 0 being none
     * Also initilalises field variables
     * @param originalList
     * @return
     */
    public double[][] setUp(ArrayList<List<Integer>> originalList){
        int graphSize = originalList.size();
        adjacencyMatrix = new double[graphSize][graphSize];

        //for each graph index
        for(List<Integer> set : originalList){
            //find the from and to
            int idFrom = originalList.indexOf(set);
            if(set.size() != 0){
                for(int i = 0; i < set.size(); i++){
                    int idTo = set.get(i);
                    //mark the vertex as 1 if they have an edge
                    //0 otherwise
                    adjacencyMatrix[idFrom][idTo] = 1.0;
                }
            }
        }   
        //Initialising field variables
        start = 0;
        N = adjacencyMatrix.length;
        end = (1 << N) - 1;
        return adjacencyMatrix;
    }

    public List<Integer> getPath() {
        if (!ranSolver) solve();
        return path;
      }
    
    public void solve() {
    
        // Run the solver
        int state = 1 << start;

        //Creating a memo to remember which vertex we have visited already
        Double[][] memo = new Double[N][1 << N];

        //Creating a previous to back track and construct hamiltonian path
        Integer[][] prev = new Integer[N][1 << N];
        iterateAdjacencyMatrix(start, state, memo, prev);
    
        // Regenerate path
        int index = start;
        while (true) {
          path.add(index);
          Integer nextIndex = prev[index][state];

          //no edge, move onto next
          if (nextIndex == null) break;

          //if after applying mask is still 0, add it into our collection of visited vertexes
          int nextState = state | (1 << nextIndex);
          state = nextState;
          index = nextIndex;
        }
        path.add(start);
        ranSolver = true;
      }
      /**
       * Traverses through the adjacency matrix and find where there is a 1 (edge) between edges in the ijth entry of matrix
       * Construct a hamiltonian path consisting of index of vertices
       * @param i
       * @param state
       * @param memo
       * @param prev
       * @return
       */
      private double iterateAdjacencyMatrix(int i, int state, Double[][] memo, Integer[][] prev) {
    
        // Done this tour. Return cost of going back to start node.
        if (state == end) {
            return adjacencyMatrix[i][start];
        }
        // Return cached answer if already computed.
        if (memo[i][state] != null) {
            return memo[i][state];
        }
    
        double minCost = Double.POSITIVE_INFINITY;
        int index = -1;
        for (int next = 0; next < N; next++) {
    
          // Skip if the next node has already been visited.
          if ((state & (1 << next)) != 0) continue;
    
          int nextState = state | (1 << next);
          double newCost = adjacencyMatrix[i][next] + iterateAdjacencyMatrix(next, nextState, memo, prev);
          if (newCost < minCost) {
            minCost = newCost;
            index = next;
          }
          //break we have found path.
          break;
        }
    
        prev[i][state] = index;
        return memo[i][state] = minCost;
      }

    //Question 3: Strongly Connected Components
    //Write a method that finds every ‘strongly connected component’ of pages.
    // Korsaraju's Algorithm is helpful
    //https://www.youtube.com/watch?v=V8qIqJxCioo - Kosaraju's Algorithm for Strongly Connected Components (SCC)

    /**
     * Returns set of vertices that can reach each other, if it terminates, that is ONE SCC
     * 
     * @return The Strong Connected Component
     */
    @Override
    public String[][] getStronglyConnectedComponents(){
        //storing results
        ArrayList<Stack<Integer>> result = new ArrayList<>();

        //Kosaraju's algorithm to find strongly connected components
        kosarajuAlgorithm(result);

        //2d array to store strongly connected components and the vertices
        String[][] strong_components = new String[result.size()][]; 
        
        //result
        return kosarajuResult(result,strong_components);
    }

    /**
     * Format Kosaraju's Algorithm
     * 
     * divide vertices of graph to SCC
     * 
     * @param result ArrayList to store result
     * @param components 2D Array to store formatted result
     * @return the SCC
     */

    private String[][] kosarajuResult(ArrayList<Stack<Integer>> result, String[][] strong_components){
        //add components
        for(int unique = 0; unique < result.size(); unique++){
            //find number
            strong_components[unique] = new String[result.get(unique).size()];

            //set the size of array
            int size = result.get(unique).size();
            
            //Add SCC
            for (int connect = 0; connect < size; connect++){
                //find url of vertices
                strong_components[unique][connect] = idURL.get(result.get(unique).get(connect));
            }
        }
        return strong_components; 
    }

    /**
     * Performing Kosaraju's algorithm
     * 
     * Compute the set of all vertices that can reach other by performing a DFS (Depth First Search), starting at any vertex
     * Doing this with the transpose of the graph is also necessary to find which vertexes actually "CONNECT" to each other, can go out and come back
     * Kosaraju's Algorithm uses DFS order through the original graph to ensure DFS through tranpose graph only explores said intersaction
     * 
     * @param result
     */

    private void kosarajuAlgorithm(ArrayList<Stack<Integer>> result){
        //mark all as non visited
        boolean [] visited = new boolean[idURL.size()];

        //create stack to perform algorithm
        Stack<Integer> stack = new Stack<>();

        for(int i = 0; i < idURL.size(); i++){
            if(!visited[i]){
                depthFirstSearch(i,visited, true, stack);
            }
        }
        //mark all as non visited(tranpose of original graph)
        visited = new boolean[idURL.size()];

        while(!stack.isEmpty()){
            int curr = stack.pop();
            if(!visited[curr]){
                Stack<Integer> component = new Stack<>();
                depthFirstSearch(curr, visited, false, component);
                result.add(component);
            }
        }
    }

    /**
     * Performing DFS
     * DFS through the transpose graph starting from the last vertex and any vertices it visits must be an SCC
     * Because for it to be in the stack, it must have connected from a previous vertex, if the tranpose connects it back to the previous vertex, it is a SCC
     * DFS again from the next highest vertex in the stack that is not yet visited to find all SCC
     * repeat until found all SCC
     * 
     * @param curr      Starting Position
     * @param visited   boolean array for visits
     * @param original  original or tranposed list
     * @param stack     stack to store results
     */

     private void depthFirstSearch(int curr, boolean[] visited, boolean original, Stack<Integer> stack){
        //mark curr as visited
        visited[curr] = true;

        //to use original list or transposed list
        ArrayList<List<Integer>> currList = original ? originalList : transposedList;
        
        //DFS
        for (int next : currList.get(curr)){
            if(!visited[next]){
                depthFirstSearch(next, visited, original, stack);
            }
        }
        stack.add(curr);
     }
     
    //Question 4 Finding Graph Centers
     //Write a method that finds all the centers of the Wikipedia page graph. 
     //A vertex is considered to be the center of a graph if the maximum shortest path from 
     //that vertex to any other vertex is the minimum possible


    /**
     * BFS will be good for this as it finds the neighbours of each vertex 
     * and consider it for a member of center and returning the length of edge
     * 
     * @return the center
     */
    @Override
    public String[] getCenters(){
        return center().toArray(new String[0]);
    }

    /**
     * The center (or Jordan center[1]) of a graph is the set of all vertices of minimum eccentricity,[2] - Wikipedia (explanation source given to us)
     * 
     * One way to do it is using all pairs shortest path and finding which one has the smallest maximum distance
     * @return List of center
     */
    private ArrayList<String> center(){
        ArrayList<String> result = new ArrayList<>();

        int min = idURL.size();
        int url = 0;
        int urlSize = idURL.size();

        //for each url
        while(url < urlSize){
            int eccentricity = -1;
            //BFS to find the shortest path to each vertex
            for (int vertex : BreathFirstSearch(url)){
                //find the further vertex from url
                if(eccentricity < vertex){
                    //BFS returns -1 if can't reach, if it is not then assign eccentricity to vertex
                    if(vertex > -1){
                        eccentricity = vertex;
                    }
                }
            }
        //vertices with min eccentricity
        if(eccentricity < min){
            //check if its smaller because of -1 if not then update minimum
            if(eccentricity > -1){
                min = eccentricity;
                //since new min need to update the result 
                result = new ArrayList<>();
            }
        }
        if(eccentricity == min){
            result.add(idURL.get(url));
        }
        url++;
        }
        return result;
    }

}
