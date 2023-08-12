import java.io.*;

public class CITS2200ProjectTester {
	public static void loadGraph(CITS2200Project project, String path) {
		// The graph is in the following format:
		// Every pair of consecutive lines represent a directed edge.
		// The edge goes from the URL in the first line to the URL in the second line.
		try {
			try (BufferedReader reader = new BufferedReader(new FileReader(path))) {
                while (reader.ready()) {
                	String from = reader.readLine();
                	String to = reader.readLine();
                	System.out.println("Adding edge from " + from + " to " + to);
                	project.addEdge(from, to);
                }
            }
		} catch (Exception e) {
			System.out.println("There was a problem:");
			System.out.println(e.toString());
		}
	}

	public static void main(String[] args) {
		// Change this to be the path to the graph file.
		String pathToGraphFile = "CITS2200/src/testdata.txt";
		
		CITS2200Project proj = new MyCITS2200Project();
		//loading the vertexes and edges
		loadGraph(proj, pathToGraphFile);
		

		//Test class for shortest path
		// int shortest = proj.getShortestPath("/wiki/Guinea-Bissau", "/wiki/Tajikistan");
		// System.out.println("Shortest Path is " + shortest);

		//Test class for hamiltonian path
		String[] hampath = proj.getHamiltonianPath();
		System.out.println(hampath.length);
		System.out.println("hampath is " + hampath[0]);

		//Test class for graph center
		// String[] center = proj.getCenters();
		// System.out.println("Center is " + center[0]);
		
		// //Test class for strongly connected components
		// String[][] scc = proj.getStronglyConnectedComponents();
		// System.out.println("Number of Strongest Connected Components is " + scc.length);
		
		}
	}