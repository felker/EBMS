public class CreateScatteringMatrix{
    public static void main(String[] args){
	int m = Integer.parseInt(args[0]);
	double matrix[][] = new double[m][m];

	//create pdf
	for (int i = 0; i < m; ++i)
	    for (int j = i; j < m; ++j)
		matrix[i][j] = 1./(m-i);
	//		System.out.println("PDF");
		//		System.out.printf("%d %d\n", m, m);

	//convert to cdf for mini-app
	for (int i = 0; i < m; ++i){
	    for (int j = 1; j < m; ++j){
		matrix[i][j] += matrix[i][j-1];
	    }
	    matrix[i][m-1] = 1.0;
	}

	System.out.printf("%d %d\n",m,m);
	for (int i = 0; i < m; ++i){
	    for (int j = 0; j < m; ++j){
		System.out.printf("%5.3f ",matrix[i][j]);
	    }
	    System.out.println();
	}

    }



}