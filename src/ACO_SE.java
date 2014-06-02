import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class ACO_SE {

	static int numberOfCells, numberOfMachines, numberOfParts,
			maximumNumberOfIteration;

	static double maximumEfficiency, initialThreshHold, visibilityRatio,
			evaporationRatio, minimumEfficiency;

	static int[][] partMachineMatrix;

	static double[][] commonPheromoneMatrix;

	static boolean DEBUG = false;

	static int[] bestTour = null;

	static {

		BufferedReader bin;
		try {
			bin = new BufferedReader(new FileReader("matrix.txt"));
			numberOfCells = Integer.parseInt(bin.readLine());
			numberOfMachines = Integer.parseInt(bin.readLine());
			numberOfParts = Integer.parseInt(bin.readLine());
			initialThreshHold = Double.parseDouble(bin.readLine());
			visibilityRatio = Double.parseDouble(bin.readLine());
			evaporationRatio = Double.parseDouble(bin.readLine());
			maximumNumberOfIteration = Integer.parseInt(bin.readLine());
			partMachineMatrix = new int[numberOfMachines][];
			for (int i1 = 0; i1 < numberOfMachines; i1++)
				partMachineMatrix[i1] = new int[numberOfParts];
			for (int i = 0; i < numberOfMachines; i++) {
				String[] row = bin.readLine().split(" ");
				int j = 0;
				for (String element : row) {
					partMachineMatrix[i][j] = Integer.parseInt(element);
					System.out.print(partMachineMatrix[i][j] + " ");
					j++;
				}
				System.out.println();
			}
			bin.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		commonPheromoneMatrix = new double[numberOfMachines + numberOfParts][];
		for (int i1 = 0; i1 < numberOfMachines + numberOfParts; i1++)
			commonPheromoneMatrix[i1] = new double[numberOfCells];

		for (int i = 0; i < numberOfMachines + numberOfParts; i++)
			for (int j = 0; j < numberOfCells; j++)
				commonPheromoneMatrix[i][j] = 1;
	}

	static void addToCommonExperience(double[][] pheromoneMatrix) {
		for (int i = 0; i < commonPheromoneMatrix.length; i++)
			for (int j = 0; j < commonPheromoneMatrix[i].length; j++)
				commonPheromoneMatrix[i][j] += pheromoneMatrix[i][j];
	}

	public static void main(String[] args) throws NumberFormatException,
			IOException {

		int numberOfAntGroups = Integer.parseInt(args[0]);
		int numberOfAntsInEachGroup = Integer.parseInt(args[1]);
		int numberOfIterationsForEachGroup = Integer.parseInt(args[2]);
		int numberOfAntUsingCommonMatrix = Integer.parseInt(args[3]);
		if (args != null && args.length > 0
				&& args[4].equalsIgnoreCase("debug"))
			DEBUG = true;

		String name = "sharedExp.txt";
		int fname = 0;
		while (new File(fname + "_" + name).exists())
			fname++;
		FileWriter fr = new FileWriter(fname + "_" + name);
		ArrayList<ACO> groups = new ArrayList<ACO>();
		for (int i = 0; i < numberOfAntGroups; i++) {
			ACO aco = new ACO(numberOfCells, numberOfMachines, numberOfParts,
					numberOfAntsInEachGroup, numberOfIterationsForEachGroup,
					initialThreshHold, visibilityRatio, evaporationRatio,
					partMachineMatrix, DEBUG);
			groups.add(aco);
		}
		for (int i = 0; i < maximumNumberOfIteration / numberOfAntGroups; i++) {
			System.out.println("PROGRESS ITERATION " + i);
			int j = 0;
			for (ACO aco : groups) {
				System.out.println("PROGRESS GROUP " + i + "," + j++);
				aco.start();
				addToCommonExperience(aco.pheromoneMatrix);
				if (maximumEfficiency < aco.maximumEfficiency
						&& minimumEfficiency < aco.minimumEfficiency) {
					maximumEfficiency = aco.maximumEfficiency;
					bestTour = aco.bestTour;
				}
				fr.write(aco.maximumEfficiency + "\t");
				for (int xi = 0; xi < numberOfMachines + numberOfParts; xi++)
					fr.write(aco.bestTour[xi] + " ");
				fr.write("\r\n");
			}
			{
				System.out.println("PROGRESS GROUP " + i + ",afterGroup");
				ACO aco = new ACO(numberOfCells, numberOfMachines,
						numberOfParts, numberOfAntUsingCommonMatrix,
						numberOfIterationsForEachGroup, initialThreshHold,
						visibilityRatio, evaporationRatio, partMachineMatrix,
						DEBUG);
				aco.updatePheromoneMatrix(commonPheromoneMatrix);
				aco.start();
				if (maximumEfficiency < aco.maximumEfficiency
						&& minimumEfficiency < aco.minimumEfficiency) {
					maximumEfficiency = aco.maximumEfficiency;
					bestTour = aco.bestTour;
					commonPheromoneMatrix = aco.pheromoneMatrix;
				}
				fr.write(aco.maximumEfficiency + "\t");
				for (int xi = 0; xi < numberOfMachines + numberOfParts; xi++)
					fr.write(aco.bestTour[xi] + " ");
				fr.write("\r\n");
			}
			for (ACO aco : groups) {
				aco.updatePheromoneMatrix(commonPheromoneMatrix);
			}
		}
		fr.flush();
		fr.close();
	}
}

class ACO {

	int numberOfCells, numberOfMachines, numberOfParts,
			maximumNumberOfIteration, numberOfAnts;
	double maximumEfficiency = -1, initialThreshHold, visibilityRatio,
			evaporationRatio, minimumEfficiency = -1;
	double[][] pheromoneMatrix;
	int[][] partMachineMatrix;
	int[] bestTour;
	boolean DEBUG = false;

	public void updatePheromoneMatrix(double[][] pherm) {
		for (int i = 0; i < pheromoneMatrix.length; i++)
			for (int j = 0; j < pheromoneMatrix[i].length; j++)
				pheromoneMatrix[i][j] += pherm[i][j];
	}

	ACO(int numberOfCells, int numberOfMachines, int numberOfParts,
			int numberOfAnts, int maximumNumberOfIteration,
			double initialThreshHold, double visibilityRatio,
			double evaporationRatio, int[][] partMachineMatrix, boolean DEBUG) {

		this.numberOfCells = numberOfCells;
		this.numberOfMachines = numberOfMachines;
		this.numberOfParts = numberOfParts;
		this.numberOfAnts = numberOfAnts;
		this.maximumNumberOfIteration = maximumNumberOfIteration;
		this.initialThreshHold = initialThreshHold;
		this.visibilityRatio = visibilityRatio;
		this.evaporationRatio = evaporationRatio;
		this.partMachineMatrix = partMachineMatrix;
		this.DEBUG = DEBUG;

		bestTour = new int[numberOfMachines + numberOfParts];

		for (int i1 = 0; i1 < numberOfMachines + numberOfParts; i1++)
			bestTour[i1] = 0;
		pheromoneMatrix = new double[numberOfMachines + numberOfParts][];
		for (int i1 = 0; i1 < numberOfMachines + numberOfParts; i1++)
			pheromoneMatrix[i1] = new double[numberOfCells];
		for (int i = 0; i < numberOfMachines + numberOfParts; i++)
			for (int j = 0; j < numberOfCells; j++)
				pheromoneMatrix[i][j] = 1;
	}

	void show_pm(int[] btour) {
		int[] mac = new int[numberOfMachines];
		int[] par = new int[numberOfParts];
		int pos = 0;
		for (int i = 0; i < numberOfCells; i++) {
			for (int j = 0; j < numberOfMachines; j++)
				if (btour[j] == i)
					mac[pos++] = j;
		}
		pos = 0;
		for (int i = 0; i < numberOfCells; i++) {
			for (int j = numberOfMachines; j < numberOfMachines + numberOfParts; j++)
				if (btour[j] == i)
					par[pos++] = j - numberOfMachines;
		}
		System.out.println("generated machine-part matrix : ");
		for (int i = 0; i < numberOfMachines; i++) {
			for (int j = 0; j < numberOfParts; j++)
				System.out.print(partMachineMatrix[mac[i]][par[j]] + " ");
			System.out.println();
		}
	}

	public void start() throws IOException {
		minimumEfficiency = initialThreshHold;
		maximumEfficiency = 0;
		for (int i = 0; i < maximumNumberOfIteration / numberOfAnts; i++) {
			int[] tour = null;
			double eff1 = -1;
			for (int j = 0; j < numberOfAnts; j++) {
				int[] tour_i = generate();
				double eff_i = evaluate(tour_i);
				if (eff_i > eff1) {
					tour = tour_i;
					eff1 = eff_i;
				}
			}

			if (eff1 > maximumEfficiency && eff1 > minimumEfficiency) {
				maximumEfficiency = (double) eff1;
				for (int xi = 0; xi < numberOfMachines + numberOfParts; xi++)
					bestTour[xi] = tour[xi];
				updatepher(tour, eff1);
			}
			if (DEBUG) {
				System.out.println("Tour Efficiency : " + eff1 + "\t");
				for (int xi = 0; xi < numberOfMachines + numberOfParts; xi++)
					System.out.print(tour[xi] + " ");
			}
			evaporate();
		}
		System.out.println("Best Tour Efficiency: " + maximumEfficiency + "\t");
		for (int xi = 0; xi < numberOfMachines + numberOfParts; xi++)
			System.out.print(bestTour[xi] + " ");
	}

	double evaluate(int[] tour) {
		int q1 = 0, n1 = 0, n2 = 0;
		for (int i = numberOfMachines; i < numberOfMachines + numberOfParts; i++) {
			for (int j = 0; j < numberOfMachines; j++) {
				if (tour[i] == tour[j]) {
					q1++;
					if (partMachineMatrix[j][i - numberOfMachines] == 1)
						n1++;
				} else if (partMachineMatrix[j][i - numberOfMachines] == 0)
					n2++;
			}
		}
		if (q1 == 0 || numberOfMachines * numberOfParts == q1)
			return 0;
		double qq = (double) q1 / (numberOfMachines * numberOfParts);
		double nn1 = (double) n1 / q1, nn2 = (double) n2
				/ (numberOfMachines * numberOfParts - q1);
		double eff = qq * nn1 + (1 - qq) * nn2;
		return eff;
	}

	void updatepher(int[] tour, double eff) throws IOException {
		for (int i = 0; i < numberOfMachines + numberOfParts; i++) {
			int tt = tour[i];
			pheromoneMatrix[i][tt] += eff;
		}
	}

	void evaporate() throws IOException {

		for (int i = 0; i < numberOfMachines + numberOfParts; i++) {
			for (int j = 0; j < numberOfCells; j++) {
				pheromoneMatrix[i][j] = pheromoneMatrix[i][j]
						* (1 - evaporationRatio);
			}
		}

		double sum = 0;
		for (int i = 0; i < numberOfMachines + numberOfParts; i++) {
			for (int j = 0; j < numberOfCells; j++) {
				sum += pheromoneMatrix[i][j];
			}
		}

		if (sum > 0) {
			for (int i = 0; i < numberOfMachines + numberOfParts; i++) {
				for (int j = 0; j < numberOfCells; j++) {
					pheromoneMatrix[i][j] = pheromoneMatrix[i][j] / sum;
				}
			}
		}
	}

	int[] generate() throws IOException {
		int[] tour = new int[numberOfMachines + numberOfParts];
		for (int i = 0; i < numberOfMachines + numberOfParts; i++) {
			double sum = 0;
			/********** normallizing *************/
			double[] pher = new double[numberOfCells];
			for (int j = 0; j < numberOfCells; j++)
				sum = sum + pheromoneMatrix[i][j];
			if (sum == 0)
				continue;
			for (int j = 0; j < numberOfCells; j++)
				pher[j] = pheromoneMatrix[i][j] / sum;
			int max = 0;
			for (int j = 1; j < numberOfCells; j++) {
				if (pher[max] < pher[j])
					max = j;
				if (pher[max] == pher[j]) {
					boolean ran1 = new Random().nextBoolean();
					if (ran1)
						max = j;
				}
			}
			int ranch = 0;
			for (int j = 1; j < numberOfCells; j++)
				pher[j] += pher[j - 1];

			double ran1 = new Random().nextDouble();

			if (ran1 < pher[0]) {
				ranch = 0;
			} else
				for (int j = 1; j < numberOfCells; j++) {
					if (ran1 >= pher[j - 1] && ran1 < pher[j]) {
						ranch = j;
						break;
					}
				}
			double ran = new Random().nextDouble();
			if (ran > visibilityRatio)
				max = ranch;
			tour[i] = max;
		}
		return tour;
	}

}