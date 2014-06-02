import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

public class ACO_Tabu {
	int k, m, n, q;
	double nmax, t, v, y, nmin;
	double[][] km;
	int[][] pm;
	int[] btour;
	static boolean DEBUG = false;

	FileWriter fr_log = new FileWriter("tabu_log.txt");

	int[] generate() throws IOException {
		int[] tour = new int[m + n];
		for (int i = 0; i < m + n; i++) {
			double sum = 0;
			/********** normallizing *************/
			double[] pher = new double[k];
			for (int j = 0; j < k; j++)
				sum = sum + km[i][j];
			if (sum == 0)
				continue;
			for (int j = 0; j < k; j++)
				pher[j] = km[i][j] / sum;
			int max = 0;
			for (int j = 1; j < k; j++) {
				if (pher[max] < pher[j])
					max = j;
				if (pher[max] == pher[j]) {
					boolean ran1 = new Random().nextBoolean();
					if (ran1)
						max = j;
				}
			}
			int ranch = 0;
			for (int j = 1; j < k; j++)
				pher[j] += pher[j - 1];

			if (DEBUG) {
				fr_log.write("\n");
				fr_log.write("Cumulative Normallized Pher Matrix Row");
				fr_log.write("\n");
				for (int j = 0; j < k; j++) {
					fr_log.write(pher[j] + " ");
				}
				fr_log.write("\n");
				fr_log.write("********************************");
				fr_log.flush();
			}

			double ran1 = new Random().nextDouble();

			if (ran1 < pher[0]) {
				ranch = 0;
			} else
				for (int j = 1; j < k; j++) {
					if (ran1 >= pher[j - 1] && ran1 < pher[j]) {
						ranch = j;
						break;
					}
				}
			double ran = new Random().nextDouble();
			if (ran > v)
				max = ranch;
			if (DEBUG) {
				fr_log.write("\n");
				fr_log.write("Chosen Value : " + max + " [" + ran1 + ","
						+ ranch + "," + ran + "]");
				fr_log.write("\n");
				fr_log.write("********************************");
				fr_log.flush();
			}
			tour[i] = max;
		}

		double efft1 = evaluate(tour);

		double eff2 = -1;
		if (DEBUG) {
			fr_log.write("\n");
			fr_log.write("Efficiency before tabu " + efft1);
			fr_log.write("\n");
			fr_log.flush();
		}
		int[][] tabu = new int[k][];
		for (int u2 = 0; u2 < k; u2++)
			tabu[u2] = new int[2];
		int fill = 0;
		int[] toure = new int[m + n];
		for (int u2 = 0; u2 < m + n; u2++)
			toure[u2] = tour[u2];
		for (int u = 0; u < k; u++) {
			int ran1 = new Random().nextInt(m + n);
			int ran2 = new Random().nextInt(k);
			int flag = 0;
			for (int u1 = 0; u1 < fill; u1++) {
				if (tabu[u1][0] == ran1 && tabu[u1][1] == ran2) {
					flag = 1;
					break;
				}
			}
			if (flag == 1) {
				u--;
				continue;
			}
			tabu[fill][0] = ran1;
			tabu[fill][1] = ran2;
			fill++;
			int[] tour1 = new int[m + n];
			for (int u1 = 0; u1 < m + n; u1++)
				tour1[u1] = tour[u1];
			tour1[ran1] = ran2;

			eff2 = evaluate(tour1);

			if (DEBUG) {
				fr_log.write("\n");
				fr_log.write("Efficiency after tabu " + eff2);
				fr_log.write("\n");
				fr_log.flush();
			}

			if (efft1 < eff2) {
				for (int u1 = 0; u1 < m + n; u1++)
					toure[u1] = tour1[u1];
				efft1 = evaluate(toure);
			}
		}
		return toure;
	}

	void useTabuOnBestResult() throws IOException {

		int[][] tabu = new int[k][];
		for (int u2 = 0; u2 < k; u2++)
			tabu[u2] = new int[2];
		int fill = 0;
		int[] toure = new int[m + n];
		for (int u2 = 0; u2 < m + n; u2++)
			toure[u2] = btour[u2];
		for (int u = 0; u < k; u++) {
			int ran1 = new Random().nextInt(m + n);
			int ran2 = new Random().nextInt(k);
			int flag = 0;
			for (int u1 = 0; u1 < fill; u1++) {
				if (tabu[u1][0] == ran1 && tabu[u1][1] == ran2) {
					flag = 1;
					break;
				}
			}
			if (flag == 1) {
				u--;
				continue;
			}
			tabu[fill][0] = ran1;
			tabu[fill][1] = ran2;
			fill++;
			int[] tour1 = new int[m + n];
			for (int u1 = 0; u1 < m + n; u1++)
				tour1[u1] = btour[u1];
			tour1[ran1] = ran2;

			double eff2 = evaluate(tour1);

			if (DEBUG) {
				fr_log.write("\n");
				fr_log.write("Efficiency after tabu only on best" + eff2);
				fr_log.write("\n");
				fr_log.flush();
			}

			if (eff2 > nmax && eff2 > nmin) {
				for (int u1 = 0; u1 < m + n; u1++)
					btour[u1] = tour1[u1];
				nmax = eff2;
				updatepher(btour, eff2);
			}
		}
	}

	double evaluate(int[] tour) {
		int q1 = 0, n1 = 0, n2 = 0;
		for (int i = m; i < m + n; i++) {
			for (int j = 0; j < m; j++) {
				if (tour[i] == tour[j]) {
					q1++;
					if (pm[j][i - m] == 1)
						n1++;
				} else if (pm[j][i - m] == 0)
					n2++;
			}
		}
		if (q1 == 0 || m * n == q1)
			return 0;
		double qq = (double) q1 / (m * n);
		double nn1 = (double) n1 / q1, nn2 = (double) n2 / (m * n - q1);
		double eff = qq * nn1 + (1 - qq) * nn2;
		return eff;
	}

	void updatepher(int[] tour, double eff) throws IOException {
		for (int i = 0; i < m + n; i++) {
			int tt = tour[i];
			km[i][tt] += eff;
		}

		if (DEBUG) {
			fr_log.write("\n");
			fr_log.write("Update Pher");
			fr_log.write("\n");
			for (int i = 0; i < m + n; i++) {
				for (int j = 0; j < k; j++) {
					fr_log.write(km[i][j] + " ");
				}
				fr_log.write("\n");
			}
			fr_log.write("********************************");
			fr_log.flush();
		}
	}

	void evaporate() throws IOException {

		for (int i = 0; i < m + n; i++) {
			for (int j = 0; j < k; j++) {
				km[i][j] = km[i][j] * (1 - y);
			}
		}

		double sum = 0;
		for (int i = 0; i < m + n; i++) {
			for (int j = 0; j < k; j++) {
				sum += km[i][j];
			}
		}

		if (sum > 0) {
			for (int i = 0; i < m + n; i++) {
				for (int j = 0; j < k; j++) {
					km[i][j] = km[i][j] / sum;
				}
			}
		}

		if (DEBUG) {
			fr_log.write("\n");
			fr_log.write("Evaporate");
			fr_log.write("\n");
			for (int i = 0; i < m + n; i++) {
				for (int j = 0; j < k; j++) {
					fr_log.write(km[i][j] + " ");
				}
				fr_log.write("\n");
			}
			fr_log.write("********************************");
			fr_log.flush();
		}
	}

	void start() throws IOException {
		String name = "tabu.txt";
		int fname = 0;
		while (new File(fname + "_" + name).exists())
			fname++;
		FileWriter fr = new FileWriter(fname + "_" + name);
		nmin = t;
		nmax = 0;
		for (int i = 0; i < q; i++) {
			int[] tour = generate();
			double eff1 = evaluate(tour);

			if (DEBUG) {
				fr_log.write("\n");
				fr_log.write("Got tabu eff " + eff1);
				fr_log.write("Old eff " + nmax);
				fr_log.write("\n");
				fr_log.flush();
			}
			fr.write(eff1 + "\t");
			for (int xi = 0; xi < m + n; xi++)
				fr.write(btour[xi] + " ");
			fr.write("\n");
			if (eff1 > nmax && eff1 > nmin) {
				nmax = (double) eff1;
				for (int xi = 0; xi < m + n; xi++)
					btour[xi] = tour[xi];
				updatepher(tour, eff1);
			} else {
				useTabuOnBestResult();
			}
			evaporate();
		}
		System.out.println("best tour");
		for (int xi = 0; xi < m + n; xi++)
			System.out.print(btour[xi] + " ");
		System.out.println("max eff= : " + nmax);
		show_pm(btour);
		fr.write(nmax + "");
		fr.flush();
		fr.close();
	}

	void show_pm(int[] btour) {
		int[] mac = new int[m];
		int[] par = new int[n];
		int pos = 0;
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < m; j++)
				if (btour[j] == i)
					mac[pos++] = j;
		}
		pos = 0;
		for (int i = 0; i < k; i++) {
			for (int j = m; j < m + n; j++)
				if (btour[j] == i)
					par[pos++] = j - m;
		}
		System.out.println("generated machine-part matrix : ");
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++)
				System.out.print(pm[mac[i]][par[j]] + " ");
			System.out.println();
		}
	}

	public ACO_Tabu() throws NumberFormatException, IOException {

		BufferedReader bin = new BufferedReader(new FileReader("matrix.txt"));
		k = Integer.parseInt(bin.readLine());
		m = Integer.parseInt(bin.readLine());
		n = Integer.parseInt(bin.readLine());
		t = Double.parseDouble(bin.readLine());
		v = Double.parseDouble(bin.readLine());
		y = Double.parseDouble(bin.readLine());
		q = Integer.parseInt(bin.readLine());
		pm = new int[m][];
		for (int i1 = 0; i1 < m; i1++)
			pm[i1] = new int[n];
		for (int i = 0; i < m; i++) {
			String[] row = bin.readLine().split(" ");
			int j = 0;
			for (String element : row) {
				pm[i][j] = Integer.parseInt(element);
				System.out.print(pm[i][j] + " ");
				j++;
			}
			System.out.println();
		}
		btour = new int[m + n];
		bin.close();
		km = new double[m + n][];
		for (int i1 = 0; i1 < m + n; i1++)
			km[i1] = new double[k];
		for (int i1 = 0; i1 < m + n; i1++)
			btour[i1] = 0;
		for (int i = 0; i < m + n; i++)
			for (int j = 0; j < k; j++)
				km[i][j] = 1;
	}

	public static void main(String[] args) throws NumberFormatException,
			IOException {
		if (args != null && args.length > 0
				&& args[0].equalsIgnoreCase("debug"))
			DEBUG = true;
		ACO_Tabu ant = new ACO_Tabu();
		ant.start();
		ant.fr_log.close();
	}
}