package ru.lab;

public class Lab1 {
    private static double[] partialSums;

    private static double lnSeries(double x, int threads, int nTerms) throws InterruptedException {
        double z = (x - 1.0) / (x + 1.0);
        partialSums = new double[threads];
        int step = nTerms / threads;
        Thread[] workers = new Thread[threads];

        for (int t = 0; t < threads; t++) {
            final int threadIndex = t;
            final int start = t * step;
            final int end = (t == threads - 1) ? nTerms - 1 : (t + 1) * step - 1;

            workers[t] = new Thread(() -> {
                double localSum = 0.0;

                if (threadIndex == 0) {
                    try {
                        Thread.sleep(1); // 1 мс
                    } catch (InterruptedException ignored) {}
                }

                for (int n = start; n <= end; n++) {
                    double term = Math.pow(z, 2 * n + 1) / (2 * n + 1);
                    localSum += term;
                }
                partialSums[threadIndex] = localSum;
            });
            workers[t].start();
        }

        for (Thread t : workers) t.join();

        double total = 0.0;
        for (double s : partialSums) total += s;

        return 2 * total;
    }

    private static int estimateTerms(double x, double eps) {
        double z = (x - 1.0) / (x + 1.0);
        int n = 0;
        while (2 * Math.pow(Math.abs(z), 2 * n + 1) / (1 - z * z) > eps) n++;
        return n;
    }



    public void run() throws InterruptedException {
        double x = 0.0001;
        int[] threadCounts = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 36};
        double[] epsValues = {1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
        int repeats = 1000;

        warmUpJVM();

        for (double eps : epsValues) {
            System.out.printf("eps = %.0e%n", eps);
            for (int threads : threadCounts) {
                long totalTime = 0;
                int nTerms = estimateTerms(x, eps);
                double res = -999;
                for (int i = 0; i < repeats; i++) {
                    long start = System.nanoTime();
                    res = lnSeries(x, threads, nTerms);
                    long end = System.nanoTime();
                    totalTime += (end - start);
                }
                double avgTimeMs = totalTime / (repeats * 1_000_000.0);
                System.out.printf("%2d;%.4f;%3d;%.16f %n", threads, avgTimeMs, nTerms, res);
            }
            System.out.println();
        }
    }

    private static void warmUpJVM() throws InterruptedException {
        double x = 0.0001;
        double eps = 1e-6;
        int threads = 4;
        int nTerms = estimateTerms(x, eps);
        for (int i = 0; i < 500; i++) {
            lnSeries(x, threads, nTerms);
        }
        System.gc();
    }
}
