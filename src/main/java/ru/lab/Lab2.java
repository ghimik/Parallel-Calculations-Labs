package ru.lab;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.BiFunction;
import java.util.function.Function;

public class Lab2 {

    public static class Integrator {

        // левые прямоугольники
        public static double rectangle(Function<Double, Double> f, double a, double b, int n) {
            double h = (b - a) / n;
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                double x = a + i * h;
                sum += f.apply(x);
            }
            return sum * h;
        }

        public static double trapezoid(Function<Double, Double> f, double a, double b, int n) {
            double h = (b - a) / n;
            double sum = 0.5 * (f.apply(a) + f.apply(b));
            for (int i = 1; i < n; i++) {
                double x = a + i * h;
                sum += f.apply(x);
            }
            return sum * h;
        }

        public static double simpson(Function<Double, Double> f, double a, double b, int n) {
            if (n % 2 != 0) n++;
            double h = (b - a) / n;
            double sum = f.apply(a) + f.apply(b);

            for (int i = 1; i < n; i++) {
                double x = a + i * h;
                sum += f.apply(x) * (i % 2 == 0 ? 2 : 4);
            }
            return sum * h / 3.0;
        }
    }

    public static class ParallelIntegrator {

        private static double parallelSum( double a, double b, int n,
                                          int threads, BiFunction<Double, Double, Double> methodChunk) {
            ExecutorService pool = Executors.newFixedThreadPool(threads);
            List<Future<Double>> futures = new ArrayList<>();
            int chunk = n / threads;
            double h = (b - a) / n;

            for (int t = 0; t < threads; t++) {
                int start = t * chunk;
                int end = (t == threads - 1) ? n : (t + 1) * chunk;
                futures.add(pool.submit(() -> methodChunk.apply((double) start, (double) end)));
            }

            double sum = 0;
            try {
                for (Future<Double> fPart : futures)
                    sum += fPart.get();
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                pool.shutdown();
            }
            return sum * h;
        }

        public static double rectangle(Function<Double, Double> f, double a, double b, int n, int threads) {
            double h = (b - a) / n;
            return parallelSum(a, b, n, threads, (start, end) -> {
                double local = 0;
                for (int i = start.intValue(); i < end.intValue(); i++) {
                    double x = a + i * h;
                    local += f.apply(x);
                }
                return local;
            });
        }

        public static double trapezoid(Function<Double, Double> f, double a, double b, int n, int threads) {
            double h = (b - a) / n;
            double edge = 0.5 * (f.apply(a) + f.apply(b));
            double body = parallelSum(a, b, n, threads, (start, end) -> {
                double local = 0;
                for (int i = start.intValue() + 1; i < end.intValue(); i++) {
                    double x = a + i * h;
                    local += f.apply(x);
                }
                return local;
            });
            return h * (edge + body);
        }

        public static double simpson(Function<Double, Double> f, double a, double b, int n, int threads) {
            if (n % 2 != 0) n++;
            double h = (b - a) / n;
            double edge = f.apply(a) + f.apply(b);

            double body = parallelSum(a, b, n, threads, (start, end) -> {
                double local = 0;
                for (int i = start.intValue() + 1; i < end.intValue(); i++) {
                    double x = a + i * h;
                    local += f.apply(x) * (i % 2 == 0 ? 2 : 4);
                }
                return local;
            });
            return (edge + body) * h / 3.0;
        }
    }

    public static class IntegrationUtils {

        public static Function<Double, Double> funcA() {
            return x -> x / Math.sqrt(x * x - 16);
        }

        public static Function<Double, Double> funcB() {
            return x -> {
                double ln = Math.log(x);
                return Math.pow(x, 4) * ln * ln;
            };
        }

        public static Function<Double, Double> funcC() {
            return x -> 1.0 / (3 * Math.pow(Math.cos(x), 2) + 2 * Math.pow(Math.sin(x), 2));
        }

        public static Function<Double, Double> funcD() {
            return x -> Math.sqrt(x) / (Math.pow(Math.pow(x, 3), 0.25) + 1);
        }

        public static Function<Double, Double> funcE() {
            return x -> (Math.pow(x, 3) - 2 * Math.pow(x, 2) - 12 * x - 7)
                    / (Math.pow(x, 3) - Math.pow(x, 2));
        }

    }

    public static void main(String[] args) {
        final double[] ACCURACIES = {1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13};
        final int[] THREAD_COUNTS = {2, 4};
        final int REPEATS = 100;

        List<Function<Double, Double>> functions = List.of(
                IntegrationUtils.funcA(), // x / sqrt(x² - 16)
                IntegrationUtils.funcB(), // x⁴ * (ln x)²
                IntegrationUtils.funcC(), // 1 / (3cos²x + 2sin²x)
                IntegrationUtils.funcD(), // sqrt(x) / ⁴√(x³) + 1
                IntegrationUtils.funcE()  // (x³ - 2x² - 12x - 7) / (x³ - x²)
        );

        double[][] intervals = {
                {5.0, 10.0},            // для funcA: x > 4 (область определения sqrt(x² - 16))
                {1.0, 3.0},             // для funcB: ln(x) требует x > 0, избегаем x=0
                {0.0, Math.PI / 2},     // для funcC: функция периодическая, выбран интервал без особенностей
                {0.0, 2.0},             // для funcD: x >= 0, подкоренное выражение > 0
                {2.0, 5.0}              // для funcE: исключаем x=0 и x=1 (особенности в знаменателе)
        };

        String[] functionNames = {"A", "B", "C", "D", "E"};
        String[] methodNames = {"Прямоугольники", "Трапеции", "Симпсон"};

        // Для каждого метода интегрирования
        for (int methodIdx = 0; methodIdx < 3; methodIdx++) {
            System.out.printf("\n=== МЕТОД: %s ===\n", methodNames[methodIdx]);

            // Для каждой функции
            for (int funcIdx = 0; funcIdx < functions.size(); funcIdx++) {
                Function<Double, Double> f = functions.get(funcIdx);
                double a = intervals[funcIdx][0];
                double b = intervals[funcIdx][1];

                System.out.printf("\n--- Функция %s на [%.1f, %.1f] ---\n",
                        functionNames[funcIdx], a, b);

                // Для каждой точности
                for (double eps : ACCURACIES) {
                    System.out.printf("\nТочность: %.0e\n", eps);

                    int optimalN = findOptimalN(f, a, b, eps, methodIdx);
                    System.out.printf("Оптимальное N: %,d\n", optimalN);

                    // Замеряем время для последовательной версии
                    int finalMethodIdx1 = methodIdx;
                    long seqTime = measureTime(() -> {
                        for (int i = 0; i < REPEATS; i++) {
                            calculateSequential(f, a, b, optimalN, finalMethodIdx1);
                        }
                    }, REPEATS);

                    System.out.printf("Последовательный: %.9f мс\n", seqTime / 1_000_000.0);

                    // Замеряем время для параллельных версий
                    for (int threads : THREAD_COUNTS) {
                        int finalMethodIdx = methodIdx;
                        long parTime = measureTime(() -> {
                            for (int i = 0; i < REPEATS; i++) {
                                calculateParallel(f, a, b, optimalN, finalMethodIdx, threads);
                            }
                        }, REPEATS);

                        System.out.printf("Параллельный (%d потоков): %.9f мс\n",
                                threads, parTime / 1_000_000.0);
                    }
                }
            }
        }
    }

    private static int findOptimalN(Function<Double, Double> f, double a, double b,
                                    double eps, int methodIdx) {
        int N = 100;
        for (int i = 0; i < 20; i++) {
            double result1 = calculateSequential(f, a, b, N, methodIdx);
            double result2 = calculateSequential(f, a, b, N * 2, methodIdx);
            if (Math.abs(result1 - result2) < eps) return N * 2;
            N *= 2;
            // if (N > 10_000_000) break;
        }
        return N;
    }

    private static double calculateSequential(Function<Double, Double> f, double a,
                                              double b, int n, int methodIdx) {
        switch (methodIdx) {
            case 0: return Integrator.rectangle(f, a, b, n);
            case 1: return Integrator.trapezoid(f, a, b, n);
            case 2: return Integrator.simpson(f, a, b, n);
            default: throw new IllegalArgumentException();
        }
    }

    private static double calculateParallel(Function<Double, Double> f, double a,
                                            double b, int n, int methodIdx, int threads) {
        switch (methodIdx) {
            case 0: return ParallelIntegrator.rectangle(f, a, b, n, threads);
            case 1: return ParallelIntegrator.trapezoid(f, a, b, n, threads);
            case 2: return ParallelIntegrator.simpson(f, a, b, n, threads);
            default: throw new IllegalArgumentException();
        }
    }

    private static long measureTime(Runnable task, int repeats) {
        for (int i = 0; i < 100; i++) {
            task.run();
        }

        long start = System.nanoTime();
        task.run();
        long end = System.nanoTime();

        return (end - start) / repeats;
    }


}
