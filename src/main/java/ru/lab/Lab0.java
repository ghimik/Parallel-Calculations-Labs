package ru.lab;

public class Lab0 {

    public Runnable getExampleRunnable() {
        return () -> System.out.println("Hello from thread " + Thread.currentThread().getName());
    }

    public void run() throws InterruptedException {
        Thread[] threads = new Thread[10];

        for (int i = 0; i < 10; i++) {
            Thread thread = new Thread(getExampleRunnable());
            threads[i] = thread;
        }

        for (Thread thread : threads) {
            thread.start();
        }

        for (Thread thread : threads) {
            thread.join();
        }
    }

    public void runSorts() {

    }
}
