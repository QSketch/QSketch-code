# Qsketch - Quantile Sketch

​	The key techniques of QSketch are Value Focus, Distribution Calibration and Double Filtration. The key idea of Value Focus is to maintain the Candidate and the Representative to record values close to the target quantile. The key idea of Distribution Calibration is to insert positive (negative) infinities into the data structure by probability method. In this way, the 𝑤-quantile of the original distribution is just the 0.5-quantile of the calibrated distribution, and we convert the quantile estimation problem into the median estimation problem. The key idea of Double Filtration is to filter infrequent items in both stages of QSketch to make room for frequent items.

## Getting Started

To run the experiment of an algorithm, first go to the corresponding folder, for example, if you want to test per-key Qsketch, 

```
cd Qsketch-main
```

Then, to compile the code:

```
make
```

To run Qsketch on the benchmarks:

```
./main
```

To run test algorithm for choosing parameters:

```
./test
```