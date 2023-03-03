# DNA-Sequencing

DNA sequence is treated as a "language", otherwise known as k-mer counting.
The long biological sequence briken down into length 6
overlapping words (hexamers).
[data.csv](data.csv) is obtained by applying BAG
of WORDS on the hexamers using count vectorizer. Our aim is to perform Bayes’ classification on the
feature vectors thus obtained and parallelize the classification task.

## HPC APIs used for parallelization
<ul>
  <li> OpenMP </li>
  <li> MPI </li>
  <li> CUDA </li>
</ul>

## Profiling of serial code
### Code
[se.cpp](se.cpp)

### Functional profiling- gprof:
From the flat profile we can observe that the highest percentage of the total
execution time the program spent is in matrix and vector operations. These functions that takes significant time are covar, shift, dot, dot1, mean, transpose and
transpose1.  
The call graph shows how much time was spent in each function and its children.
From this information, we can find functions that, while they themselves may not
have used much time, called other functions that did use unusual amounts of time.
Both fit and transform functions are called one time only. Referring covar function by
index in the call graph, it is clear that it takes the highest percentage of the total
execution time and called from main function.

### Line profiling- gcov:
From line profiling we can observe different hotspots in the program. Line 67 with
loop at line 65 in shift function is a hotspot but low at hotness scale. Line 51 with
loop at line 49 in mean function has medium hotness. Line 113 with loop at line 111
in covar function has high hotness. Also, line 157 with loop at line 155 in dot function
has high hotness.  
The main functions call the expensive functions a number of times using for loops.
High percentage of branches are taken as the condition of loops only fails once.

## Hardware resource profiling- likwid
The machine on which we are running our program has a single socket with six cores
and there are two threads per core. Each core has private L1 cache of size 32 kB and
private L2 cache of size 256 kB. L3 cache is shared with size of 12 MB.

## OpenMP
### Code
[new_omp.cpp](new_omp.cpp)
### Thread vs Time
Thread | Time
--- | ---
1 | 23.7732
2 | 11.9414
4 | 6.36412
6 | 6.12247
8 | 6.04035
10 | 5.52894
12 | 5.14293
16 | 7.06707
24 | 7.66674
32 | 8.63845
64 | 11.9971
128 | 18.8138

### Speedup vs Processors
Processors | Speedup
--- | ---
1 | 1
2 | 1.990821847
4 | 3.735504673
6 | 3.882942669
8 | 3.935732201
10 | 4.299775364
12 | 4.622501181
16 | 3.363940077
24 | 3.100822514
32 | 2.752021485
64 | 1.981578882
128 | 1.263604376

### Parallelization fraction
Thread | Time | f(Parallelization fraction)
--- | --- | ---
2 | 11.9414 | 0.995389767
4 | 6.36412 | 0.976398073
6 | 6.12247 | 0.890956035
8 | 6.04035 | 0.852477339
10 | 5.52894 | 0.852699679
12 | 5.14293 | 0.854909348
16 | 7.06707 | 0.74957818
24 | 7.66674 | 0.706961657
32 | 8.63845 | 0.657167219
64 | 11.9971 | 0.503214638
128 | 18.8138 | 0.210255683

## MPI
### Code
[new_mpi.cpp](new_mpi.cpp)

## CUDA
### Code
[new_cuda.ipynb](new_omp.ipynb)
### Thread vs Time
Thread | Time
--- | ---
1 | 22.0034
2 | 6.86
4 | 2.45065
6 | 2.4258
8 | 1.74086
10 | 1.89988
12 | 1.75416
16 | 1.5273
24 | 1.55547
32 | 1.52225
64 | 0.798186
128 | 0.804836

### Speedup vs Processors
Processors | Speedup
--- | ---
1 | 1
2 | 3.207492711
4 | 8.978597515
6 | 9.070574656
8 | 12.63938513
10 | 11.5814683
12 | 12.54355361
16 | 14.40673083
24 | 14.14582088
32 | 14.45452455
64 | 27.56675762
128 | 27.33898583

### Parallelization fraction
Thread | Time | f(Parallelization fraction)
--- | --- | ---
2 | 6.86 | 1.376460002
4 | 2.45065 | 1.184832041
6 | 2.4258 | 1.067704082
8 | 1.74086 | 1.052436831
10 | 1.89988 | 1.015172403
12 | 1.75416 | 1.003939391
16 | 1.5273 | 0.9926272
24 | 1.55547 | 0.969712428
32 | 1.52225 | 0.960843881
64 | 0.798186 | 0.979021637
128 | 0.804836 | 0.9710082
