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

