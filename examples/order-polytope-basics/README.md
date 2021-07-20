## Compilation
In folder examples, first run cmake, to create the makefile:

```bash
cmake .
```

Then, in folder examples/order-polytope-basics compile and build using the makefile:

```bash
make
```

## Usage:
```bash
 ./order_polytope poset_data.txt
```  
where `poset_data.txt` is the file containing the poset data as follows:
- first line of file tells the number of elements in the poset data
- next `m` lines represent the order relations of the poset as a pair of indices, for example:  
a line containing `i j` means A<sub>i</sub> <= A<sub>j</sub>. ( 0 <= `i`, `j` < n)
    