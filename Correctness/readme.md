# To check the correctness
These are to check the correctness of our codes.

The known results for types $G_2$, $F_4$ and $E_6$ can be found in $Table 1$ our paper.



The results from our codes are as follows :

## $G_2$

```{python}
not connected : number of components are 4
Number of orbit : 3
Number of vtx of coset ver. : 8
Type 1 is pure : True
f-vector : [1, 2, 1]
Betti number : {0: 1, 1: 0}
--------------------------
```

## $F_4$

```{python}
not connected : number of components are 2
Number of orbit : 12
Number of vtx of coset ver. : 141
Number of vtx of vtx elimination : 28
Type 1 vtx elimination is pure : True
f-vector of vtx elimination : [1, 14, 24]
Betti number : {0: 1, 1: 11}
--------------------------
not connected : number of components are 16
Number of orbit : 3
Number of vtx of coset ver. : 80
Number of vtx of vtx elimination : 16
Type 2 vtx elimination is pure : True
f-vector of vtx elimination : [1, 1]
Betti number : {0: 1}
--------------------------
```

## $E_6$

```{python}
not connected : number of components are 2
Number of orbit : 36
Number of vtx of coset ver. : 704
Number of vtx of vtx elimination : 64
Type 1 vtx elimination is pure : True
f-vector of vtx elimination : [1, 32, 150, 180]
Betti number : {0: 1, 1: 0, 2: 61}
--------------------------
connected
Number of orbit : 27
Number of vtx of coset ver. : 576
Number of vtx of vtx elimination : 32
Type 2 vtx elimination is pure : True
f-vector of vtx elimination : [1, 32, 80]
Betti number : {0: 1, 1: 49}
--------------------------
```
