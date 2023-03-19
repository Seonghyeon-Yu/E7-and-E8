# The-Betti-numbers-of-type-E7-and-E8

Python and SageMath codes for our paper :

_"The Betti numbers of real toric varieties associated to Weyl chambers of types_ $E_{7}$ _and_ $E_{8}$"

# Environment

We have worked on _"SageMath 9.3 Notebook"_ ; 

if you want to run these codes, need to be installed _"SageMath 9.3"_.

We used the Python packages : __numpy, math, pickle5, itertools__

and used the SageMath package : __SimplicialComplex__

# Codes

## Basic Tools

### Type(simple_root)

```python
return [["representatives of orbit types of reduced row space of characteristic matrix"], ["number of orbit types"]]
```


### reflection(u,v)

```python
return ["the integral primitive vector directed to the vector where v is reflected by the hyperplane perpendicular to u"]
```



### coset_representation(simple_root, fundamental_coweight)

``` python
return ["coset representations of Weyl group"]
```



### vtx(simple_root, fundamental_coweight)

```python
return ["coweight lattice coordinates of vertices of Coxeter complex"]
```



### vtxpre(simple_root, fundamental_coweight)

```python
return [["vertices classified by fundamental coweight"]]
```



### char_coordinate(vector, simple_root)

````python
return ["coordinate of vector in characteristic matrix"]
````



### find_reduced_rowspace(simple_root)

```python
return [["elements of reduced rowspace classified by orbit types"]]
```



### find_simplex(simple_root, fundamental_coweight)

```python
return ["simplices of Coexter complex"]
```



## Basic version

This is the basic version model to basic unrefined set for each orbit types .

### maximal_simplex_basic(simple_root, fundamental_coweight)

```python
pickle.dump(["maximal simplices of a type"], 'simplices of {name} of type {type number} basic ver.pkl', 'wb')
```



## Coset version

This is the 'coset method' version model.

The results are same as the results of basic version, but we reduce the time complexity using by mathematical theory.

### maximal_simplex_coset(simple_root, fundamental_coweight)

```python
for num in range(1, coset_representation):
	pickle.dump(["maximal simplices of a type"], 'simplices of {name} of type {type number} No.{num} coset ver.pkl', 'wb')


pickle.dump(["maximal simplices of a type"], 'simplices of {name} of type {type number} coset ver.pkl', 'wb')
```







#Data

