# Homepage

PetIGA is a software framework that approximates the solution of partial differential equations using isogeometric analysis. It is an extension of PETSc, adding the NURBS discretization capability and the integration of forms. This framework can be used to solve linear, nonlinear, time-dependent, or time-dependent nonlinear problems. In this framework, the user has to provide the evaluation of the linear form (right-hand side, or residual of a nonlinear problem) at a Gauss point, as well as the bilinear form (left-hand side or Jacobian of the nonlinear residual). The framework is designed so that researchers can focus on the physics of the problem and ignore issues of parallelism and performance.


## Features of PetIGA


* Flexibility and Robustness: PetIGA supports various nonlinear and linear problems in solid and fluid mechanics, demonstrating strong scaling on up to 4096 cores.

* Parallel Scalability: The framework is well-suited for large-scale applications, ensuring efficient parallel scalability.

* Support for Numerical Differentiation: PetIGA includes support for numerical differentiation, which is essential for evaluating the Jacobians.

* Open-Source and Extensible: The framework is open-source and designed to be reusable, promoting code reuse and flexibility in scientific research.

## Objectives

This documentation provides an overview of the PetIGA framework, including its features, capabilities, and applications. It is intended for beginners who want to understand the usage of PetIGA. The documentation covers some simple tutorials and examples to help users get started with the framework and apply it to solve real-world problems in solid and fluid mechanics.

### Codeblocks

Some `code` goes here.

### Plain codeblock

A plain codeblock:

```
Some code here
def myfunction()
// some comment
```

#### Code for a specific language

Some more code with the `py` at the start:

``` py
import tensorflow as tf
def whatever()
```

#### With a title

``` py title="bubble_sort.py"
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```

#### With line numbers

``` py linenums="1"
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```

#### Highlighting lines

``` py hl_lines="2 3"
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```

## Icons and Emojs

:smile: 

:fontawesome-regular-face-laugh-wink:

:fontawesome-brands-twitter:{ .twitter }

:octicons-heart-fill-24:{ .heart }