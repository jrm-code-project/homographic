# homographic
Homographic functions, Möbius transforms, and exact real arithmetic.

## Overview
This project provides tools for working with homographic functions, Möbius transformations, and exact real arithmetic. It includes implementations of linear fractional transformations, binary expressions, and streams for computation.

## Files and Functionality

### `lft.lisp`
This file implements Linear Fractional Transformations (LFTs), also known as homographic functions.  LFTs are central to the computation of exact real arithmetic.

### `binary-expression.lisp`
This file provides utilities for working with binary expressions. It includes functions for parsing, evaluating, and manipulating binary representations of numbers.

### `cf-stream.lisp`
This file implements continued fraction streams, which are used to represent real numbers as infinite sequences. These streams are essential for exact real arithmetic and provide a compact representation of irrational numbers.

### `lft-stream.lisp`
This file extends the functionality of LFTs by integrating them with streams. It enables the computation of real numbers using lazy evaluation and stream-based transformations.

### `lft-arith.lisp`
This file provides arithmetic operations for LFTs, including addition, subtraction, multiplication, and division. These operations are implemented to work seamlessly with the lazy evaluation model.

## Usage
To use the provided files, load the `linear-fractional-transformation.asd` system definition file into your Lisp environment. This will load all the necessary modules. You can then call the relevant functions. Each file is modular and can be used independently or in combination with others for advanced computations.a

## Examples
Examples of usage can be found in the respective files. Refer to the comments and function definitions for detailed explanations.
