# Copilot instructions for `homographic`

## What this is
Common Lisp library implementing exact real arithmetic via Linear Fractional
Transformations (LFTs / Möbius transforms) and continued-fraction streams,
following Peter Potts' thesis and Gosper's algorithms for arithmetic on
continued fractions. All source lives in `lisp/`.

## Build / load
This is an ASDF system (`lisp/linear-fractional-transformation.asd`), not a
standalone project with a test runner or lint config — there is none in this
repo. It depends on several sibling systems that are **not** part of this
repo: `fold`, `generic-arithmetic`, `named-let`, `promise`, `series`,
`stream`, `utilities`. These must be reachable through ASDF's source registry
(e.g. as sibling checkouts registered with `asdf:*central-registry*` or
`~/common-lisp/`) before this system will load. To load it from a REPL:

```lisp
(asdf:load-system "linear-fractional-transformation")
```

There are no automated tests in this repo; validate changes by loading the
system and exercising functions interactively (e.g. via SBCL/`sly`/`slime`).

## Architecture (read these in dependency order)
The `.asd` file's `:components` list is also the dependency/read order:
1. **`package.lisp`** — defines the single package `LINEAR-FRACTIONAL-TRANSFORMATION`
   (nickname `LFT`). It shadows `COMMON-LISP` arithmetic operators with
   `GENERIC-ARITHMETIC` versions, shadows `LET`/`LET*`/`DEFUN`/`FUNCALL` with
   `NAMED-LET`/`SERIES` versions, and shadows `STREAM`. When adding exported
   symbols, add them to the `:export` list here.
2. **`lft.lisp`** — core `LFT` (2x2 matrix, `(ax+b)/(cx+d)`) and `BILFT`
   (bilinear, two-variable) classes, their constructors (`make-lft`,
   `make-bilft`), composition, evaluation (`evaluate-lft`/`evaluate-bilft`),
   and matrix-comparison helpers used to decide when a digit can be emitted.
3. **`lft-stream.lisp`** — `LFT-STREAM`: a lazy stream of LFTs (built on the
   `STREAM` package's generic stream), representing a real number as an
   infinite composition of transformations. Converts rationals/floats to/from
   streams and implements the "refine" (consume more input) / "emit digit"
   algorithm at the heart of exact real arithmetic.
4. **`cf-stream.lisp`** — `CF-STREAM`: continued-fraction streams, convertible
   to/from `LFT-STREAM`.
5. **`binary-expression.lisp`** — `BINARY-EXPRESSION`: lazy binary operations
   (two `LFT-STREAM` operands plus a `BILFT`) that refine whichever operand
   side yields more precision, per Gosper's algorithm.
6. **`lft-arith.lisp`** — generic arithmetic methods (`+`, `*`, `exp`, `log`,
   `sqrt`, trig/hyperbolic functions, etc.) dispatching on `LFT-STREAM`,
   rationals, and floats, implemented in terms of the lower layers above.

Dispatch is almost entirely through CLOS `defmethod`/`defgeneric` on numeric
and stream types (`cl:rational`, `single-float`, `double-float`,
`lft-stream`, `(eql 'infinity)`, etc.) rather than `etypecase`/`ecase`.

## Conventions used throughout this codebase
- **Predicates end in `?`**, not `-p` (e.g. `unimodular?`, `disjoint?`,
  `2x2-matrix-range?`, `empty-stream?`).
- **`%`-prefixed functions are internal/low-level**, puncturing an
  abstraction or carrying non-obvious preconditions (e.g. `%make-lft`,
  `%rat->lft-stream`, `%sqrt-rat`). Don't export these or add `%` to
  ordinary helpers.
- **CLOS methods for a generic are grouped as `:method` clauses inside the
  `defgeneric` form** when a generic is first introduced locally (see
  `->cf-stream`, `->lft-stream`, `funcall-lft`, `funcall-bilft` in
  `lft.lisp`); later specializations elsewhere are separate `defmethod`
  forms.
- **`check-type` at the top of functions** validates integer/rational
  argument preconditions (see `evaluate-lft`, `evaluate-bilft`,
  `format-term`).
- Symbol/package designators in `defpackage`/`in-package` use literal
  uppercase strings (e.g. `(in-package "LINEAR-FRACTIONAL-TRANSFORMATION")`),
  not unquoted symbols.
- Classes use `:car`/`:delayed-cdr` initargs and `delay`/`force` (from the
  `PROMISE` system) for laziness, matching the `STREAM` package's stream
  protocol — follow this pattern when adding new lazy stream types rather
  than inventing ad hoc laziness.
- No `loop` macro and no imperative iteration; this repo follows the
  project-wide Lisp style of tail-recursive named `let`, `fold-left`, and
  CLOS dispatch over branching (see the `lisp-functional` skill for the full
  rule set if available in your environment).
