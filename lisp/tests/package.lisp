;;; -*- Lisp -*-

;;; LINEAR-FRACTIONAL-TRANSFORMATION-TESTS is a separate ASDF
;;; system/package from LINEAR-FRACTIONAL-TRANSFORMATION (rather than
;;; tests living inside the main package) so that the test suite -- and
;;; its dependency on FIVEAM -- is never loaded into a production image
;;; via a plain `(asdf:load-system "linear-fractional-transformation")`;
;;; it is only pulled in by the
;;; `linear-fractional-transformation/tests` system definition, which
;;; `(asdf:test-system "linear-fractional-transformation")` loads on
;;; demand.
;;;
;;; Tests here deliberately reach into LINEAR-FRACTIONAL-TRANSFORMATION's
;;; internal (non-exported) symbols via the `lft::foo` package-prefix
;;; syntax throughout tests.lisp (e.g. `lft::unimodular?`,
;;; `lft::evaluate-lft`, `lft::canonicalize-lft-coefficients`) rather than
;;; only exercising the package's public API (MAKE-LFT, MAKE-BILFT, and
;;; friends). This is an intentional white-box-testing convention, not an
;;; oversight: most of the interesting logic in lft.lisp -- coefficient
;;; canonicalization, matrix/tensor comparison predicates, digit
;;; truncation, sign/range predicates -- is deliberately unexported
;;; plumbing, and testing it directly lets the suite pin down its exact
;;; behavior without first building a full continued-fraction stream for
;;; every case. The tradeoff is that the test suite is tightly coupled to
;;; internal names/shapes and will need updating on any internal rename,
;;; even when LFT's public contract hasn't changed.

(defpackage "LINEAR-FRACTIONAL-TRANSFORMATION-TESTS"
  (:nicknames "LFT-TESTS")
  (:shadowing-import-from "FUNCTION" "COMPOSE" "INVERSE")
  (:shadowing-import-from "GENERIC-ARITHMETIC"
                          "+" "-" "*" "/" "<" "<=" "=" ">" ">=" "/="
                          "1+" "1-"
                          "ABS"
                          "ACOS"
                          "ACOSH"
                          "ASIN"
                          "ASINH"
                          "ATAN"
                          "ATANH"
                          "CEILING"
                          "CIS"
                          "COMPLEX"
                          "CONJUGATE"
                          "COS"
                          "COSH"
                          "DECODE-FLOAT"
                          "DENOMINATOR"
                          "EXP"
                          "EXPT"
                          "FCEILING"
                          "FFLOOR"
                          "FLOAT"
                          "FLOAT-DIGITS"
                          "FLOAT-PRECISION"
                          "FLOAT-RADIX"
                          "FLOAT-SIGN"
                          "FLOOR"
                          "FROUND"
                          "FTRUNCATE"
                          "IMAGPART"
                          "INTEGER-DECODE-FLOAT"
                          "LOG"
                          "MAX"
                          "MIN"
                          "MOD"
                          "NUMERATOR"
                          "PHASE"
                          "RANDOM"
                          "RATIONAL"
                          "RATIONALIZE"
                          "REALPART"
                          "RECIPROCAL"
                          "REM"
                          "ROUND"
                          "SCALE-FLOAT"
                          "SIGNUM"
                          "SIN"
                          "SINH"
                          "SQRT"
                          "TAN"
                          "TANH"
                          "TRUNCATE"
                          "ZEROP")
  (:shadowing-import-from "LINEAR-FRACTIONAL-TRANSFORMATION" "PI")
  (:use "COMMON-LISP" "FIVEAM" "FUNCTION" "GENERIC-ARITHMETIC"
        "LINEAR-FRACTIONAL-TRANSFORMATION")
  (:export "LFT-SUITE"))
