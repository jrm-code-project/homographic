;;; -*- Lisp -*-

(defpackage "LINEAR-FRACTIONAL-TRANSFORM"
  (:nicknames "LFT")
  (:shadowing-import-from "SERIES"
                          "DEFUN"
                          "FUNCALL"
                          ;; "LET"
                          "LET*"
                          "MULTIPLE-VALUE-BIND")
  (:shadowing-import-from "STREAM" "SCAN-STREAM" "STREAM")
  (:shadowing-import-from "UTILITIES" "LET")
  (:shadow "PI")
  (:use "COMMON-LISP" "PROMISE" "SERIES" "STREAM" "UTILITIES")
  (:export
   "->CF-STREAM"
   "2PI"
   "BIG-K-STREAM"
   "BILFT"
   "CF-STREAM"
   "COMPOSE"
   "CONS-CF-STREAM"
   "CONS-LFT-STREAM"
   "EULER-GOMPERTZ"
   "INFINITY"
   "INVERSE"
   "LFT"
   "LFT-STREAM"
   "LIMIT-STREAM->CF-STREAM"
   "MAKE-BILFT"
   "MAKE-LFT"
   "NEGATE"
   "PHI"
   "PI"
   "PI/2"
   "RECIPROCAL"
   "SQRT-TWO"
   "X*"
   "X+"
   "X-"
   "X/"
   "X-CBRT"
   "X-COSH"
   "X-EXP"
   "X-EXPT"   
   "X-LOG"
   "X-SINH"
   "X-SQRT"
   "X-TANH"
   ))

