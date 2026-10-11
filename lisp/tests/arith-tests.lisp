;;; -*- Lisp -*-

;;; PACKAGE-BOUNDARY POLICY: see tests/package.lisp and tests/tests.lisp
;;; for why this suite lives in LFT-TESTS and reaches into
;;; LINEAR-FRACTIONAL-TRANSFORMATION's internal (non-exported) symbols
;;; via the `lft::foo` prefix (e.g. `lft::lft-stream->double`,
;;; `lft::->lft-stream`, `lft::log2`, `lft::pi/3`).
;;;
;;; This file exercises lft-arith.lisp, which implements transcendental
;;; functions (exp, log, expt, trigonometric, inverse trigonometric, and
;;; hyperbolic functions) over both exact rationals and infinite
;;; continued-fraction streams (LFT-STREAM objects).
;;;
;;; NOTE ON RESULT SHAPES: many of these functions return a plain
;;; rational/integer for exact inputs (e.g. (exp 0) => 1, (cos 0) => 1,
;;; (asin 0) => 0) but return an LFT-STREAM for inexact/irrational
;;; results. The APPROX= helper below handles both uniformly by only
;;; calling LFT::LFT-STREAM->DOUBLE on actual LFT-STREAM objects.
;;;
;;; NOTE ON A KNOWN, DELIBERATELY-UNFIXED LIMITATION: converting an
;;; LFT-STREAM that was built from an *exact* rational (e.g.
;;; (lft::->lft-stream 1), or the stream that (sqrt 1) would produce)
;;; back to a double via LFT::LFT-STREAM->DOUBLE hangs indefinitely; this
;;; is a pre-existing defect in stream refinement/termination that is
;;; broader than lft-arith.lisp and is out of scope here. Tests
;;; deliberately avoid constructing exact-valued streams directly (e.g.
;;; they call (acos 0) rather than (acos (lft::->lft-stream 0))).

(in-package "LFT-TESTS")

(in-suite lft-suite)

(defconstant +epsilon+ 1d-9
  "Absolute tolerance used when comparing computed double-float
approximations of lft-arith.lisp results against reference values.")

(defun approx= (actual expected &optional (epsilon +epsilon+))
  "Return true if ACTUAL (a rational, integer, or LFT-STREAM) converts to
a double-float within EPSILON of the double-float EXPECTED."
  (let ((actual-double (if (numberp actual)
                            (coerce actual 'double-float)
                            (lft::lft-stream->double actual))))
    (< (abs (- actual-double expected)) epsilon)))

;;; -------------------------------------------------------------------
;;; exp
;;; -------------------------------------------------------------------

(test exp-exact
  (is (= 1 (exp 0))))

(test exp-general
  (is (approx= (exp 1) 2.718281828459045d0))
  (is (approx= (exp -1) 0.36787944117144233d0))
  (is (approx= (exp 2) 7.38905609893065d0)))

;;; -------------------------------------------------------------------
;;; log (including the divide-down-to-exactly-1 fix for powers of two)
;;; -------------------------------------------------------------------

(test log-exact
  (is (= 0 (log 1))))

(test log-general
  (is (approx= (log 2) 0.6931471805599453d0))
  (is (approx= (log 3) 1.0986122886681098d0))
  (is (approx= (log 4) 1.3862943611198906d0))
  (is (approx= (log 1/2) -0.6931471805599453d0)))

(test log-powers-of-two-regression
  ;; LOG-RAT's DIVIDE-DOWN helper used to call %LOG-RAT on exactly 1 when
  ;; dividing an exact power of two all the way down, which signalled a
  ;; SIMPLE-TYPE-ERROR because %LOG-RAT's domain excludes 1. Any power of
  ;; two should log cleanly now.
  (is (approx= (log 8) 2.0794415416798357d0))
  (is (approx= (log 1/4) -1.3862943611198906d0)))

;;; -------------------------------------------------------------------
;;; expt / rat-pow
;;; -------------------------------------------------------------------

(test expt-integer-exponents
  (is (= 8 (expt 2 3)))
  (is (= 1 (expt 0 0)))
  (is (= 0 (expt 0 5)))
  (is (= 1/2 (expt 2 -1))))

(test expt-rational-exponent
  (is (approx= (expt 2 1/2) 1.4142135623730951d0)))

;;; -------------------------------------------------------------------
;;; sqrt / named constants
;;; -------------------------------------------------------------------

(test sqrt-general
  (is (approx= (sqrt 2) 1.4142135623730951d0)))

(test named-constants
  (is (approx= pythagoras 1.4142135623730951d0))
  (is (approx= theodorus 1.7320508075688772d0))
  (is (approx= sqrt-five 2.23606797749979d0))
  (is (approx= pi 3.141592653589793d0)))

;;; -------------------------------------------------------------------
;;; tan
;;; -------------------------------------------------------------------

(test tan-exact
  (is (= 0 (tan 0))))

(test tan-general
  (is (approx= (tan 1) 1.5574077246549023d0)))

;;; -------------------------------------------------------------------
;;; cos / sin
;;; -------------------------------------------------------------------

(test cos-sin-exact
  (is (= 1 (cos 0)))
  (is (= 0 (sin 0))))

(test cos-sin-general
  (is (approx= (cos 1) 0.5403023058681398d0))
  (is (approx= (sin 1) 0.8414709848078965d0)))

;;; -------------------------------------------------------------------
;;; atan (one- and two-argument forms, including the zero-numerator fix)
;;; -------------------------------------------------------------------

(test atan-general
  (is (approx= (atan 1) 0.7853981633974483d0))      ; pi/4
  (is (approx= (atan 1 1) 0.7853981633974483d0))
  (is (approx= (atan -1 1) -0.7853981633974483d0)))

(test atan-zero-numerator-regression
  ;; ATAN-RAT used to fall through to (%RAT-ATAN (/ Y X)) = (%RAT-ATAN 0)
  ;; when Y was exactly 0 and X was positive, but %RAT-ATAN's domain
  ;; excludes 0, signalling a type error. atan(0, positive) should be 0;
  ;; atan(0, negative) should still correctly yield pi via the existing
  ;; MINUSP X branch.
  (is (= 0 (atan 0 1)))
  (is (approx= (atan 0 -1) 3.141592653589793d0)))

;;; -------------------------------------------------------------------
;;; acos / asin (sign-correctness + division-by-zero fixes)
;;; -------------------------------------------------------------------

(test acos-exact-boundaries
  (is (= 0 (acos 1)))
  (is (approx= (acos -1) 3.141592653589793d0))
  (is (approx= (acos 0) 1.5707963267948966d0)))        ; pi/2

(test acos-general
  (is (approx= (acos 1/2) 1.0471975511965979d0))
  (is (approx= (acos -1/2) 2.0943951023931957d0)))

(test asin-exact-boundaries
  (is (= 0 (asin 0))))

(test asin-general
  (is (approx= (asin 1/2) 0.5235987755982989d0)))      ; pi/6

;;; -------------------------------------------------------------------
;;; hyperbolic functions
;;; -------------------------------------------------------------------

(test hyperbolic-exact
  (is (= 0 (sinh 0)))
  (is (= 1 (cosh 0)))
  (is (= 0 (tanh 0)))
  (is (= 0 (atanh 0)))
  (is (= 1 (lft::sech 0))))

(test hyperbolic-general
  (is (approx= (sinh 1) 1.175201177597046d0))
  (is (approx= (cosh 1) 1.5430806875228882d0))
  (is (approx= (tanh 1) 0.7615941762924194d0))
  (is (approx= (atanh 1/2) 0.5493061443340548d0))
  (is (approx= (lft::coth 1) 1.313035249710083d0))
  (is (approx= (lft::csch 1) 0.8509181141853333d0)))

;;; -------------------------------------------------------------------
;;; log-base
;;; -------------------------------------------------------------------

(test log-base-general
  (is (approx= (lft::log-base 8 2) 3.0d0)))

;;; -------------------------------------------------------------------
;;; square / cbrt (via lft-stream, since there is no rational CBRT method)
;;; -------------------------------------------------------------------

(test square-general
  (is (= 9 (square 3)))
  (is (approx= (square pythagoras) 2.0d0)))

(test cbrt-via-stream
  (is (approx= (lft::cbrt (sqrt 2)) 1.122462048309373d0)))
