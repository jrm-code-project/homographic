;;; -*- Lisp -*-

;;; PACKAGE-BOUNDARY POLICY: see tests/package.lisp and tests/tests.lisp
;;; for why this suite lives in LFT-TESTS and reaches into
;;; LINEAR-FRACTIONAL-TRANSFORMATION's internal (non-exported) symbols
;;; via the `lft::foo` prefix.
;;;
;;; Unlike the other test files (tests.lisp, arith-tests.lisp,
;;; lft-stream-tests.lisp, binary-expression-tests.lisp,
;;; cf-stream-tests.lisp), each of which exercises a single source file
;;; in relative isolation, this file is deliberately END-TO-END and
;;; CROSS-FUNCTIONAL: every test here chains together machinery from at
;;; least three of {lft.lisp, lft-stream.lisp, binary-expression.lisp,
;;; cf-stream.lisp, lft-arith.lisp} in a single test, the way a real
;;; caller would -- round-tripping values between the LFT-STREAM and
;;; CF-STREAM representations, composing BINARY-EXPRESSION-based
;;; arithmetic (+, -, *) with LFT-ARITH's transcendental functions, and
;;; checking mathematical identities (Pythagorean, inverse-function,
;;; hyperbolic) survive the full pipeline from an irrational LFT-STREAM
;;; input down to a final double-float or exact rational answer.
;;;
;;; This file depends on (and is loaded after) arith-tests.lisp,
;;; lft-stream-tests.lisp, binary-expression-tests.lisp, and
;;; cf-stream-tests.lisp so it can reuse their helpers (APPROX=,
;;; STREAM-APPROX=, CF-STREAM->LIST) and named LFT-STREAM fixtures
;;; (SQRT-TWO, PHI) rather than redefining them.

(in-package "LFT-TESTS")

(in-suite lft-suite)

;;; -------------------------------------------------------------------
;;; Exact rational round trips: ->LFT-STREAM -> ->CF-STREAM ->
;;; CF-STREAM->RATIONAL must reproduce the original value exactly
;;; (lft.lisp + lft-stream.lisp + cf-stream.lisp).
;;; -------------------------------------------------------------------

(test round-trip-rational-through-lft-stream-and-cf-stream
  (dolist (rat '(22/7 -22/7 5 -5 0 1/3))
    (is (= rat (lft::cf-stream->rational
                (lft::->cf-stream (lft::->lft-stream rat)))))))

;;; -------------------------------------------------------------------
;;; Irrational round trip: an LFT-STREAM converted to CF-STREAM and back
;;; to LFT-STREAM (via a different code path than it started on) still
;;; converges to the same double (lft-stream.lisp + cf-stream.lisp).
;;; -------------------------------------------------------------------

(test round-trip-irrational-through-cf-stream
  (is (approx= (lft::lft-stream->double
                (lft::cf-stream->lft-stream (lft::->cf-stream sqrt-two)))
               1.4142135623730951d0))
  (is (approx= (lft::lft-stream->double
                (lft::cf-stream->lft-stream (lft::->cf-stream lft::pythagoras)))
               1.4142135623730951d0))
  (is (approx= (lft::lft-stream->double
                (lft::cf-stream->lft-stream (lft::->cf-stream lft::euler-gompertz)))
               0.5963473623231941d0)))

;;; -------------------------------------------------------------------
;;; BINARY-EXPRESSION's + dispatch result round-tripped through
;;; CF-STREAM (binary-expression.lisp + cf-stream.lisp).
;;; -------------------------------------------------------------------

(test binary-expression-sum-round-trips-through-cf-stream
  (is (approx= (lft::lft-stream->double
                (lft::cf-stream->lft-stream (lft::->cf-stream (+ sqrt-two phi))))
               3.03224755112299d0)))

;;; -------------------------------------------------------------------
;;; An exact-rational BINARY-EXPRESSION sum both prints correctly (the
;;; print-hang fix) and round-trips to the exact same rational via
;;; CF-STREAM (binary-expression.lisp + lft-stream.lisp + cf-stream.lisp).
;;; -------------------------------------------------------------------

(test exact-binary-expression-sum-prints-and-round-trips-exactly
  (let ((sum (+ (lft::->lft-stream 1/2) (lft::->lft-stream 1/3))))
    (is (search "0.8333333333333" (format nil "~a" sum)))
    (is (= 5/6 (lft::cf-stream->rational (lft::->cf-stream sum))))))

;;; -------------------------------------------------------------------
;;; Sign-crossing arithmetic: a negative BINARY-EXPRESSION difference
;;; prints with a leading "-", reports negative via MINUSP-style
;;; comparison, and round-trips to the same negative double through
;;; CF-STREAM (binary-expression.lisp + lft-stream.lisp + cf-stream.lisp).
;;; -------------------------------------------------------------------

(test negative-binary-expression-difference-end-to-end
  (let ((difference (- sqrt-two phi)))
    (is (< (lft::lft-stream->double difference) 0))
    (is (search "-0.203820426376799799402" (format nil "~a" difference)))
    (is (approx= (lft::lft-stream->double
                  (lft::cf-stream->lft-stream (lft::->cf-stream difference)))
                 -0.2038204263767998d0))))

;;; -------------------------------------------------------------------
;;; Pythagorean trig identity across BINARY-EXPRESSION's square (via
;;; BILFT) and LFT-ARITH's SIN/COS, on an irrational LFT-STREAM argument
;;; (binary-expression.lisp + lft-arith.lisp).
;;; -------------------------------------------------------------------

(test sin-squared-plus-cos-squared-identity
  (is (approx= (+ (square (sin sqrt-two)) (square (cos sqrt-two))) 1.0d0)))

;;; -------------------------------------------------------------------
;;; Hyperbolic Pythagorean identity: cosh^2(x) - sinh^2(x) = 1
;;; (lft-arith.lisp's hyperbolic functions + binary-expression square).
;;; -------------------------------------------------------------------

(test cosh-squared-minus-sinh-squared-identity
  (is (approx= (- (square (cosh sqrt-two)) (square (sinh sqrt-two))) 1.0d0)))

;;; -------------------------------------------------------------------
;;; Inverse-function round trips through LFT-ARITH's transcendental
;;; functions on LFT-STREAM arguments (log/exp, atan/tan, asin/sin,
;;; acos/cos, atanh/tanh all cross-checked together).
;;; -------------------------------------------------------------------

(test transcendental-inverse-round-trips
  (is (approx= (log (exp sqrt-two)) 1.4142135623730951d0))
  (is (approx= (atan (tan 1)) 1.0d0))
  (is (approx= (square (sqrt sqrt-two)) 1.4142135623730951d0))
  (is (approx= (asin (sin 1/2)) 0.5d0))
  (is (approx= (acos (cos 1/2)) 0.5d0))
  (is (approx= (tanh (atanh 1/2)) 0.5d0)))

;;; -------------------------------------------------------------------
;;; PI-based identities: sin(pi/2) = 1, cos(pi/2) = 0 (lft-arith.lisp's
;;; PI/2 constant -- itself built from the Ramanujan/Chudnovsky-style
;;; BIG-K-STREAM machinery in lft-stream.lisp -- combined with SIN/COS).
;;; -------------------------------------------------------------------

(test pi-over-2-trig-identity
  (is (approx= (sin lft::pi/2) 1.0d0))
  (is (approx= (cos lft::pi/2) 0.0d0)))

;;; -------------------------------------------------------------------
;;; exp(x + 1) = exp(x) * e: BINARY-EXPRESSION's + (on an LFT-STREAM and
;;; an exact integer) combined with LFT-ARITH's EXP, cross-checked
;;; against plain double-float arithmetic.
;;; -------------------------------------------------------------------

(test exp-of-sum-equals-product-of-exps
  (is (approx= (exp (+ sqrt-two 1))
               (* (lft::lft-stream->double (exp sqrt-two)) (exp 1.0d0))
               1d-6)))
