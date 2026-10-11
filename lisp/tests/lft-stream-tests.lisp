;;; -*- Lisp -*-

;;; PACKAGE-BOUNDARY POLICY: see tests/package.lisp and tests/tests.lisp
;;; for why this suite lives in LFT-TESTS and reaches into
;;; LINEAR-FRACTIONAL-TRANSFORMATION's internal (non-exported) symbols
;;; via the `lft::foo` prefix.
;;;
;;; This file exercises lft-stream.lisp, which builds the infinite
;;; continued-fraction-stream machinery (LFT-STREAM objects) on top of
;;; lft.lisp's single LFT workhorse predicates/conversions: comparison
;;; against rationals, sign/zero decomposition, digit emission, and
;;; conversion back to double-float, plus the low-level rational-domain
;;; stream generators (%EXP-RAT, %LOG-RAT, %RAT-POW, %RAT-ATAN,
;;; %RAT-TAN) that lft-arith.lisp builds its public transcendental
;;; functions on top of.
;;;
;;; NOTE ON CPS RECEIVER SIGNATURES: LFT-STREAM-GREATER-THAN-RAT,
;;; LFT-STREAM-LESS-THAN-RAT, LFT-STREAM-MINUS-P, LFT-STREAM-PLUS-P, and
;;; LFT-STREAM-NON-NEGATIVE-P all call their RECEIVE-TRUE/RECEIVE-FALSE
;;; (or IF-ALWAYS-.../IF-NEVER-...) continuations with the (possibly
;;; refined) LFT-STREAM itself, not a plain boolean -- so tests use
;;; receiver lambdas that discard their argument and return a sentinel
;;; keyword to record which branch fired.
;;;
;;; NOTE ON NAMED CONSTANTS: SQRT-TWO and PHI (both exported, used bare
;;; here) are CIRCULAR-LFT-STREAM instances -- infinite, non-terminating
;;; streams -- making them convenient, always-irrational fixtures that
;;; never hit the exact-rational-stream termination edge cases covered
;;; separately in tests.lisp's LFT-STREAM->DOUBLE-EXACT-RATIONAL-
;;; REGRESSION test.

(in-package "LFT-TESTS")

(in-suite lft-suite)

(defun stream-approx= (actual expected &optional (epsilon 1d-9))
  "Return true if ACTUAL (a rational, integer, or LFT-STREAM) converts to
a double-float within EPSILON of the double-float EXPECTED."
  (let ((actual-double (if (numberp actual)
                            (coerce actual 'double-float)
                            (lft::lft-stream->double actual))))
    (< (abs (- actual-double expected)) epsilon)))

;;; -------------------------------------------------------------------
;;; Named circular-stream constants
;;; -------------------------------------------------------------------

(test circular-lft-stream-constants
  (is (stream-approx= sqrt-two 1.4142135623730951d0))
  (is (stream-approx= phi 1.618033988749895d0)))

;;; -------------------------------------------------------------------
;;; lft-stream-greater-than-rat / lft-stream-less-than-rat
;;; -------------------------------------------------------------------

(test lft-stream-greater-than-rat
  (is (eq :true (lft::lft-stream-greater-than-rat
                 sqrt-two 1
                 (lambda (s) (declare (ignore s)) :true)
                 (lambda (s) (declare (ignore s)) :false))))
  (is (eq :false (lft::lft-stream-greater-than-rat
                  sqrt-two 2
                  (lambda (s) (declare (ignore s)) :true)
                  (lambda (s) (declare (ignore s)) :false)))))

(test lft-stream-less-than-rat
  (is (eq :true (lft::lft-stream-less-than-rat
                 sqrt-two 2
                 (lambda (s) (declare (ignore s)) :true)
                 (lambda (s) (declare (ignore s)) :false))))
  (is (eq :false (lft::lft-stream-less-than-rat
                  sqrt-two 1
                  (lambda (s) (declare (ignore s)) :true)
                  (lambda (s) (declare (ignore s)) :false)))))

;;; -------------------------------------------------------------------
;;; lft-stream-minus-p / lft-stream-plus-p / lft-stream-non-negative-p
;;; -------------------------------------------------------------------

(test lft-stream-plus-p
  (is (eq :always (lft::lft-stream-plus-p
                   sqrt-two
                   (lambda (s) (declare (ignore s)) :always)
                   (lambda (s) (declare (ignore s)) :never)))))

(test lft-stream-minus-p
  (is (eq :always (lft::lft-stream-minus-p
                   (negate sqrt-two)
                   (lambda (s) (declare (ignore s)) :always)
                   (lambda (s) (declare (ignore s)) :never)))))

(test lft-stream-non-negative-p
  (is (eq :yes (lft::lft-stream-non-negative-p
                sqrt-two
                (lambda (s) (declare (ignore s)) :yes)
                (lambda (s) (declare (ignore s)) :no))))
  (is (eq :no (lft::lft-stream-non-negative-p
               (negate sqrt-two)
               (lambda (s) (declare (ignore s)) :yes)
               (lambda (s) (declare (ignore s)) :no)))))

;;; -------------------------------------------------------------------
;;; lft-stream-zero-p / zerop
;;; -------------------------------------------------------------------

(test lft-stream-zero-p
  (is (eq t (lft::lft-stream-zero-p (lft::->lft-stream 0))))
  (is (eq nil (lft::lft-stream-zero-p sqrt-two))))

(test zerop-lft-stream
  (is-true (zerop (lft::->lft-stream 0)))
  (is-false (zerop sqrt-two)))

;;; -------------------------------------------------------------------
;;; Arithmetic on lft-stream via funcall-lft / add2 / multiply2 /
;;; negate / reciprocal
;;; -------------------------------------------------------------------

(test lft-stream-arithmetic
  (is (stream-approx= (+ 1 sqrt-two) 2.414213562373095d0))
  (is (stream-approx= (+ sqrt-two 1) 2.414213562373095d0))
  (is (stream-approx= (* 2 sqrt-two) 2.8284271247461903d0))
  (is (stream-approx= (* sqrt-two 2) 2.8284271247461903d0))
  (is (stream-approx= (negate sqrt-two) -1.4142135623730951d0))
  (is (stream-approx= (reciprocal sqrt-two) 0.7071067811865476d0)))

;;; -------------------------------------------------------------------
;;; lft-stream-decompose-digit / lft-stream-decompose-sign /
;;; canonicalize-lft-stream
;;; -------------------------------------------------------------------

(test lft-stream-decompose-digit
  (multiple-value-bind (digit remainder) (lft::lft-stream-decompose-digit sqrt-two)
    ;; The first continued-fraction digit of sqrt(2) is the LFT 2x+1 (over
    ;; 1), and the remainder stream still converges to the same value
    ;; sqrt(2) decomposes into (the digit undoes to reconstruct it).
    (is (= 2 (slot-value digit 'lft::a)))
    (is (= 1 (slot-value digit 'lft::b)))
    (is (= 0 (slot-value digit 'lft::c)))
    (is (= 1 (slot-value digit 'lft::d)))
    (is (stream-approx= remainder 0.20710678118654757d0))))

(test lft-stream-decompose-sign
  ;; A positive stream decomposes with the identity sign (LFT-SIGN-
  ;; POSITIVE, x) and an unchanged remainder.
  (multiple-value-bind (sign remainder) (lft::lft-stream-decompose-sign sqrt-two)
    (is (eq sign lft::lft-sign-positive))
    (is (stream-approx= remainder 1.4142135623730951d0)))
  ;; A negative stream decomposes with LFT-SIGN-NEGATIVE and a remainder
  ;; equal to the reciprocal of its magnitude.
  (multiple-value-bind (sign remainder) (lft::lft-stream-decompose-sign (negate sqrt-two))
    (is (eq sign lft::lft-sign-negative))
    (is (stream-approx= remainder 0.7071067811865476d0)))
  ;; Exactly zero still decomposes successfully (as positive-signed zero).
  (multiple-value-bind (sign remainder) (lft::lft-stream-decompose-sign (lft::->lft-stream 0))
    (is (eq sign lft::lft-sign-positive))
    (is (stream-approx= remainder 0.0d0))))

(test canonicalize-lft-stream
  ;; CANONICALIZE-LFT-STREAM's head is the decomposed sign LFT.
  (is (eq lft::lft-sign-positive (lft::stream-car (lft::canonicalize-lft-stream sqrt-two)))))

(test print-object-on-exact-integer-lft-stream-does-not-hang
  ;; Regression test: printing a plain exact-integer-valued LFT-STREAM
  ;; used to hang because LFT-STREAM->LFT-DIGIT-STREAM alternated
  ;; forever between LFT-DIGIT-HIGH/LFT-DIGIT-LOW once the stream's
  ;; final term collapsed to an EXACT-CONSTANT-VALUE sitting exactly on
  ;; an integer boundary. Fixed by detecting EXACT-CONSTANT-VALUE and
  ;; switching to the already-terminating (->LFT-STREAM value) encoding.
  (is (search "7" (format nil "~a" (lft::->lft-stream 7))))
  (is (search "-3" (format nil "~a" (lft::->lft-stream -3))))
  (is (search "1.25" (format nil "~a" (lft::->lft-stream 5/4)))))

;;; -------------------------------------------------------------------
;;; ->lft-stream (all numeric-type methods)
;;; -------------------------------------------------------------------

(test ->lft-stream-numeric-types
  (is (stream-approx= (lft::->lft-stream 3) 3.0d0))
  (is (stream-approx= (lft::->lft-stream 3/2) 1.5d0))
  (is (stream-approx= (lft::->lft-stream 1.5f0) 1.5d0))
  (is (stream-approx= (lft::->lft-stream 2.25d0) 2.25d0)))

;;; -------------------------------------------------------------------
;;; Low-level rational-domain stream generators used by lft-arith.lisp
;;; -------------------------------------------------------------------

(test exp-rat
  (is (stream-approx= (lft::%exp-rat 1) 2.718281828459045d0)))

(test log-rat
  (is (stream-approx= (lft::%log-rat 3/2) 0.4054651081081644d0)))

(test rat-pow
  (is (stream-approx= (lft::%rat-pow 2 1/2) 1.4142135623730951d0)))

(test rat-atan
  (is (stream-approx= (lft::%rat-atan 1/2) 0.4636476090008061d0)))

(test rat-tan
  (is (stream-approx= (lft::%rat-tan 1/2) 0.5463024898437905d0)))

;;; -------------------------------------------------------------------
;;; Exported stream-building primitives / constants
;;; -------------------------------------------------------------------

(test euler-gompertz-constant
  (is (stream-approx= euler-gompertz 0.5963473623231941d0)))

(test big-k-stream-sanity
  ;; BIG-K-STREAM (cons-stream 1 ...) (cons-stream 1 ...) builds the LFT
  ;; 0 1 / 1 1 (i.e. 1/(1+...)) at its head -- a basic sanity check that
  ;; it produces a well-formed, convertible LFT-STREAM rather than
  ;; testing a specific continued fraction's convergence target.
  (let ((s (lft::big-k-stream (lft::ones) (lft::ones))))
    (is (typep s 'lft::lft-stream))
    (is (numberp (lft::lft-stream-truncate s)))))
