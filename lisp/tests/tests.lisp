;;; -*- Lisp -*-

;;; PACKAGE-BOUNDARY POLICY: this suite deliberately reaches into
;;; LINEAR-FRACTIONAL-TRANSFORMATION's internal (non-exported) symbols
;;; throughout via the `lft::foo` prefix -- see tests/package.lisp for why
;;; LINEAR-FRACTIONAL-TRANSFORMATION-TESTS is its own ASDF system/package
;;; in the first place -- not an oversight to be cleaned up.

(in-package "LFT-TESTS")

(def-suite lft-suite
  :description "Test suite for lft.lisp (linear-fractional-transformation)")

(in-suite lft-suite)

;;; -------------------------------------------------------------------
;;; Low-level comparison / determinant predicates
;;; -------------------------------------------------------------------

(test 2x1v<
  (is-true  (lft::2x1v< 1 2 3 4))   ; 1/2 < 3/4
  (is-false (lft::2x1v< 3 4 1 2))   ; 3/4 < 1/2 is false
  (is-false (lft::2x1v< 1 2 1 2)))  ; equal, not strictly less

(test 2x2m<2x1v
  ;; both 1/1 and 2/1 are less than 5/1
  (is-true  (lft::2x2m<2x1v 1 2 1 1 5 1))
  ;; 2/1 is not less than 1/1, so the conjunction fails
  (is-false (lft::2x2m<2x1v 1 2 1 1 1 1)))

(test 2x2m<2x2m
  (is-true  (lft::2x2m<2x2m 1 2 1 1 5 6 1 1))
  (is-false (lft::2x2m<2x2m 1 2 1 1 1 6 1 1)))

(test unimodular?
  (is-true  (lft::unimodular? 1 0 0 1))    ; identity, det 1
  (is-true  (lft::unimodular? 0 1 1 0))    ; swap, det -1
  (is-true  (lft::unimodular? 2 3 1 2))    ; det 4-3 = 1
  (is-false (lft::unimodular? 2 0 0 1))    ; det 2
  (is-false (lft::unimodular? 2 2 2 2)))   ; det 0 (singular)

(test disjoint?
  ;; range [1/1, 2/1] vs [3/1, 4/1] -- disjoint
  (is-true  (lft::disjoint? 1 2 1 1 3 4 1 1))
  ;; range [1/1, 3/1] vs [2/1, 4/1] -- overlapping
  (is-false (lft::disjoint? 1 3 1 1 2 4 1 1)))

(test 2x2-matrix-range?
  (is-true  (lft::2x2-matrix-range? 1 2 0 1))
  (is-false (lft::2x2-matrix-range? 1 -1 1 1)))

;;; -------------------------------------------------------------------
;;; evaluate-lft / evaluate-bilft
;;; -------------------------------------------------------------------

(test evaluate-lft-ordinary
  (is (= 7/2 (lft::evaluate-lft 2 3 0 2 2 1)))     ; (2*2+3)/(0*2+2*1) = 7/2
  (is (= 4   (lft::evaluate-lft 1 2 0 1 2 1))))    ; x + 2 at x=2

(test evaluate-lft-pole-unimodular
  ;; unimodular matrix (det -1) representing 1/x: pole at x=0 maps to infinity
  (is (eq 'lft:infinity (lft::evaluate-lft 0 1 1 0 0 1))))

(test evaluate-lft-pole-singular
  ;; singular (det 0) matrix: c non-zero, limit is a/c
  (is (= 0 (lft::evaluate-lft 0 0 1 0 1 0)))
  ;; singular, c zero, d non-zero: limit is b/d
  (is (= 3 (lft::evaluate-lft 0 6 0 2 1 0))))

(test evaluate-lft-division-by-zero
  (signals division-by-zero (lft::evaluate-lft 0 0 0 0 1 1)))

(test evaluate-bilft-ordinary
  ;; (xy + x + y + 1)/1 at x=2,y=3 -> (6+2+3+1)/1 = 12
  (is (= 12 (lft::evaluate-bilft 1 1 1 1 0 0 0 1 2 1 3 1))))

;;; -------------------------------------------------------------------
;;; Formatting
;;; -------------------------------------------------------------------

(test format-lft-equation
  (is (string= "2x + 3" (with-output-to-string (s) (lft::format-lft-equation s 2 3 0 1))))
  (is (string= "1 + 0/x" (with-output-to-string (s) (lft::format-lft-equation s 1 0 1 0))))
  (is (string= "2 - 3/x" (with-output-to-string (s) (lft::format-lft-equation s 2 -3 1 0)))))

(test format-bilft-equation
  (is (string= "xy + 2x + 3y + 4"
               (with-output-to-string (s) (lft::format-bilft-equation s 1 2 3 4 0 0 0 1)))))

(test print-object-lft
  (is (string= "#<LFT 2x + 3>" (format nil "~a" (make-lft 2 3 0 1)))))

(test print-object-bilft
  (is (string= "#<BILFT x + y>" (format nil "~a" lft::bilft-add))))

;;; -------------------------------------------------------------------
;;; Construction and coefficient canonicalization
;;; -------------------------------------------------------------------

(defun lft-coefficients (lft)
  (list (slot-value lft 'lft::a) (slot-value lft 'lft::b)
        (slot-value lft 'lft::c) (slot-value lft 'lft::d)))

(test make-lft-with-integers
  (is (equal '(1 2 3 4) (lft-coefficients (make-lft 1 2 3 4)))))

(test make-lft-with-rationals
  ;; coefficients with rational entries are cleared to a common integer form
  (is (equal '(1 2 3 4) (lft-coefficients (make-lft 1/1 2 3 4))))
  (is (equal '(1 4 0 2) (lft-coefficients (make-lft 1/2 2 0 1)))))

(test make-lft-with-floats
  (is (equal '(1 2 3 4) (lft-coefficients (make-lft 1.0 2.0 3.0 4.0)))))

(test make-lft-canonicalizes-sign
  ;; a negative C (or zero C with negative D) flips the sign of every coefficient
  (is (equal '(-1 -2 1 -3) (lft-coefficients (make-lft 1 2 -1 3)))))

(test make-lft-reduces-common-powers-of-two
  ;; non-unimodular, all-even coefficients are halved until no longer all even
  (is (equal '(1 2 1 3) (lft-coefficients (make-lft 2 4 2 6)))))

(test make-lft-skips-reduction-when-unimodular
  ;; determinant +/-1 coefficients are left untouched even if not in lowest terms
  ;; relative to each other (there is nothing further to reduce)
  (is (equal '(2 3 1 2) (lft-coefficients (make-lft 2 3 1 2)))))

(test list->lft
  (is (equal '(1 2 3 4) (lft-coefficients (lft::list->lft '(1 2 3 4))))))

;;; -------------------------------------------------------------------
;;; funcall / evaluation through MAKE-LFT instances
;;; -------------------------------------------------------------------

(test funcall-lft-on-integer-and-rational
  (let ((f (make-lft 1 1 0 1))) ; x + 1
    (is (= 4 (funcall f 3)))
    (is (= 5/2 (funcall f 3/2)))))

(test funcall-lft-on-infinity
  (let ((f (make-lft 0 1 1 0))) ; 1/x
    (is (= 0 (funcall f 'lft:infinity)))))

(test funcall-lft-on-float
  (let ((f (make-lft 1 1 0 1)))
    (is (= 4.0d0 (funcall f 3.0d0)))))

;;; -------------------------------------------------------------------
;;; Inverse / negate / reciprocal / composition
;;; -------------------------------------------------------------------

(test inverse-lft-round-trips
  (let* ((f (make-lft 2 3 1 2))      ; unimodular
         (f-inv (inverse f)))
    (is (= 5 (funcall f-inv (funcall f 5))))
    (is (= 7 (funcall f (funcall f-inv 7))))))

(test inverse-lft-is-tame-adjugate
  ;; The tame inverse is the adjugate [d -b; -c a]; here that is (2 -3 -1 2),
  ;; but canonicalization negates all four coefficients because the
  ;; adjugate's C coefficient (-1) is negative.
  (is (equal '(-2 3 1 -2) (lft-coefficients (inverse (make-lft 2 3 1 2))))))

(test negate-lft
  (let ((f (make-lft 1 1 0 1)))
    (is (= -4 (funcall (negate f) 3)))))

(test reciprocal-lft
  (let ((f (make-lft 1 1 0 1)))
    (is (= 1/5 (funcall (reciprocal f) 4)))))

(test compose2-lft
  (let ((f (make-lft 1 1 0 1))   ; x + 1
        (g (make-lft 2 0 0 1)))  ; 2x
    (is (= 7 (funcall (compose2 f g) 3)))     ; f(g(3)) = 2*3+1 = 7
    (is (= 7 (funcall (compose f g) 3)))))    ; COMPOSE is an alias via FUNCTION

;;; -------------------------------------------------------------------
;;; Generic arithmetic built on LFT (ADD2, MULTIPLY2, SQUARE, CUBE, EXPT)
;;; -------------------------------------------------------------------

(test square-and-cube-lft
  (let ((f (make-lft 1 1 0 1))) ; x + 1
    (is (= 5 (funcall (square f) 3)))    ; f(f(3)) = 3+1+1 = 5
    (is (= 6 (funcall (cube f) 3)))))    ; f(f(f(3))) = 6

(test expt-lft
  (let ((f (make-lft 2 0 0 1))) ; 2x
    (is (= 8 (funcall (expt f 3) 1)))        ; (2x)^3 applied as composition, at x=1 -> 8
    ;; GENERIC-ARITHMETIC's EXPT special-cases exponent 0 (via an
    ;; (EQL 0) method with higher precedence than LFT.lisp's own
    ;; ((base lft) (exponent integer)) method) and always returns the
    ;; literal 1, regardless of base type.
    (is (= 1 (expt f 0)))
    (is (= 4 (funcall (expt f -1) 8)))))     ; inverse of 2x is x/2

(test add2-lft-rational
  (let ((f (make-lft 1 1 0 1))) ; x + 1
    (is (= 9 (funcall (add2 f 5) 3)))
    (is (= 9 (funcall (add2 5 f) 3)))))

(test multiply2-lft-rational
  (let ((f (make-lft 1 1 0 1))) ; x + 1
    (is (= 20 (funcall (multiply2 f 5) 3)))
    (is (= 20 (funcall (multiply2 5 f) 3)))))

(test add2-lft-lft-yields-bilft
  ;; adding two LFTs of independent variables yields a BILFT of (x, y)
  (let ((result (add2 (make-lft 1 1 0 1) (make-lft 2 0 0 1)))) ; (x+1) + (2y)
    (is (typep result 'lft:bilft))
    (is (= 12 (funcall result 3 4)))))

(test multiply2-lft-lft-yields-bilft
  (let ((result (multiply2 (make-lft 1 1 0 1) (make-lft 2 0 0 1)))) ; (x+1) * (2y)
    (is (= 32 (funcall result 3 4)))))

;;; -------------------------------------------------------------------
;;; Composing LFTs with BILFTs
;;; -------------------------------------------------------------------

(test compose-lft-bilft
  (let* ((g (make-lft 1 1 0 1))               ; x + 1
         (bi (make-bilft 1 0 0 0 0 0 0 1)))   ; x*y
    (is (= 13 (funcall (lft::compose-lft-bilft g bi) 3 4)))
    (is (= 13 (funcall (lft::funcall-lft g bi) 3 4)))))

(test compose-bilft-lft-x-and-y
  (let ((bi (make-bilft 1 0 0 0 0 0 0 1))   ; x*y
        (h (make-lft 2 0 0 1)))             ; 2x
    (is (= 24 (funcall (lft::compose-bilft-lft-x bi h) 3 4)))
    (is (= 24 (funcall (lft::compose-bilft-lft-y bi h) 3 4)))))

(test transpose-bilft
  (let ((tb (lft::transpose-bilft (make-bilft 1 2 3 4 5 6 7 8))))
    (is (equal '(1 3 2 4 5 7 6 8)
               (list (slot-value tb 'lft::a) (slot-value tb 'lft::b)
                     (slot-value tb 'lft::c) (slot-value tb 'lft::d)
                     (slot-value tb 'lft::e) (slot-value tb 'lft::f)
                     (slot-value tb 'lft::g) (slot-value tb 'lft::h))))))

;;; -------------------------------------------------------------------
;;; Range / pole predicates
;;; -------------------------------------------------------------------

(test range?
  (is-true  (lft::range? (make-lft 1 2 0 1)))
  (is-false (lft::range? (make-lft 1 -1 1 1))))

(test pole-negative-p
  (is-true  (lft::pole-negative-p (make-lft 1 0 0 1)))
  (is-true  (lft::pole-negative-p (make-lft 1 1 1 0))))

;;; -------------------------------------------------------------------
;;; Sign / comparison predicates on an LFT's range for non-negative x
;;; -------------------------------------------------------------------

(test lft-zero-p
  (is (eq :always    (lft::lft-zero-p (make-lft 0 0 1 0)
                                       (lambda () :always) (lambda () :never) (lambda () :sometimes))))
  (is (eq :never      (lft::lft-zero-p (make-lft 1 1 0 1)
                                        (lambda () :always) (lambda () :never) (lambda () :sometimes)))))

(test lft-minus-p
  (is (eq :always (lft::lft-minus-p (make-lft -1 -1 0 1)
                                     (lambda () :always) (lambda () :never) (lambda () :unknown)))))

(test lft-plus-p
  (is (eq :always (lft::lft-plus-p (make-lft 1 1 0 1)
                                    (lambda () :always) (lambda () :never) (lambda () :unknown)))))

(test lft-non-negative-p
  (is (eq :nn (lft::lft-non-negative-p (make-lft 0 1 0 1)
                                        (lambda () :nn) (lambda () :nnn) (lambda () :unknown)))))

(test lft-less-than-rat
  (is (eq :lt (lft::lft-less-than-rat (make-lft 0 1 1 10) 1
                                       (lambda () :lt) (lambda () :not-lt) (lambda () :unknown)))))

(test lft-greater-than-rat
  (is (eq :not-gt (lft::lft-greater-than-rat (make-lft 0 1 1 10) 1
                                              (lambda () :gt) (lambda () :not-gt) (lambda () :unknown)))))

;;; -------------------------------------------------------------------
;;; Rational-combinator helper LFTs
;;; -------------------------------------------------------------------

(test lft-add-rat
  (is (= 8 (funcall (lft::lft-add-rat 5) 3))))

(test lft-subtract-rat
  (is (= -2 (funcall (lft::lft-subtract-rat 5) 3))))

(test lft-multiply-by-rat
  (is (= 15 (funcall (lft::lft-multiply-by-rat 5) 3))))

(test lft-divide-by-rat
  (is (= 3/5 (funcall (lft::lft-divide-by-rat 5) 3))))

;;; -------------------------------------------------------------------
;;; Digit emission / numeric coercion
;;; -------------------------------------------------------------------

(test lft-truncate-success
  ;; (3x+6)/(x+2): floor(3/1) = floor(6/2) = 3, so truncation succeeds with digit 3
  (lft::lft-truncate (make-lft 3 6 1 2)
                      (lambda (digit remainder) (declare (ignore remainder)) (is (= 3 digit)))
                      (lambda () (fail "expected lft-truncate to succeed"))))

(test lft-truncate-failure
  ;; (3x+1)/(x+1): floor(3/1)=3, floor(1/1)=1, disagree, so truncation fails
  (is (eq :fail (lft::lft-truncate (make-lft 3 1 1 1)
                                    (lambda (digit remainder) (declare (ignore digit remainder)) :ok)
                                    (lambda () :fail)))))

(test lft-truncate-zero-a-zero-c-constant-regression
  ;; Same exact-constant special case as LFT->DOUBLE: A=C=0 means
  ;; (a*x+b)/(c*x+d) = b/d for every x, which TRUNCATE(A,C) could not
  ;; even compute (division by zero) before this fix. See
  ;; LFT->DOUBLE-ZERO-A-ZERO-C-CONSTANT-REGRESSION below for the full
  ;; explanation.
  (lft::lft-truncate (make-lft 0 5 0 2)
                      (lambda (digit remainder)
                        (is (= 2 digit))
                        (is (= 1/2 (lft::lft->double remainder #'identity (lambda () (fail "expected remainder to resolve"))))))
                      (lambda () (fail "expected lft-truncate to succeed")))
  (is (eq :fail (lft::lft-truncate (make-lft 0 1 0 0)
                                    (lambda (digit remainder) (declare (ignore digit remainder)) :ok)
                                    (lambda () :fail)))))

(test lft->single-success
  (is (= 3.0 (lft::lft->single (make-lft 3 6 1 2) #'identity (lambda () :fail)))))

(test lft->single-failure-without-positive-c-d
  (is (eq :fail (lft::lft->single (make-lft 1 1 0 1) #'identity (lambda () :fail)))))

(test lft->single-zero-a-zero-c-constant-regression
  ;; Same exact-constant special case as LFT->DOUBLE (see
  ;; LFT->DOUBLE-ZERO-A-ZERO-C-CONSTANT-REGRESSION below for the full
  ;; explanation of why A=C=0 must resolve directly to B/D).
  (is (= 2.5 (lft::lft->single (make-lft 0 5 0 2) #'identity (lambda () :fail))))
  (is (eq :fail (lft::lft->single (make-lft 0 1 0 0) #'identity (lambda () :fail)))))

(test lft->double-success
  (is (= 3.0d0 (lft::lft->double (make-lft 3 6 1 2) #'identity (lambda () :fail)))))

(test lft->double-zero-a-zero-c-constant-regression
  ;; When A and C are both zero, (a*x+b)/(c*x+d) = b/d for every x -- the
  ;; LFT is already an exact constant and should resolve even though C is
  ;; not strictly positive. Before this fix, LFT->DOUBLE rejected this
  ;; case and asked for another refinement; REFINE-LFT-STREAM has no more
  ;; digits to give for a terminating stream, so the degenerate LFT kept
  ;; collapsing toward all-zero coefficients and LFT-STREAM->DOUBLE
  ;; looped forever on exact-rational streams (e.g. (->lft-stream 1)).
  (is (= 2.5d0 (lft::lft->double (make-lft 0 5 0 2) #'identity (lambda () :fail))))
  ;; A = C = D = 0 is genuinely unrecoverable (0/0 for every x) and must
  ;; still fail.
  (is (eq :fail (lft::lft->double (make-lft 0 1 0 0) #'identity (lambda () :fail)))))

(test lft-stream->double-exact-rational-regression
  ;; LFT-STREAM->DOUBLE used to hang indefinitely when converting a
  ;; stream built from an exact rational (e.g. via ->LFT-STREAM, or
  ;; produced internally by functions like SQRT/ACOS on exact inputs)
  ;; because REFINE-LFT-STREAM has no further digits to refine a
  ;; terminating stream with. See LFT->DOUBLE-ZERO-A-ZERO-C-CONSTANT-
  ;; REGRESSION above for the underlying fix.
  (is (= 1.0d0 (lft::lft-stream->double (lft::->lft-stream 1))))
  (is (= 0.0d0 (lft::lft-stream->double (lft::->lft-stream 0))))
  (is (= 1.5d0 (lft::lft-stream->double (lft::->lft-stream 3/2))))
  (is (= -5.0d0 (lft::lft-stream->double (lft::->lft-stream -5)))))

;;; -------------------------------------------------------------------
;;; BILFT sign / zero / truncation predicates
;;; -------------------------------------------------------------------

(test poles-negative-p
  (is-false (lft::poles-negative-p lft::bilft-add))
  (is-true  (lft::poles-negative-p (make-bilft 0 0 0 0 1 0 0 1))))

(test bilft-zero-p
  (is (eq :zero (lft::bilft-zero-p (make-bilft 0 0 0 0 1 0 0 1)
                                    (lambda () :zero) (lambda () :not-zero) (lambda () :unknown))))
  (is (eq :not-zero (lft::bilft-zero-p (make-bilft 1 1 1 1 1 1 1 1)
                                        (lambda () :zero) (lambda () :not-zero) (lambda () :unknown)))))

(test bilft-minus-p
  (is (eq :plus (lft::bilft-minus-p (make-bilft 1 1 1 1 1 1 1 1)
                                     (lambda () :minus) (lambda () :plus) (lambda () :unknown))))
  (is (eq :minus (lft::bilft-minus-p (make-bilft -1 -1 -1 -1 1 1 1 1)
                                      (lambda () :minus) (lambda () :plus) (lambda () :unknown)))))

(test bilft-truncate
  (is (= 2 (lft::bilft-truncate (make-bilft 2 2 2 2 1 1 1 1)
                                 (lambda (digit remainder) (declare (ignore remainder)) digit)
                                 (lambda () :fail)))))

;; To run all tests:
;; (fiveam:run! 'lft-tests:lft-suite)
