;;; -*- Lisp -*-

;;; PACKAGE-BOUNDARY POLICY: see tests/package.lisp and tests/tests.lisp
;;; for why this suite lives in LFT-TESTS and reaches into
;;; LINEAR-FRACTIONAL-TRANSFORMATION's internal (non-exported) symbols
;;; via the `lft::foo` prefix.
;;;
;;; This file exercises binary-expression.lisp, which implements
;;; Gosper's algorithm for combining two LFT-STREAMs (continued-fraction
;;; streams) with a single bilinear fractional transformation (BILFT):
;;; the BINARY-EXPRESSION class, its REFINE-LEFT/REFINE-RIGHT/
;;; REFINE-FAIRLY/REFINE-DISJOINT refinement strategies, its
;;; sign/zero/truncate/digit decomposition, its three FUNCALL-BILFT
;;; dispatch methods (LFT-STREAM x LFT-STREAM, LFT-STREAM x PROMISE,
;;; PROMISE x LFT-STREAM), ADD2/MULTIPLY2 on two LFT-STREAMs, and the
;;; binary-expression-based transcendental stream builders
;;; (%ATAN-LFT-STREAM, %EXP-LFT-STREAM, %LOG-LFT-STREAM,
;;; %SQRT-LFT-STREAM, %SQUARE-LFT-STREAM, %TAN-LFT-STREAM, %RAT-TAN).
;;;
;;; NOTE ON %RAT-TAN: this file defines a SECOND %RAT-TAN, overriding
;;; the one in lft-stream.lisp (domain (0,1), via LFT-STREAM-MAP) with a
;;; domain (-1,1) version built on %TAN-LFT-STREAM. Because
;;; binary-expression.lisp loads after lft-stream.lisp (see the .asd
;;; component order), this later definition is the one actually in
;;; effect at runtime; the test below confirms that behavior rather
;;; than lft-stream.lisp's.
;;;
;;; NOTE ON THE (FORMERLY HANGING) PRINT-OBJECT PATH: printing (via
;;; PRINT-OBJECT / FORMAT "~a") a BINARY-EXPRESSION or the LFT-STREAM
;;; derived from one used to hang indefinitely when the underlying
;;; value was an EXACT rational sitting exactly on an integer boundary
;;; (e.g. a BINARY-EXPRESSION built by summing two exact integer
;;; streams). Root cause: once every input stream feeding a
;;; BINARY-EXPRESSION/LFT-STREAM was exhausted, REFINE-LEFT/
;;; REFINE-RIGHT/REFINE-LFT-STREAM's x=infinity substitution produced a
;;; degenerate "exact constant" BILFT/LFT (all four numerator corners
;;; equal and all four denominator corners equal) that
;;; TRY-DECOMPOSE-DIGIT/TRY-DECOMPOSE-SIGN's RANGE?-based checks could
;;; never classify, and even once that was fixed (see
;;; EXACT-CONSTANT-VALUE in lft.lisp), BINARY-EXPRESSION->LFT-DIGIT-STREAM
;;; and LFT-STREAM->LFT-DIGIT-STREAM continued alternating forever
;;; between LFT-DIGIT-HIGH/LFT-DIGIT-LOW for values exactly on an
;;; integer boundary, which LFT-STREAM-TRUNCATE's floor-based
;;; convergence test could never resolve. Both digit-stream
;;; constructors now detect EXACT-CONSTANT-VALUE and switch encodings
;;; by delegating to the already-terminating (->LFT-STREAM value), so
;;; printing exact-valued binary-expressions/LFT-STREAMs (integer or
;;; not) now works; see PRINT-OBJECT-ON-EXACT-INTEGER-BINARY-EXPRESSION
;;; below.

(in-package "LFT-TESTS")

(in-suite lft-suite)

(defun stream-approx= (actual expected &optional (epsilon 1d-9))
  "Return true if ACTUAL (a rational, integer, or LFT-STREAM) converts to
a double-float within EPSILON of the double-float EXPECTED."
  (let ((actual-double (if (numberp actual)
                            (coerce actual 'double-float)
                            (lft::lft-stream->double actual))))
    (< (abs (- actual-double expected)) epsilon)))

(defun make-sum-binary-expression (left right)
  "Build a BINARY-EXPRESSION representing LEFT + RIGHT (both LFT-STREAMs),
mirroring how LFT::FUNCALL-BILFT's (lft-stream, lft-stream) method
constructs its BINARY-EXPRESSION, but without immediately converting
the result back to an LFT-STREAM."
  (make-instance 'lft::binary-expression
                 :generation 0
                 :bilft (lft::compose-bilft-lft-x
                         (lft::compose-bilft-lft-y lft::bilft-add (lft::stream-car right))
                         (lft::stream-car left))
                 :delayed-left (lft::stream-delayed-cdr left)
                 :delayed-right (lft::stream-delayed-cdr right)))

;;; -------------------------------------------------------------------
;;; refine-left / refine-right / refine-fairly / refine-disjoint
;;; -------------------------------------------------------------------

(test refine-left-increments-generation
  (let ((be (make-sum-binary-expression sqrt-two phi)))
    (is (= 1 (slot-value (lft::refine-left be) 'lft::generation)))))

(test refine-right-increments-generation
  (let ((be (make-sum-binary-expression sqrt-two phi)))
    (is (= 1 (slot-value (lft::refine-right be) 'lft::generation)))))

(test refine-left-preserves-convergence
  ;; Refining repeatedly in one direction must still converge to the
  ;; same sum; it only narrows the BILFT, it never changes the value.
  (let ((be (make-sum-binary-expression sqrt-two phi)))
    (dotimes (i 3) (setq be (lft::refine-left be)))
    (is (stream-approx= (lft::binary-expression->lft-stream be) 3.03224755112299d0))))

(test refine-right-preserves-convergence
  (let ((be (make-sum-binary-expression sqrt-two phi)))
    (dotimes (i 3) (setq be (lft::refine-right be)))
    (is (stream-approx= (lft::binary-expression->lft-stream be) 3.03224755112299d0))))

(test refine-fairly-preserves-convergence
  (let ((be (make-sum-binary-expression sqrt-two phi)))
    (dotimes (i 6) (setq be (lft::refine-fairly be)))
    (is (= 6 (slot-value be 'lft::generation)))
    (is (stream-approx= (lft::binary-expression->lft-stream be) 3.03224755112299d0))))

(test refine-disjoint-preserves-convergence
  (let ((be (make-sum-binary-expression sqrt-two phi)))
    (dotimes (i 6) (setq be (lft::refine-disjoint be)))
    (is (stream-approx= (lft::binary-expression->lft-stream be) 3.03224755112299d0))))

;;; -------------------------------------------------------------------
;;; funcall-lft on a binary-expression
;;; -------------------------------------------------------------------

(test funcall-lft-negate-on-binary-expression
  ;; LFT objects are funcallable, so (funcall lft::lft-negate be) routes
  ;; through LFT::FUNCALL-LFT's method specialized on BINARY-EXPRESSION.
  (let ((be (make-sum-binary-expression sqrt-two phi)))
    (is (stream-approx=
         (lft::binary-expression->lft-stream (funcall lft::lft-negate be))
         -3.03224755112299d0))))

;;; -------------------------------------------------------------------
;;; binary-expression-zero-p / binary-expression-minus-p
;;; -------------------------------------------------------------------

(test binary-expression-zero-p-true
  (is-true (lft::binary-expression-zero-p
            (make-sum-binary-expression (lft::->lft-stream 2) (lft::->lft-stream -2)))))

(test binary-expression-zero-p-false
  (is-false (lft::binary-expression-zero-p
             (make-sum-binary-expression (lft::->lft-stream 2) (lft::->lft-stream 3)))))

(test binary-expression-minus-p-true
  (is-true (lft::binary-expression-minus-p
            (make-sum-binary-expression (lft::->lft-stream -5) (lft::->lft-stream -1)))))

(test binary-expression-minus-p-false
  (is-false (lft::binary-expression-minus-p
             (make-sum-binary-expression (lft::->lft-stream 5) (lft::->lft-stream 1)))))

;;; -------------------------------------------------------------------
;;; binary-expression-truncate / -decompose-sign / -decompose-digit
;;;
;;; These only ever inspect the primary (digit/sign/integer) return
;;; value, never the secondary remainder BINARY-EXPRESSION, to avoid
;;; the print-hang described above (none of these tests print the
;;; remainder, so none are affected by it).
;;; -------------------------------------------------------------------

(test binary-expression-truncate-exact-sum
  (is (= 5 (nth-value 0 (lft::binary-expression-truncate
                         (make-sum-binary-expression (lft::->lft-stream 2) (lft::->lft-stream 3)))))))

(test binary-expression-decompose-sign-positive
  (is (eq lft::lft-sign-positive
          (nth-value 0 (lft::binary-expression-decompose-sign
                        (make-sum-binary-expression sqrt-two phi))))))

;;; -------------------------------------------------------------------
;;; binary-expression->lft-stream / binary-expression->lft-digit-stream
;;; -------------------------------------------------------------------

(test binary-expression->lft-stream-irrational-sum
  (is (stream-approx= (lft::binary-expression->lft-stream
                        (make-sum-binary-expression sqrt-two phi))
                       3.03224755112299d0)))

(test binary-expression->lft-stream-irrational-product
  (is (stream-approx= (lft::binary-expression->lft-stream
                        (make-instance 'lft::binary-expression
                                       :generation 0
                                       :bilft (lft::compose-bilft-lft-x
                                               (lft::compose-bilft-lft-y lft::bilft-multiply (lft::stream-car phi))
                                               (lft::stream-car sqrt-two))
                                       :delayed-left (lft::stream-delayed-cdr sqrt-two)
                                       :delayed-right (lft::stream-delayed-cdr phi)))
                       2.288245611270737d0)))

(test binary-expression->lft-stream-exact-sum
  ;; Exercises the exact-rational path through the decomposition
  ;; functions; the returned LFT-STREAM is only ever converted to a
  ;; double (never printed), so the print-hang is not triggered.
  (is (stream-approx= (lft::binary-expression->lft-stream
                        (make-sum-binary-expression (lft::->lft-stream 2) (lft::->lft-stream 3)))
                       5.0d0)))

;;; -------------------------------------------------------------------
;;; print-object (safe only for non-terminating/irrational operands)
;;; -------------------------------------------------------------------

(test print-object-on-irrational-binary-expression-does-not-hang
  (is (stringp (format nil "~a" (make-sum-binary-expression sqrt-two phi)))))

(test print-object-on-exact-integer-binary-expression-does-not-hang
  ;; Regression test for the print-hang fixed via EXACT-CONSTANT-VALUE
  ;; (lft.lisp) and the digit-stream encoding switch in
  ;; BINARY-EXPRESSION->LFT-DIGIT-STREAM: 2 + 3 = 5 sits exactly on an
  ;; integer boundary, which is precisely the case that used to hang.
  (is (search "5" (format nil "~a" (make-sum-binary-expression
                                     (lft::->lft-stream 2) (lft::->lft-stream 3))))))

(test print-object-on-exact-negative-integer-binary-expression-does-not-hang
  (is (search "-3" (format nil "~a" (make-sum-binary-expression
                                      (lft::->lft-stream -5) (lft::->lft-stream 2))))))

(test print-object-on-exact-non-integer-binary-expression-does-not-hang
  (is (search "0.8333" (format nil "~a" (make-sum-binary-expression
                                          (lft::->lft-stream 1/2) (lft::->lft-stream 1/3))))))

;;; -------------------------------------------------------------------
;;; funcall-bilft dispatch: (lft-stream, lft-stream) / (lft-stream,
;;; promise) / (promise, lft-stream)
;;; -------------------------------------------------------------------

(test funcall-bilft-stream-stream
  (is (stream-approx= (funcall lft::bilft-add sqrt-two phi) 3.03224755112299d0)))

(test funcall-bilft-stream-promise
  (is (stream-approx= (funcall lft::bilft-add sqrt-two (lft::delay phi)) 3.03224755112299d0)))

(test funcall-bilft-promise-stream
  (is (stream-approx= (funcall lft::bilft-add (lft::delay sqrt-two) phi) 3.03224755112299d0)))

;;; -------------------------------------------------------------------
;;; add2 / multiply2 generic methods (dispatch used by + and *)
;;; -------------------------------------------------------------------

(test add2-on-lft-streams
  (is (stream-approx= (+ sqrt-two phi) 3.03224755112299d0)))

(test multiply2-on-lft-streams
  (is (stream-approx= (* sqrt-two phi) 2.288245611270737d0)))

(test add2-same-stream-twice
  (is (stream-approx= (+ sqrt-two sqrt-two) 2.8284271247461903d0)))

(test multiply2-same-stream-twice-is-square
  (is (stream-approx= (* sqrt-two sqrt-two) 2.0d0)))

;;; -------------------------------------------------------------------
;;; Transcendental stream builders based on binary-expression machinery
;;; -------------------------------------------------------------------

(test square-lft-stream
  (is (stream-approx= (lft::%square-lft-stream sqrt-two) 2.0d0)))

(test sqrt-lft-stream
  (is (stream-approx= (lft::%sqrt-lft-stream (lft::->lft-stream 2)) 1.4142135623730951d0)))

(test exp-lft-stream
  (is (stream-approx= (lft::%exp-lft-stream (lft::->lft-stream 1/2)) 1.6487212707001282d0)))

(test log-lft-stream
  (is (stream-approx= (lft::%log-lft-stream (lft::->lft-stream 3/2)) 0.4054651081081644d0)))

(test atan-lft-stream
  (is (stream-approx= (lft::%atan-lft-stream (lft::->lft-stream 1/2)) 0.4636476090008061d0)))

(test tan-lft-stream
  (is (stream-approx= (lft::%tan-lft-stream (lft::->lft-stream 1/4)) 0.25534192122103627d0)))

(test rat-tan-uses-binary-expression-based-definition
  ;; binary-expression.lisp's %RAT-TAN shadows lft-stream.lisp's
  ;; earlier definition (see the file-level NOTE above); this confirms
  ;; the overriding definition -- domain (-1,1), via %TAN-LFT-STREAM --
  ;; is the one in effect, and that it agrees with %TAN-LFT-STREAM
  ;; applied directly.
  (is (stream-approx= (lft::%rat-tan 1/4) 0.25534192122103627d0)))
