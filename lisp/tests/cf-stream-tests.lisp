;;; -*- Lisp -*-

;;; PACKAGE-BOUNDARY POLICY: see tests/package.lisp and tests/tests.lisp
;;; for why this suite lives in LFT-TESTS and reaches into
;;; LINEAR-FRACTIONAL-TRANSFORMATION's internal (non-exported) symbols
;;; via the `lft::foo` prefix.
;;;
;;; This file exercises cf-stream.lisp, which builds CF-STREAM objects --
;;; infinite (or finite, rational-terminated) streams of plain integer
;;; continued-fraction coefficients -- as an alternative, more
;;; traditional representation alongside lft-stream.lisp's LFT-STREAM
;;; (streams of 1x2 LFT digits). It covers conversion in both directions
;;; between CF-STREAM and LFT-STREAM, conversion to/from CL rationals and
;;; coefficient lists, canonicalization, truncation, RECIPROCAL/NEGATE,
;;; prefix/suffix comparison, and flattening a stream of CF-STREAMs.
;;;
;;; NOTE ON BUGS FOUND AND FIXED WHILE WRITING THESE TESTS: writing this
;;; suite surfaced several real bugs in cf-stream.lisp, all now fixed
;;; (see the corresponding comments in cf-stream.lisp for detail):
;;;   - CF-STREAM->RATIONAL recursed one level too far into a stream's
;;;     empty tail, trying to take the RECIPROCAL of the 'INFINITY
;;;     sentinel for the trailing term of every finite (rational)
;;;     cf-stream, rather than treating the last coefficient as a base
;;;     case.
;;;   - NEGATE (on CF-STREAM) recursed into THE-EMPTY-STREAM (NIL) at the
;;;     end of a finite stream; since NEGATE is specialized on the
;;;     CF-STREAM class and NIL is not a CF-STREAM, this failed with "no
;;;     applicable method" rather than reaching NEGATE's own
;;;     empty-stream check. This made CF-STREAM->LFT-STREAM (which calls
;;;     NEGATE on any cf-stream with a negative leading coefficient) fail
;;;     for every negative rational/cf-stream once fully forced.
;;;   - LFT-STREAM->CF-STREAM called LFT-STREAM-ZERO-P on the fractional
;;;     part returned by LFT-STREAM-TRUNCATE without checking whether
;;;     that fractional part was already THE-EMPTY-STREAM (which
;;;     LFT-STREAM-TRUNCATE returns, instead of an LFT-STREAM, for an
;;;     exact integer with no remainder), causing STREAM-CAR to be
;;;     called on NIL.
;;;   - CF-STREAM-TRUNCATE-DELAYED/CF-STREAM-TRUNCATE did not guard
;;;     against forcing past the actual end of a finite cf-stream.
;;;   - CF-STREAM-FLATTEN-APPEND passed its arguments to
;;;     STREAM-FOLD-RIGHT-DELAYED in the wrong order (STREAM and FINAL
;;;     swapped relative to that function's (f stream final) signature,
;;;     and relative to the correct, pre-existing STREAM-FLATTEN-APPEND
;;;     in the STREAM package it was modeled on), which made it return
;;;     its input unflattened instead of appending it.
;;;   - ->CF-STREAM (the generic numeric-type dispatch function) had no
;;;     method for plain INTEGER, only RATIO/SINGLE-FLOAT/DOUBLE-FLOAT/
;;;     LFT-STREAM, unlike the analogous ->LFT-STREAM which does cover
;;;     INTEGER.

(in-package "LFT-TESTS")

(in-suite lft-suite)

(defun cf-stream->list (cf-stream &optional (limit 20))
  "Return the first LIMIT coefficients of CF-STREAM as an ordinary list
(fewer, if CF-STREAM terminates before LIMIT elements)."
  (if (or (zerop limit) (lft::empty-stream? cf-stream))
      '()
      (cons (lft::stream-car cf-stream)
            (cf-stream->list (lft::stream-cdr cf-stream) (- limit 1)))))

;;; -------------------------------------------------------------------
;;; coefficients->cf-stream / rational->cf-stream / ->cf-stream
;;; -------------------------------------------------------------------

(test coefficients->cf-stream
  (is (equal '(2 3 1 4) (cf-stream->list (lft::coefficients->cf-stream '(2 3 1 4)))))
  (is (lft::empty-stream? (lft::coefficients->cf-stream '()))))

(test rational->cf-stream
  (is (equal '(3) (cf-stream->list (lft::rational->cf-stream 3))))
  (is (equal '(-3) (cf-stream->list (lft::rational->cf-stream -3))))
  (is (equal '(1 2) (cf-stream->list (lft::rational->cf-stream 3/2))))
  (is (equal '(0 3) (cf-stream->list (lft::rational->cf-stream 1/3))))
  (is (equal '(3 7) (cf-stream->list (lft::rational->cf-stream 22/7))))
  (is (equal '(-3 -7) (cf-stream->list (lft::rational->cf-stream -22/7)))))

(test ->cf-stream-numeric-types
  (is (equal '(3) (cf-stream->list (lft::->cf-stream 3))))
  (is (equal '(1 2) (cf-stream->list (lft::->cf-stream 3/2))))
  (is (equal '(1 2) (cf-stream->list (lft::->cf-stream 1.5f0))))
  (is (equal '(2 4) (cf-stream->list (lft::->cf-stream 2.25d0))))
  ;; LFT-STREAM method: sqrt(2)'s first few CF coefficients are the
  ;; well-known periodic [1; 2, 2, 2, ...].
  (is (equal '(1 2 2 2 2) (cf-stream->list (lft::->cf-stream sqrt-two) 5))))

;;; -------------------------------------------------------------------
;;; cf-stream->rational (regression: used to always error, see header)
;;; -------------------------------------------------------------------

(test cf-stream->rational
  (is (= 22/7 (lft::cf-stream->rational (lft::rational->cf-stream 22/7))))
  (is (= -22/7 (lft::cf-stream->rational (lft::rational->cf-stream -22/7))))
  (is (= 3 (lft::cf-stream->rational (lft::rational->cf-stream 3))))
  (is (= 43/19 (lft::cf-stream->rational (lft::coefficients->cf-stream '(2 3 1 4)))))
  (is (eq 'infinity (lft::cf-stream->rational lft::the-empty-stream))))

;;; -------------------------------------------------------------------
;;; cf-stream->lft-stream / lft-stream->cf-stream (round trips)
;;; -------------------------------------------------------------------

(test cf-stream->lft-stream
  (is (= 3 (lft::lft-stream-truncate (lft::cf-stream->lft-stream (lft::rational->cf-stream 22/7)))))
  ;; Regression: negative cf-streams used to error once fully forced
  ;; (CF-STREAM->LFT-STREAM calls NEGATE on them), via NEGATE's bug above.
  (is (< (abs (- -22/7 (lft::lft-stream->double
                        (lft::cf-stream->lft-stream (lft::rational->cf-stream -22/7)))))
         1d-9)))

(test lft-stream->cf-stream
  (is (equal '(3 7) (cf-stream->list (lft::lft-stream->cf-stream (lft::->lft-stream 22/7)))))
  ;; Regression: this used to error (STREAM-CAR on NIL) for any exact
  ;; integer lft-stream, because the fractional part returned by
  ;; LFT-STREAM-TRUNCATE is THE-EMPTY-STREAM, not an LFT-STREAM.
  (is (equal '(5) (cf-stream->list (lft::lft-stream->cf-stream (lft::->lft-stream 5)))))
  (is (equal '(-5) (cf-stream->list (lft::lft-stream->cf-stream (lft::->lft-stream -5)))))
  (is (equal '(-3 -7) (cf-stream->list (lft::lft-stream->cf-stream (lft::->lft-stream -22/7)))))
  (is (equal '(1 2 2 2 2) (cf-stream->list (lft::lft-stream->cf-stream sqrt-two) 5))))

;;; -------------------------------------------------------------------
;;; print-object
;;; -------------------------------------------------------------------

(test print-object-on-cf-stream
  (is (search "3.142857142857142857142" (format nil "~a" (lft::rational->cf-stream 22/7))))
  ;; Regression: used to error (NEGATE's empty-tail bug above).
  (is (search "-3.142857142857142857142" (format nil "~a" (lft::rational->cf-stream -22/7))))
  (is (search "5" (format nil "~a" (lft::rational->cf-stream 5))))
  (is (search "-5" (format nil "~a" (lft::rational->cf-stream -5))))
  (is (search "1.414213562373095048801" (format nil "~a" (lft::->cf-stream sqrt-two)))))

;;; -------------------------------------------------------------------
;;; canonicalize-cf-stream
;;; -------------------------------------------------------------------

(test canonicalize-cf-stream
  ;; A finite coefficient list with no trailing 1 is unchanged.
  (is (equal '(2 4) (cf-stream->list (lft::canonicalize-cf-stream (lft::coefficients->cf-stream '(2 4))))))
  ;; A trailing 1 is folded into the preceding coefficient: [2;3,4,1] ==
  ;; [2;3,5] (both equal 2 + 1/(3 + 1/(4+1)) = 2 + 1/(3 + 1/5)).
  (is (equal '(2 3 5) (cf-stream->list (lft::canonicalize-cf-stream (lft::coefficients->cf-stream '(2 3 4 1))))))
  (is (lft::empty-stream? (lft::canonicalize-cf-stream lft::the-empty-stream))))

;;; -------------------------------------------------------------------
;;; cf-stream-truncate / cf-stream-truncate-delayed
;;; -------------------------------------------------------------------

(test cf-stream-truncate
  (is (equal '(2 3 1) (cf-stream->list (lft::cf-stream-truncate (lft::coefficients->cf-stream '(2 3 1 4 5)) 3))))
  (is (equal '() (cf-stream->list (lft::cf-stream-truncate (lft::coefficients->cf-stream '(2 3 1 4 5)) 0))))
  ;; Regression: requesting more elements than the cf-stream actually has
  ;; used to error (forcing past the end); it must simply stop early.
  (is (equal '(2 3) (cf-stream->list (lft::cf-stream-truncate (lft::coefficients->cf-stream '(2 3)) 10)))))

;;; -------------------------------------------------------------------
;;; reciprocal / negate
;;; -------------------------------------------------------------------

(test reciprocal-cf-stream
  ;; A leading zero coefficient is stripped (1/(0 + x) == x).
  (is (equal '(3 1) (cf-stream->list (lft::reciprocal (lft::coefficients->cf-stream '(0 3 1))))))
  ;; A nonzero leading coefficient gets a new leading zero prepended
  ;; (1/x == (0; x)).
  (is (equal '(0 2 3 1) (cf-stream->list (lft::reciprocal (lft::coefficients->cf-stream '(2 3 1)))))))

(test negate-cf-stream
  ;; Regression: used to error once fully forced (empty-tail bug above).
  (is (equal '(-2 -3 -1) (cf-stream->list (lft::negate (lft::coefficients->cf-stream '(2 3 1))))))
  (is (equal '(-3) (cf-stream->list (lft::negate (lft::coefficients->cf-stream '(3)))))))
;;; NB: NEGATE is not called directly on an already-empty cf-stream
;;; (representing "infinity") by any code path in cf-stream.lisp --
;;; CF-STREAM->LFT-STREAM, its only caller, always guards with
;;; EMPTY-STREAM? first -- so that is intentionally left untested/
;;; unsupported here, consistent with NEGATE's CF-STREAM specializer
;;; never matching THE-EMPTY-STREAM (which is just NIL).

;;; -------------------------------------------------------------------
;;; cf-stream-prefix / cf-stream-suffix
;;; -------------------------------------------------------------------

(test cf-stream-prefix
  (is (equal '(2 3) (cf-stream->list (lft::cf-stream-prefix
                                       (lft::coefficients->cf-stream '(2 3 1))
                                       (lft::coefficients->cf-stream '(2 3 4))))))
  (is (equal '() (cf-stream->list (lft::cf-stream-prefix
                                    (lft::coefficients->cf-stream '(1 3 1))
                                    (lft::coefficients->cf-stream '(2 3 4)))))))

(test cf-stream-suffix
  (is (equal '(4) (cf-stream->list (lft::cf-stream-suffix
                                     (lft::coefficients->cf-stream '(2 3 1))
                                     (lft::coefficients->cf-stream '(2 3 4))))))
  ;; CF-STREAM-SUFFIX returns RIGHT unchanged once LEFT/RIGHT diverge --
  ;; including immediately, on the very first element.
  (is (equal '(2 3 4) (cf-stream->list (lft::cf-stream-suffix
                                         (lft::coefficients->cf-stream '(1 3 1))
                                         (lft::coefficients->cf-stream '(2 3 4)))))))

;;; -------------------------------------------------------------------
;;; cf-stream-flatten-append (regression: used to return its input
;;; unflattened due to swapped arguments, see header)
;;; -------------------------------------------------------------------

(test cf-stream-flatten-append
  (is (equal '(1 2 3 4)
             (cf-stream->list
              (lft::cf-stream-flatten-append
               (lft::cons-stream (lft::coefficients->cf-stream '(1 2))
                                  (lft::cons-stream (lft::coefficients->cf-stream '(3 4))
                                                     lft::the-empty-stream)))))))
