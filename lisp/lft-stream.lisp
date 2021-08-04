;;; -*- Lisp -*-

(in-package "LINEAR-FRACTIONAL-TRANSFORM")

(defclass lft-stream (stream:stream)
  ())

(defmacro cons-lft-stream (lft tail)
  `(make-instance 'lft-stream
                  :car ,lft
                  :delayed-cdr (delay ,tail)))

(defun show-lft-stream (lft-stream)
  (do ((lft-stream lft-stream (stream-cdr lft-stream))
       (count 0 (+ count 1)))
      ((or (empty-stream? lft-stream) (>= count 20)) nil)
    (format t "~s~%" (stream-car lft-stream))))

(defconstant phi
  (if (boundp 'phi)
      phi
      (let ((stream nil))
        (setq stream (cons-lft-stream (make-lft 1 1
                                                1 0)
                                      stream))
        stream)))

(defconstant sqrt-two
  (if (boundp 'sqrt-two)
      sqrt-two
      (let ((stream nil))
        (setq stream (cons-lft-stream (make-lft 1 2
                                                1 1)
                                      stream))
        stream)))

(defmethod funcall-lft (lft (arg lft-stream))
  (cons-lft-stream lft arg))

(defmethod add2 ((left rational) (right lft-stream))
  (funcall (lft-add-rat left) right))

(defmethod add2 ((left lft-stream) (right rational))
  (funcall (lft-add-rat right) left))

(defmethod mul2 ((left rational) (right lft-stream))
  (funcall (lft-multiply-by-rat left) right))

(defmethod mul2 ((left lft-stream) (right rational))
  (funcall (lft-multiply-by-rat right) left))

(defmethod negate ((object lft-stream))
  (funcall lft-negate object))

(defmethod reciprocal ((object lft-stream))
  (funcall lft-reciprocal object))

(defun refine-lft-stream (lft-stream)
  (let ((tail (stream-cdr lft-stream)))
    (if (empty-stream? tail)
        (make-instance 'lft-stream
                       :car (funcall (stream-car lft-stream) (make-lft 1 0 0 0))
                       :delayed-cdr (delay the-empty-stream))
        (make-instance 'lft-stream
                       :car (funcall (stream-car lft-stream) (stream-car tail))
                       :delayed-cdr (stream-delayed-cdr tail)))))

(defun lft-stream-greater-than-rat (lft-stream rat receive-true receive-false)
  (lft-greater-than-rat
   (stream-car lft-stream) rat
   (lambda () (funcall receive-true lft-stream))
   (lambda () (funcall receive-false lft-stream))
   (lambda () (lft-stream-greater-than-rat (refine-lft-stream lft-stream) rat receive-true receive-false))))

(defun lft-stream-less-than-rat (lft-stream rat receive-true receive-false)
  (lft-less-than-rat
   (stream-car lft-stream) rat
   (lambda () (funcall receive-true lft-stream))
   (lambda () (funcall receive-false lft-stream))
   (lambda () (lft-stream-less-than-rat (refine-lft-stream lft-stream) rat receive-true receive-false))))

(defun lft-stream-minus-p (lft-stream receive-true receive-false)
  (lft-minus-p
   (stream-car lft-stream)
   (lambda () (funcall receive-true lft-stream))
   (lambda () (funcall receive-false lft-stream))
   (lambda () (lft-stream-minus-p (refine-lft-stream lft-stream) receive-true receive-false))))

(defun lft-stream-truncate (lft-stream)
  (lft-truncate
   (stream-car lft-stream)
   (lambda (integer-part fractional-part)
     (values integer-part
             (lft-zero-p fractional-part
                         (lambda () the-empty-stream)
                         (lambda () (make-instance 'lft-stream
                                                   :car fractional-part
                                                   :delayed-cdr (stream-delayed-cdr lft-stream)))
                         (lambda () (make-instance 'lft-stream
                                                   :car fractional-part
                                                   :delayed-cdr (stream-delayed-cdr lft-stream))))))
   (lambda () (lft-stream-truncate (refine-lft-stream lft-stream)))))

(defun lft-stream-zero-p (lft-stream)
  (lft-zero-p (stream-car lft-stream)
              (lambda () (values t lft-stream))
              (lambda () (values nil lft-stream))
              (lambda () (lft-stream-zero-p (refine-lft-stream lft-stream)))))

(defun lft-stream->radix-stream (lft-stream radix)
  (if (empty-stream? lft-stream)
      the-empty-stream
      (multiple-value-bind (integer-part fractional-part) (lft-stream-truncate lft-stream)
        (cons-stream integer-part
                     (if (empty-stream? fractional-part)
                         the-empty-stream
                         (lft-stream->radix-stream (x* fractional-part radix) radix))))))

(defvar *print-lft-stream-length* 20)

(defun print-lft-stream-as-decimal (lft-stream cl:stream)
  (labels ((print-digits (digit-stream count)
             (unless (empty-stream? digit-stream)
               (format cl:stream "~d~@[.~]" (abs (stream-car digit-stream)) (zerop count))
               (if (> count *print-lft-stream-length*)
                   (format cl:stream "...")
                   (print-digits (stream-cdr digit-stream) (+ count 1))))))
    (print-digits
     (lft-stream->radix-stream
      (lft-stream-minus-p
       lft-stream
       (lambda (refined)
         (format cl:stream "-")
         (negate refined))
       #'identity)
      10)
     0)))

(defvar *print-lft-stream-as-decimal* t)

(defmethod print-object ((object lft-stream) cl:stream)
  (print-unreadable-object (object cl:stream :type t)
    (when *print-lft-stream-as-decimal*
      (print-lft-stream-as-decimal object cl:stream))))

(defun lft-compose-refine? (left right if-success if-failure)
  (let ((attempt (funcall left right)))
    (if (refine? attempt)
        (funcall if-success attempt)
        (funcall if-failure))))

(defconstant lft-digit-minus-one
  (if (boundp 'lft-digit-minus-one)
      lft-digit-minus-one
      (make-lft 1 0
                1 2)))

(defun lft-digit-minus-one (object)
  (funcall lft-digit-minus-one object))

(defconstant lft-digit-one
  (if (boundp 'lft-digit-one)
      lft-digit-one
      (make-lft 2 1
                0 1)))

(defun lft-digit-one (object)
  (funcall lft-digit-one object))

(defconstant lft-digit-zero
  (if (boundp 'lft-digit-zero)
      lft-digit-zero
      (make-lft 3 1
                1 3)))

(defun lft-digit-zero (object)
  (funcall lft-digit-zero object))

(defconstant inverse-lft-digit-minus-one
  (if (boundp 'inverse-lft-digit-minus-one)
      inverse-lft-digit-minus-one
      (inverse-lft lft-digit-minus-one)))

(defun inverse-lft-digit-minus-one (object)
  (funcall inverse-lft-digit-minus-one object))

(defconstant inverse-lft-digit-one
  (if (boundp 'inverse-lft-digit-one)
      inverse-lft-digit-one
      (inverse-lft lft-digit-one)))

(defun inverse-lft-digit-one (object)
  (funcall inverse-lft-digit-one object))

(defconstant inverse-lft-digit-zero
  (if (boundp 'inverse-lft-digit-zero)
      inverse-lft-digit-zero
      (inverse-lft lft-digit-zero)))

(defun inverse-lft-digit-zero (object)
  (funcall inverse-lft-digit-zero object))

(defun lft-try-factor-digit (lft if-success if-failure)
  (lft-compose-refine?
   inverse-lft-digit-one lft
   (lambda (new-lft)
     (funcall if-success lft-digit-one new-lft))
   (lambda ()
     (lft-compose-refine?
      inverse-lft-digit-minus-one lft
      (lambda (new-lft)
        (funcall if-success lft-digit-minus-one new-lft))
      (lambda ()
        (lft-compose-refine?
         inverse-lft-digit-zero lft
         (lambda (new-lft)
           (funcall if-success lft-digit-zero new-lft))
         if-failure))))))

(defun lft-stream-factor-digit (lft-stream)
  (lft-try-factor-digit
   (stream-car lft-stream)
   (lambda (digit factor)
     (values digit (make-instance 'lft-stream
                                  :car factor
                                  :delayed-cdr (stream-delayed-cdr lft-stream))))
   (lambda ()
     (lft-stream-factor-digit (refine-lft-stream lft-stream)))))

(defun lft-stream->lft-digit-stream (lft-stream)
  (multiple-value-bind (digit remainder) (lft-stream-factor-digit lft-stream)
    (cons-lft-stream digit (lft-stream->lft-digit-stream remainder))))

(defconstant lft-sign-infinity
  (if (boundp 'lft-sign-infinity)
      lft-sign-infinity
      (make-lft 1  1
                -1 1)))

(defun lft-sign-infinity (object)
  (funcall lft-sign-infinity object))

(defconstant lft-sign-negative
  (if (boundp 'lft-sign-negative)
      lft-sign-negative
      (make-lft 0 -1
                1  0)))

(defun lft-sign-negative (object)
  (funcall lft-sign-negative object))

(defconstant lft-sign-positive
  (if (boundp 'lft-sign-positive)
      lft-sign-positive
      (make-lft 1 0
                0 1)))

(defun lft-sign-positive (object)
  (funcall lft-sign-positive object))

(defconstant lft-sign-zero
  (if (boundp 'lft-sign-zero)
      lft-sign-zero
      (make-lft 1 -1
                1  1)))

(defun lft-sign-zero (object)
  (funcall lft-sign-zero object))

(defconstant inverse-lft-sign-infinity
  (if (boundp 'inverse-lft-sign-infinity)
      inverse-lft-sign-infinity
      (inverse-lft lft-sign-infinity)))

(defun inverse-lft-sign-infinity (object)
  (funcall inverse-lft-sign-infinity object))

(defconstant inverse-lft-sign-negative
  (if (boundp 'inverse-lft-sign-negative)
      inverse-lft-sign-negative
      (inverse-lft lft-sign-negative)))

(defun inverse-lft-sign-negative (object)
  (funcall inverse-lft-sign-negative object))

(defconstant inverse-lft-sign-positive
  (if (boundp 'inverse-lft-sign-positive)
      inverse-lft-sign-positive
      (inverse-lft lft-sign-positive)))

(defun inverse-lft-sign-positive (object)
  (funcall inverse-lft-sign-positive object))

(defconstant inverse-lft-sign-zero
  (if (boundp 'inverse-lft-sign-zero)
      inverse-lft-sign-zero
      (inverse-lft lft-sign-zero)))

(defun inverse-lft-sign-zero (object)
  (funcall inverse-lft-sign-zero object))

(defun lft-try-factor-sign (lft if-success if-failure)
  (lft-compose-refine?
   inverse-lft-sign-positive lft
   (lambda (new-lft)
     (funcall if-success lft-sign-positive new-lft))
   (lambda ()
     (lft-compose-refine?
      inverse-lft-sign-negative lft
      (lambda (new-lft)
        (funcall if-success lft-sign-negative new-lft))
      (lambda ()
        (lft-compose-refine?
         inverse-lft-sign-zero lft
         (lambda (new-lft)
           (funcall if-success lft-sign-zero new-lft))
         (lambda ()
           (lft-compose-refine?
            inverse-lft-sign-infinity lft
            (lambda (new-lft)
              (funcall if-success lft-sign-infinity new-lft))
            if-failure))))))))

(defun lft-stream-factor-sign (lft-stream)
  (lft-try-factor-sign
   (stream-car lft-stream)
   (lambda (sign factor)
     (values sign (make-instance 'lft-stream
                                 :car factor
                                 :delayed-cdr (stream-delayed-cdr lft-stream))))
   (lambda ()
     (lft-stream-factor-sign (refine-lft-stream lft-stream)))))

(defun canonicalize-lft-stream (lft-stream)
  (multiple-value-bind (sign remainder) (lft-stream-factor-sign lft-stream)
    (cons-lft-stream sign (lft-stream->lft-digit-stream remainder))))

(defun lft-stream->double (lft-stream)
  (lft->double
   (stream-car lft-stream)
   #'identity
   (lambda () (lft-stream->double (refine-lft-stream lft-stream)))))

(defun %rat->lft-stream (numerator denominator)
  (multiple-value-bind (integer-part remainder) (truncate numerator denominator)
    (cons-lft-stream (make-lft integer-part 1
                               1 0)
                     (if (zerop remainder)
                         the-empty-stream
                         (%rat->lft-stream denominator remainder)))))

(defun rat->lft-stream (rat)
  (etypecase rat
    (integer (%rat->lft-stream rat 1))
    (rational (%rat->lft-stream (numerator rat) (denominator rat)))
    ((or single-float double-float) (%rat->lft-stream (numerator (rational rat)) (denominator (rational rat))))))

(defun %rat-atan (x)
  (labels ((next-term (addends numerators)
             (cons-lft-stream (make-lft (stream-car addends) (stream-car numerators)
                                        1                    0)
                              (next-term (stream-cdr addends) (stream-cdr numerators)))))
    (cons-lft-stream
     (make-lft 0 x
               1 0)
     (next-term (odds) (stream-map (lambda (square) (* square (square x))) (squares))))))

(defun %exp-rat (x)
  ;; 1/2 < x <= 2
  (labels ((l (n)
             (cons-lft-stream
              (make-lft (+ (* 4 n) 2) x
                        x             0)
              (l (+ n 1)))))
    (cons-lft-stream
     (make-lft (+ 2 x) x
               (- 2 x) x)
     (l 1))))

(defun %log-rat (x)
  (if (< x 1)
      (error "%log-rat called on small x"))
  (labels ((l (n)
             (cons-lft-stream
              (funcall (make-lft 0             (- x 1)
                                 (- x 1) (+ (* n 2) 1))
                       (make-lft 0       (+ n 1)
                                 (+ n 1)       2))
              (l (+ n 1)))))
    (l 0)))

(defun %rat-pow (x y)
  (labels ((l (n)
             (cons-lft-stream
              (funcall (make-lft 0             (- x 1)
                                 (- x 1) (- (* n 2) 1))
                       (make-lft 0       (- n y)
                                 (+ n y)       2))
              (l (+ n 1)))))
    (cons-lft-stream
     (make-lft y 1
               0 1)
     (l 1))))

(defun %sqrt-rat (p q)
  (labels ((rollover (a b c)
             (let ((d (+ (* 2 (- b a)) c)))
               (if (> d 0)
                   (cons-lft-stream lft-digit-minus-one (rollover (* a 4) d c))
                   (cons-lft-stream lft-digit-one (rollover (- d) (* b 4) c))))))
    (rollover p q (- p q))))

(defmethod x-sqrt ((number rational))
  (%sqrt-rat (numerator number) (denominator number)))

(defconstant ramanujan-pi-stream
  ;; sqrt(10005)/pi
  (if (boundp 'ramanujan-pi-stream)
      ramanujan-pi-stream
      (labels ((l (n)
                 (let ((c (* (- (* 2 n) 1)
                             (- (* 6 n) 5)
                             (- (* 6 n) 1)
                             (+ (* 545140134 n)
                                13591409)))
                       (d (* (- (* 2 n) 1)
                             (- (* 6 n) 5)
                             (- (* 6 n) 1)
                             (+ n 1)))
                       (e (* 10939058860032000 n n n n)))
                   (cons-lft-stream (make-lft (- e d c) (- (+ e d) c)
                                              (+ e d c) (- e (- d c)))
                                    (l (+ n 1))))))
        (cons-lft-stream (make-lft 6795705 6795704
                                   213440  213440)
                         (l 1)))))
