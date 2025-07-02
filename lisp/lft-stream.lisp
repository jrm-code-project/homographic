;;; -*- Lisp -*-

;;; The algorithms implemented here were developed by Peter Potts in
;;; his doctoral thesis.

(in-package "LINEAR-FRACTIONAL-TRANSFORMATION")

(eval-when (:compile-toplevel :load-toplevel :execute)
(defclass lft-stream (stream::stream)
  ())

(defmacro cons-lft-stream (lft tail)
  `(make-instance 'lft-stream
                  :car ,lft
                  :delayed-cdr (delay ,tail)))

(defmacro delay-lft-stream (lft-stream)
  `(cons-lft-stream lft-identity ,lft-stream))

)

(defun show-lft-stream (lft-stream &optional (length 5))
  (do ((lft-stream lft-stream (stream-cdr lft-stream))
       (count 0 (+ count 1)))
      ((or (empty-stream? lft-stream) (>= count length))
       (unless (empty-stream? lft-stream)
         (format t "   ..."))
       (format t "}"))
    (format t "~:[   ~;{~]~s~%" (zerop count) (stream-car lft-stream))))

(defun lft-stream-map (f stream) ;; ugh.
  (if (empty-stream? stream)
      the-empty-stream
      (cons-lft-stream (funcall f (stream-car stream))
                       (lft-stream-map f (stream-cdr stream)))))

(defun circular-lft-stream (lft)
  (let ((stream nil))
    (setq stream (cons-lft-stream lft stream))
    stream))

(defparameter phi (circular-lft-stream (make-lft 1 1
                                                 1 0)))

(defparameter sqrt-two (circular-lft-stream (make-lft 1 2
                                                      1 1)))

(defmethod funcall-lft (lft (arg lft-stream))
  ;; (cons-lft-stream lft arg)
  (make-instance 'lft-stream
                 :car (funcall lft (stream-car arg))
                 :delayed-cdr (stream-delayed-cdr arg)))

(defmethod add2 ((left cl:rational) (right lft-stream))
  (funcall (lft-add-rat left) right))

(defmethod add2 ((left lft-stream) (right cl:rational))
  (funcall (lft-add-rat right) left))

(defmethod multiply2 ((left cl:rational) (right lft-stream))
  (funcall (lft-multiply-by-rat left) right))

(defmethod multiply2 ((left lft-stream) (right cl:rational))
  (funcall (lft-multiply-by-rat right) left))

(defmethod negate ((object lft-stream))
  (funcall lft-negate object))

(defmethod reciprocal ((object lft-stream))
  (funcall lft-reciprocal object))

(defun refine-lft-stream (lft-stream)
  (let ((tail (stream-cdr lft-stream)))
    (if (empty-stream? tail)
        (make-instance 'lft-stream
                       :car (funcall (stream-car lft-stream) (make-lft 0 1 0 0))
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

(defun lft-stream-plus-p (lft-stream receive-true receive-false)
  (lft-plus-p
   (stream-car lft-stream)
   (lambda () (funcall receive-true lft-stream))
   (lambda () (funcall receive-false lft-stream))
   (lambda () (lft-stream-plus-p (refine-lft-stream lft-stream) receive-true receive-false))))

(defun lft-stream-non-negative-p (lft-stream receive-true receive-false)
  (lft-non-negative-p
   (stream-car lft-stream)
   (lambda () (funcall receive-true lft-stream))
   (lambda () (funcall receive-false lft-stream))
   (lambda () (lft-stream-non-negative-p (refine-lft-stream lft-stream) receive-true receive-false))))

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
                         (lft-stream->radix-stream (* fractional-part radix) radix))))))

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

(defun decompose (left composed)
  (values left (funcall (inverse-lft left) composed)))

(defun decompose-range? (left composed if-success if-failure)
  (multiple-value-bind (left right) (decompose left composed)
    (if (range? right)
        (funcall if-success left right)
        (funcall if-failure))))

(defparameter lft-digit-low
  (if (boundp 'lft-digit-low)
      lft-digit-low
      (make-lft 1 0
                1 2)))

(defun lft-digit-low (object)
  (funcall lft-digit-low object))

(defparameter lft-digit-high
  (if (boundp 'lft-digit-high)
      lft-digit-high
      (make-lft 2 1
                0 1)))

(defun lft-digit-high (object)
  (funcall lft-digit-high object))

(defparameter lft-digit-zero
  (if (boundp 'lft-digit-zero)
      lft-digit-zero
      (make-lft 3 1
                1 3)))

(defun lft-digit-zero (object)
  (funcall lft-digit-zero object))

(defun try-decompose-digit (lft if-success if-failure)
  (decompose-range?
   lft-digit-high lft
   if-success
   (lambda ()
     (decompose-range?
      lft-digit-low lft
      if-success
      (lambda ()
        (decompose-range?
         lft-digit-zero lft
         if-success
         if-failure))))))

(defun lft-stream-decompose-digit (lft-stream)
  (try-decompose-digit
   (stream-car lft-stream)
   (lambda (digit remainder)
     (values digit (make-instance 'lft-stream
                                  :car remainder
                                  :delayed-cdr (stream-delayed-cdr lft-stream))))
   (lambda ()
     (lft-stream-decompose-digit (refine-lft-stream lft-stream)))))

(defun lft-stream->lft-digit-stream (lft-stream)
  (multiple-value-bind (digit remainder) (lft-stream-decompose-digit lft-stream)
    (cons-lft-stream digit (lft-stream->lft-digit-stream remainder))))

(defparameter lft-sign-infinity
  (if (boundp 'lft-sign-infinity)
      lft-sign-infinity
      (make-lft 1  1
                -1 1)))

(defun lft-sign-infinity (object)
  (funcall lft-sign-infinity object))

(defparameter lft-sign-negative
  (if (boundp 'lft-sign-negative)
      lft-sign-negative
      (make-lft 0 -1
                1  0)))

(defun lft-sign-negative (object)
  (funcall lft-sign-negative object))

(defparameter lft-sign-positive
  (if (boundp 'lft-sign-positive)
      lft-sign-positive
      (make-lft 1 0
                0 1)))

(defun lft-sign-positive (object)
  (funcall lft-sign-positive object))

(defparameter lft-sign-zero
  (if (boundp 'lft-sign-zero)
      lft-sign-zero
      (make-lft 1 -1
                1  1)))

(defun lft-sign-zero (object)
  (funcall lft-sign-zero object))

(defun try-decompose-sign (lft if-success if-failure)
  (decompose-range?
   lft-sign-positive lft
   if-success
   (lambda ()
     (decompose-range?
      lft-sign-negative lft
      if-success
      (lambda ()
        (decompose-range?
         lft-sign-zero lft
         if-success
         (lambda ()
           (decompose-range?
            lft-sign-infinity lft
            if-success
            if-failure))))))))

(defun lft-stream-decompose-sign (lft-stream)
  (try-decompose-sign
   (stream-car lft-stream)
   (lambda (sign remainder)
     (values sign (make-instance 'lft-stream
                                 :car remainder
                                 :delayed-cdr (stream-delayed-cdr lft-stream))))
   (lambda ()
     (lft-stream-decompose-sign (refine-lft-stream lft-stream)))))

(defun canonicalize-lft-stream (lft-stream)
  (multiple-value-bind (sign remainder) (lft-stream-decompose-sign lft-stream)
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

(defgeneric ->lft-stream (number)
  (:method ((number integer))
    (%rat->lft-stream number 1))
  (:method ((number cl:rational))
    (%rat->lft-stream (numerator number) (denominator number)))
  (:method ((number single-float))
    (->lft-stream (rational number)))
  (:method ((number double-float))
    (->lft-stream (rational number))))

;;; Reinhold Heckmann
(defun %sqrt-rat (p q)
  (let ((diff (- p q)))
    (labels ((rollover (num den)
               (let ((d (+ (* 2 (- den num)) diff)))
                 (if (> d 0)
                     (cons-lft-stream lft-digit-low (rollover (* num 4) d))
                     (cons-lft-stream lft-digit-high (rollover (- d) (* den 4)))))))
      (cons-lft-stream lft-sign-positive (rollover p q)))))

(defmethod sqrt ((number cl:rational))
  (%sqrt-rat (numerator number) (denominator number)))

(defun %exp-rat (x)
  (check-type x (cl:rational (1/2) 2))
  (cons-lft-stream
   lft-sign-positive
   (cons-lft-stream
    (make-lft (+ 2 x) x
              (- 2 x) x)
    (lft-stream-map (lambda (n)
                      (make-lft (+ (* 4 n) 2) x
                                x             0))
                    (naturals)))))

(defun %log-rat (x)
  (check-type x (cl:rational (1) cl:*))
  (cons-lft-stream
   lft-sign-positive
   (lft-stream-map
    (lambda (n)
      (funcall (make-lft 0             (- x 1)
                         (- x 1) (+ (* n 2) 1))
               (make-lft 0       (+ n 1)
                         (+ n 1)       2)))
    (integers))))

(defun %rat-pow (x y)
  (check-type x (cl:rational (1) cl:*))
  (check-type y (cl:rational (0) (1)))
  (cons-lft-stream
   lft-sign-positive
   (cons-lft-stream
    (make-lft y 1
              0 1)
    (lft-stream-map
     (lambda (n)
       (funcall (make-lft 0             (- x 1)
                          (- x 1) (- (* n 2) 1))
                (make-lft 0       (- n y)
                          (+ n y)       2)))
     (naturals)))))

(defparameter ramanujan-pi-stream
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

(defun big-k-stream (numerators denominators)
  (cons-lft-stream (make-lft 0 (stream-car numerators)
                             1 (stream-car denominators))
                   (big-k-stream (stream-cdr numerators) (stream-cdr denominators))))

(defun %rat-atan (z)
  (check-type z (cl:rational (0) 1))
  (let ((z-squared (square z)))
    (cons-lft-stream (make-lft 0 z
                               1 1)
                     (big-k-stream (stream-map (lambda (square)
                                                 (* z-squared square))
                                               (squares))
                                   (stream-cdr (odds))))))

(defun %rat-tan (x)
  (check-type x (cl:rational (0) 1))
  (lft-stream-map
   (lambda (n)
     (funcall (make-lft 0 x
                        x (+ (* n 4) 1))
              (make-lft 0 x
                        x (- (* n -4) 3))))
   (integers)))

(defparameter euler-gompertz
  (big-k-stream (cons-stream 1 (double-stream (naturals)))
                (ones)))



