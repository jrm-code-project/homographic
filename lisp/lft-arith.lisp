;;; -*- Lisp -*-

(in-package "LINEAR-FRACTIONAL-TRANSFORM")

(defconstant pi
  (if (boundp 'pi)
      pi
      (x/ (x-sqrt 10005) ramanujan-pi-stream)))

(defconstant 2pi
  (if (boundp '2pi)
      2pi
      (x* (x/ (x-sqrt 10005) ramanujan-pi-stream) 2)))

(defconstant pi/2
  (if (boundp 'pi/2)
      pi/2
      (x/ (x/ (x-sqrt 10005) ramanujan-pi-stream) 2)))

(defun atan-rat (rat)
  (cond ((minusp rat) (funcall lft-negate (atan-rat (- rat))))
        ((> rat 1) (make-instance 'binary-expression
                                  :bilft (make-bilft 0 1 -1 0
                                                     0 0  0 1)
                                  :delayed-left (delay (x/ pi 2))
                                  :delayed-right (delay (atan-rat (/ 1 rat)))))
        (t (%rat-atan rat))))

(defun exp-rat (rat)
  (cond ((minusp rat) (reciprocal (exp-rat (negate rat))))
        ((> rat 1) (%square-lft-stream (exp-rat (/ rat 2))))
        (t (%exp-rat rat))))

(defmethod x-exp ((number rational))
  (exp-rat number))

(defun exp-lft-stream (lft-stream)
  (lft-stream-minus-p
   lft-stream
   (lambda (refined)
     (funcall lft-reciprocal (exp-lft-stream (funcall lft-negate refined))))
   (lambda (refined)
     (lft-stream-greater-than-rat
      refined 1
      (lambda (refined*)
        (%square-lft-stream
         (exp-lft-stream
          (x/ refined* 2))))
      #'%exp-lft-stream))))

(defmethod x-exp ((number lft-stream))
  (exp-lft-stream number))

(defconstant log2
  (if (boundp 'log2)
      log2
      (%log-rat 2)))

(defun log-rat (rat)
  (cond ((minusp rat) (error "Cannot take the log of a negative number."))
        ((< rat 1) (negate (log-rat (/ 1 rat))))
        ((> rat 2) (funcall bilft-add
                            log2
                            (delay (log-rat (/ rat 2)))))
        (t (%log-rat rat))))

(defmethod x-log ((number rational))
  (log-rat number))

(defun log-lft-stream (lft-stream)
  (lft-stream-minus-p
   lft-stream
   (lambda (refined)
     (declare (ignore refined))
     (error "Cannot take the log of a negative number."))
   (lambda (refined)
     (lft-stream-less-than-rat
      refined 1
      (lambda (refined*)
        (negate (log-lft-stream (reciprocal refined*))))
      (lambda (refined*)
        (lft-stream-greater-than-rat
         refined* 2
         (lambda (refined**)
           (funcall bilft-add
                    log2
                    (delay (log-lft-stream (x/ refined** 2)))))
         #'%log-lft-stream))))))

(defmethod x-log ((number lft-stream))
  (log-lft-stream number))

(defmethod x-cbrt ((number lft-stream))
  (x-exp (x/ (x-log number) 3)))

(defmethod x-sqrt ((number lft-stream))
  (lft-stream-less-than-rat
   number 1
   #'%sqrt-lft-stream
   (lambda (number*)
     (x-exp (x/ (x-log number*) 2)))))

(defun tan-lft-stream (lft-stream)
  (lft-stream-minus-p
   lft-stream
   (lambda (refined)
     (cons-lft-stream lft-negate (tan-lft-stream (negate refined))))
   (lambda (refined)
     (lft-stream-less-than-rat
      refined 1
      #'%tan-lft-stream
      (lambda (refined*)
        (let ((tan-half-x (delay (tan-lft-stream (x/ refined* 2)))))
          (binary-expression->lft-stream
           (make-instance 'binary-expression
                          :generation 0
                          :bilft (make-bilft 0  1 1 0
                                             -1 0 0 1)
                          :delayed-left tan-half-x
                          :delayed-right tan-half-x))))))))

(defun cos-lft-stream (lft-stream)
  (let ((tan-half-x (delay (tan-lft-stream (x/ lft-stream 2)))))
    (binary-expression->lft-stream
     (make-instance 'binary-expression
                    :generation 0
                    :bilft (make-bilft -1 0 0 1
                                       1 0 0 1)
                    :delayed-left tan-half-x
                    :delayed-right tan-half-x))))

(defun sin-lft-stream (lft-stream)
  (let ((tan-half-x (delay (tan-lft-stream (x/ lft-stream 2)))))
    (binary-expression->lft-stream
     (make-instance 'binary-expression
                    :generation 0
                    :bilft (make-bilft 0 1 1 0
                                       1 0 0 1)
                    :delayed-left tan-half-x
                    :delayed-right tan-half-x))))

(defun log-base (number base)
  (x/ (etypecase number
        ((or integer rational) (log-rat number))
        (lft-stream (log-lft-stream number)))
      (etypecase base
        ((or integer rational) (log-rat base))
        (lft-stream (log-lft-stream base)))))

(defun expt-expression (base exponent)
  (etypecase exponent
    ((or integer rational)
     (etypecase base
       ((or integer rational) (%rat-pow base exponent))
       (lft-stream (exp-lft-stream
                    (funcall (lft-multiply-by-rat exponent)
                             (log-lft-stream base))))))
    (lft-stream
     (exp-lft-stream
      (funcall bilft-multiply
               exponent
               (etypecase base
                 ((or integer rational) (log-rat base))
                 (lft-stream (log-lft-stream base))))))))

(defun hyperbolic-function (bilft x)
  (let ((exp (x-exp x)))
    (binary-expression->lft-stream
     (make-instance 'binary-expression
                    :generation 0
                    :bilft (compose-bilft-lft-x
                            (compose-bilft-lft-y
                             bilft
                             (stream-car exp))
                            (stream-car exp))
                    :delayed-left (stream-delayed-cdr exp)
                    :delayed-right (stream-delayed-cdr exp)))))

(defun x-sinh (x)
  (hyperbolic-function (make-bilft 1 0 0 -1
                                   0 1 1  0)
                       x))

(defun x-cosh (x)
  (hyperbolic-function (make-bilft 1 0 0 1
                                   0 1 1 0)
                       x))

(defun x-tanh (x)
  (hyperbolic-function (make-bilft 1 0 0 -1
                                   1 0 0 1)
                       x))

(defun x-coth (x)
  (reciprocal (x-tanh x)))

(defun x-sech (x)
  (reciprocal (x-cosh x)))

(defun x-csch (x)
  (reciprocal (x-sinh x)))

(defun asinh-lft-stream (lft-stream)
  (let ((sqrt (x-sqrt
               (binary-expression->lft-stream
                (make-instance 'binary-expression
                               :bilft (compose-bilft-lft-y
                                       (compose-bilft-lft-x
                                        (make-bilft 1 0 0 1
                                                    0 0 0 1)
                                        (stream-car lft-stream))
                                       (stream-car lft-stream))
                               :delayed-left (stream-delayed-cdr lft-stream)
                               :delayed-right (stream-delayed-cdr lft-stream))))))
    (log-lft-stream
     (binary-expression->lft-stream
      (make-instance 'binary-expression
                     :bilft (compose-bilft-lft-x
                             (compose-bilft-lft-y
                              (make-bilft 0 1 1 0
                                          0 0 0 1)
                              (stream-car sqrt))
                             (stream-car lft-stream))
                     :delayed-left (stream-delayed-cdr lft-stream)
                     :delayed-right (stream-delayed-cdr sqrt))))))

(defun acosh-lft-stream (lft-stream)
  (let ((sqrt (x-sqrt
               (binary-expression->lft-stream
                (make-instance 'binary-expression
                               :bilft (compose-bilft-lft-y
                                       (compose-bilft-lft-x
                                        (make-bilft 1 0 0 -1
                                                    0 0 0 1)
                                        (stream-car lft-stream))
                                       (stream-car lft-stream))
                               :delayed-left (stream-delayed-cdr lft-stream)
                               :delayed-right (stream-delayed-cdr lft-stream))))))
    (log-lft-stream
     (binary-expression->lft-stream
      (make-instance 'binary-expression
                     :bilft (compose-bilft-lft-x
                             (compose-bilft-lft-y
                              (make-bilft 0 1 1 0
                                          0 0 0 1)
                              (stream-car sqrt))
                             (stream-car lft-stream))
                     :delayed-left (stream-delayed-cdr lft-stream)
                     :delayed-right (stream-delayed-cdr sqrt))))))

(defun atanh-lft-stream (lft-stream)
  (funcall (lft-multiply-by-rat 1/2)
           (log-lft-stream
            (funcall lft-sign-infinity lft-stream))))

