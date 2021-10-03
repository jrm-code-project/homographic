;;; -*- Lisp -*-

(in-package "LINEAR-FRACTIONAL-TRANSFORM")

(defparameter pi (x/ (x-sqrt 10005) ramanujan-pi-stream))
(defparameter 2pi (x* pi 2))
(defparameter pi/2 (x/ pi 2))
(defparameter pi/3 (x/ pi 3))
(defparameter pi/4 (x/ pi 4))
(defparameter pi/6 (x/ pi 6))

(defun atan-rat (rat)
  (cond ((minusp rat) (funcall lft-negate (atan-rat (- rat))))
        ((> rat 1) (make-instance 'binary-expression
                                  :bilft (make-bilft 0 1 -1 0
                                                     0 0  0 1)
                                  :delayed-left (delay pi/2)
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
     (reciprocal (exp-lft-stream (negate refined))))
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

(defparameter log2 (%log-rat 2))

(defun log-rat (rat)
  (labels ((divide-down (rat powers-of-two)
             (if (> rat 2)
                 (divide-down (/ rat 2) (+ powers-of-two 1))
                 (funcall bilft-add
                          (x* powers-of-two log2)
                          (delay (%log-rat rat))))))

    (cond ((minusp rat) (error "Cannot take the log of a negative number."))
          ((< rat 1) (negate (log-rat (/ 1 rat))))
          (t (divide-down rat 0)))))

(defmethod x-log ((number rational))
  (log-rat number))

(defun log-lft-stream (lft-stream)
  (labels ((divide-down (lft-stream powers-of-two)
             (lft-stream-greater-than-rat
              lft-stream 2
              (lambda (refined*)
                (divide-down (x/ refined* 2) (+ powers-of-two 1)))
              (lambda (refined*)
                (funcall bilft-add
                         (x* powers-of-two log2)
                         (delay (%log-lft-stream refined*)))))))

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
             (divide-down refined** 0))
           #'%log-lft-stream)))))))

(defmethod x-log ((number lft-stream))
  (log-lft-stream number))

(defmethod x-cbrt ((number lft-stream))
  (x-exp (x/ (x-log number) 3)))

(defmethod x-sqrt ((number lft-stream))
  (%sqrt-lft-stream number))

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

;;;

(defparameter pythagoras (x-sqrt 2))
(defparameter theodorus (x-sqrt 3))
(defparameter sqrt-five (x-sqrt 5))

(defun approximate-zeta-three (nterms)
  (flet ((f (k)
            (* (expt -1 k)
               (/ (* (expt (factorial (+ (* 2 k) 1)) 3)
                     (expt (factorial (* 2 k)) 3)
                     (expt (factorial k) 3)
                     (+ (* 126392 k k k k k)
                        (* 412708 k k k k)
                        (* 531578 k k k)
                        (* 336367 k k)
                        (* 104000 k)
                        12463))
                  (* (factorial (+ (* 3 k) 2))
                     (expt (factorial (+ (* 4 k) 3)) 3))))))
    (/ (big-sigma #'f 0 nterms) 24)))

(defparameter apéry (limit-stream->cf-stream (stream-map #'approximate-zeta-three (integers))))
