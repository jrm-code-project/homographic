;;; -*- Lisp -*-

(in-package "LINEAR-FRACTIONAL-TRANSFORMATION")

(defparameter pi (/ (sqrt 10005) ramanujan-pi-stream))
(defparameter 2pi (* pi 2))
(defparameter pi/2 (/ pi 2))
(defparameter pi/3 (/ pi 3))
(defparameter pi/4 (/ pi 4))
(defparameter pi/6 (/ pi 6))
(defparameter pi/8 (/ pi 8))

(defmethod square ((number lft-stream))
  (funcall (make-bilft 1 0 0 0
                       0 0 0 1)
           number
           number))

(defun exp-rat (rat)
  (cond ((minusp rat) (reciprocal (exp-rat (negate rat))))
        ((<= rat 1/2) (reciprocal (%sqrt-lft-stream (reciprocal (exp-rat (* rat 2))))))
        ((> rat 2) (%square-lft-stream (exp-rat (/ rat 2))))
        (t (%exp-rat rat))))

(defun exp-lft-stream (lft-stream)
  (lft-stream-minus-p
   lft-stream
   (lambda (refined)
     (reciprocal (exp-lft-stream (negate refined))))
   (lambda (refined)
     (lft-stream-greater-than-rat
      refined 2
      (lambda (refined*)
        (%square-lft-stream
         (exp-lft-stream
          (/ refined* 2))))
      #'%exp-lft-stream))))

(defmethod exp ((number cl:rational))
  (exp-rat number))

(defmethod exp ((number lft-stream))
  (exp-lft-stream number))

(defparameter log2 (%log-rat 2))

(defun log-rat (rat)
  (labels ((divide-down (rat powers-of-two)
             (if (> rat 2)
                 (divide-down (/ rat 2) (+ powers-of-two 1))
                 (+
                  (* powers-of-two log2)
                  (%log-rat rat)))))

    (cond ((minusp rat) (error "Cannot take the log of a negative number."))
          ((< rat 1) (negate (log-rat (/ 1 rat))))
          ((= rat 1) 0)
          ((< rat 2) (%log-rat rat))
          (t (divide-down (/ rat 2) 1)))))

(defun log-lft-stream (lft-stream)
  (labels ((divide-down (lft-stream powers-of-two)
             (lft-stream-greater-than-rat
              lft-stream 2
              (lambda (refined*)
                (divide-down (/ refined* 2) (+ powers-of-two 1)))
              (lambda (refined*)
                (funcall bilft-add
                         (* powers-of-two log2)
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

(defmethod log ((number cl:rational) &optional base)
  (log-rat number))

(defmethod log ((number lft-stream) &optional base)
  (log-lft-stream number))

(defun rat-pow (base exponent)
  (check-type base (cl:rational 0 cl:*))
  (check-type exponent (cl:rational cl:* cl:*))
  (cond ((zerop base) (cond ((minusp exponent) (error 'division-by-zero))
                            ((zerop exponent) 1) ;; ?
                            (t 0)))
        ((< base 1) (reciprocal (rat-pow (reciprocal base) exponent)))
        ((= base 1) 1)
        ((minusp exponent) (reciprocal (rat-pow base (negate exponent))))
        ((zerop exponent) 1)
        ((< exponent 1) (%rat-pow base exponent))
        ((= exponent 1) base)
        ((and (typep exponent 'integer)
              (oddp exponent))
         (* base (rat-pow base (- exponent 1))))
        (t (rat-pow (square base) (/ exponent 2)))))

(defmethod expt ((base cl:rational) (exponent single-float))
  (rat-pow base (rational exponent)))

(defmethod expt ((base cl:rational) (exponent double-float))
  (rat-pow base (rational exponent)))

(defmethod expt ((base cl:rational) (exponent cl:rational))
  (rat-pow base exponent))

(defmethod expt ((base single-float) (exponent cl:rational))
  (rat-pow (rational base) exponent))

(defmethod expt ((base double-float) (exponent cl:rational))
  (rat-pow (rational base) exponent))

(defmethod expt ((base lft-stream) (exponent cl:rational))
  (exp-lft-stream (funcall (lft-multiply-by-rat exponent) (log base))))

(defmethod expt (base (exponent lft-stream))
  (exp-lft-stream (* (log base) exponent)))

(defmethod expt ((base (eql 0)) (exponent lft-stream))
  0)

(defmethod expt ((base (eql 1)) (exponent lft-stream))
  1)

(defmethod cbrt ((number lft-stream))
  (expt number 1/3))

(defmethod sqrt ((number lft-stream))
  (%sqrt-lft-stream number))

(defun tan-rat (rat)
  (cond ((< rat 0) (negate (tan-rat (- rat))))
        ((zerop rat) 0)
        ((<= rat 1) (%rat-tan rat))
        (t
         ;; (2 tan(x/2))/(1 - tan(x/2)^2)
         (let ((tan-half (tan-rat (/ rat 2))))
           (funcall (make-bilft 0  1 1 0
                                -1 0 0 1)
                    tan-half
                    tan-half)))))

(defun tan-lft-stream (lft-stream)
  (lft-stream-minus-p
   lft-stream
   (lambda (refined)
     (negate (tan-lft-stream (negate refined))))
   (lambda (refined)
     (lft-stream-less-than-rat
      refined 1
      #'%tan-lft-stream
      (lambda (refined*)
        (let ((tan-half-x (tan-lft-stream (/ refined* 2))))
          (funcall (make-bilft 0  1 1 0
                               -1 0 0 1)
                   tan-half-x
                   tan-half-x)))))))

(defmethod tan ((theta cl:rational))
  (tan-rat theta))

(defmethod tan ((theta lft-stream))
  (tan-lft-stream theta))

(defun cos (theta)
  (let ((tan-half-x (tan (/ theta 2))))
    (funcall (make-bilft -1 0 0 1
                         1 0 0 1) ; => #<BILFT (-xy + 1)/(xy + 1)>
             tan-half-x
             tan-half-x)))

(defun sin (theta)
  (let ((tan-half-x (tan (/ theta 2))))
    (funcall (make-bilft 0 1 1 0
                         1 0 0 1) ; => #<BILFT (x + y)/(xy + 1)>
             tan-half-x
             tan-half-x)))

(defun atan-rat (y x)
  (cond ((minusp y) (negate (atan-rat (- y) x)))
        ((minusp x) (- pi (atan-rat y (- x))))
        ((> y x) (- pi/2 (atan-rat x y)))
        (t (%rat-atan (/ y x)))))

(defun atan-lft-stream-rat (lft-stream-y x)
  (lft-stream-minus-p
   lft-stream-y
   (lambda (lft-stream-y*)
     (negate (atan-lft-stream-rat (negate lft-stream-y*) x)))
   (lambda (lft-stream-y*)
     (if (minusp x)
         (- pi (atan-lft-stream-rat lft-stream-y (- x)))
         (lft-stream-less-than-rat
          lft-stream-y* x
          (lambda (lft-stream-y**)
            (%atan-lft-stream (/ lft-stream-y** x)))
          (lambda (lft-stream-y**)
            (- pi/2 (atan-rat-lft-stream x lft-stream-y**))))))))

(defun atan-rat-lft-stream (y lft-stream-x)
  (if (minusp y)
      (negate (atan-rat-lft-stream (- y) lft-stream-x))
      (lft-stream-minus-p
       lft-stream-x
       (lambda (lft-stream-x*)
         (- pi (atan-rat-lft-stream y (negate lft-stream-x*))))
       (lambda (lft-stream-x*)
         (lft-stream-less-than-rat
          lft-stream-x* y
          (lambda (lft-stream-x**)
            (- pi/2 (atan-lft-stream-rat lft-stream-x** y)))
          (lambda (lft-stream-x**)
            (%atan-lft-stream (/ y lft-stream-x**))))))))
          
(defun atan-lft-stream (y x)
  (lft-stream-minus-p
   y
   (lambda (y*)
     (negate (atan-lft-stream (negate y*) x)))
   (lambda (y*)
     (lft-stream-minus-p
      x
      (lambda (x*)
        (- pi (atan-lft-stream y* (negate x*))))
      (lambda (x*)
        (lft-stream-minus-p
         (- y* x*)
         (lambda (_)
           (declare (ignore _))
           (%atan-lft-stream (/ y* x*)))
         (lambda (_)
           (declare (ignore _))
           (- pi/2 (atan-lft-stream x y)))))))))

(defmethod atan ((y lft-stream) &optional (x 1))
  (etypecase x
    (lft-stream (atan-lft-stream y x))
    (single-float (atan-lft-stream-rat y (rational x)))
    (double-float (atan-lft-stream-rat y (rational x)))
    (cl:rational (atan-lft-stream-rat y x))))

(defmethod atan ((y single-float) &optional (x 1))
  (etypecase x
    (lft-stream (atan-rat-lft-stream (rational y) x))
    (single-float (atan-rat (rational y) (rational x)))
    (double-float (atan-rat (rational y) (rational x)))
    (cl:rational (atan-rat (rational y) x))))

(defmethod atan ((y double-float) &optional (x 1))
  (etypecase x
    (lft-stream (atan-rat-lft-stream (rational y) x))
    (single-float (atan-rat (rational y) (rational x)))
    (double-float (atan-rat (rational y) (rational x)))
    (cl:rational (atan-rat (rational y) x))))

(defmethod atan ((y cl:rational) &optional (x 1))
  (etypecase x
    (lft-stream (atan-rat-lft-stream y x))
    (single-float (atan-rat y (rational x)))
    (double-float (atan-rat y (rational x)))
    (cl:rational (atan-rat y x))))

(defun arccos-rat (rational)
  (let ((sqrat (square rational)))
    (atan (sqrt (/ (- 1 sqrat) sqrat)))))

(defun arccos-lft-stream (lft-stream)
  (atan
   (sqrt
    (funcall
     (make-bilft -1 0 0 1
                 1  0 0 0) ; => (1 - x^2)/x^2
     lft-stream
     lft-stream))))

(defmethod acos ((x cl:rational))
  (arccos-rat x))

(defmethod acos ((x lft-stream))
  (arccos-lft-stream x))

(defun arcsin-rat (rational)
  (let ((sqrat (square rational)))
    (atan (sqrt (/ sqrat (- 1 sqrat))))))

(defun arcsin-lft-stream (lft-stream)
  (atan
   (sqrt
    (funcall (make-bilft 1  0 0 0
                         -1 0 0 1) ; => #<BILFT -xy/(xy - 1)>
             lft-stream
             lft-stream))))

(defmethod asin ((x cl:rational))
  (arcsin-rat x))

(defmethod asin ((x lft-stream))
  (arcsin-lft-stream x))

(defun log-base (number base)
  (/ (etypecase number
        ((or integer cl:rational) (log-rat number))
        (lft-stream (log-lft-stream number)))
      (etypecase base
        ((or integer cl:rational) (log-rat base))
        (lft-stream (log-lft-stream base)))))

(defun hyperbolic-function (bilft x)
  (let ((exp (exp x)))
    (funcall bilft exp exp)))

(defmethod sinh (x)
  (hyperbolic-function (make-bilft 1 0 0 -1
                                   0 1 1  0) ; => #<BILFT (xy - 1)/(x + y)>
                       x))

(defmethod cosh (x)
  (hyperbolic-function (make-bilft 1 0 0 1
                                   0 1 1 0) ; => #<BILFT (xy + 1)/(x + y)>
                       x))

(defmethod tanh (x)
  (hyperbolic-function (make-bilft 1 0 0 -1
                                   1 0 0 1) ; => #<BILFT (xy - 1)/(xy + 1)>
                       x))

(defmethod coth (x)
  (reciprocal (tanh x)))

(defmethod sech (x)
  (reciprocal (cosh x)))

(defmethod csch (x)
  (reciprocal (sinh x)))

(defun asinh-lft-stream (lft-stream)
  (log-lft-stream
   (+ lft-stream
       (sqrt
        (funcall (make-bilft 1 0 0 1
                             0 0 0 1)
                 lft-stream
                 lft-stream)))))

(defun asinh-rat (rat)
  (log-lft-stream
   (+ rat
       (sqrt (+ (* rat rat) 1)))))

(defmethod asinh ((number cl:rational))
  (asinh-rat number))

(defmethod asinh ((number single-float))
  (asinh-rat (rational number)))

(defmethod asinh ((number double-float))
  (asinh-rat (rational number)))

(defmethod asinh ((number lft-stream))
  (asinh-lft-stream number))

(defun acosh-lft-stream (lft-stream)
  (log-lft-stream
   (+ lft-stream
       (sqrt
        (funcall (make-bilft 1 0 0 -1
                             0 0 0 1)
                 lft-stream
                 lft-stream)))))

(defun acosh-rat (rat)
  (log-lft-stream
   (+ rat
       (sqrt (- (* rat rat) 1)))))

(defmethod acosh ((number cl:rational))
  (acosh-rat number))

(defmethod acosh ((number single-float))
  (acosh-rat (rational number)))

(defmethod acosh ((number double-float))
  (acosh-rat (rational number)))

(defmethod acosh ((number lft-stream))
  (acosh-lft-stream number))

(defmethod acosh ((number lft-stream))
  (acosh-lft-stream number))

(defun atanh-lft-stream (lft-stream)
  (funcall (lft-multiply-by-rat 1/2)
           (log-lft-stream
            (funcall lft-sign-infinity lft-stream))))

(defun atanh-rat (rat)
  (funcall (lft-multiply-by-rat 1/2)
           (log (funcall lft-sign-infinity rat))))

(defmethod atanh ((number cl:rational))
  (atanh-rat number))

(defmethod atanh ((number single-float))
  (atanh-rat (rational number)))

(defmethod atanh ((number double-float))
  (atanh-rat (rational number)))

(defmethod atanh ((number lft-stream))
  (atanh-lft-stream number))

;;;

(defparameter pythagoras (sqrt 2))
(defparameter theodorus (sqrt 3))
(defparameter sqrt-five (sqrt 5))

;; (defun approximate-zeta-three (nterms)
;;   (flet ((f (k)
;;             (* (expt -1 k)
;;                (/ (* (expt (factorial (+ (* 2 k) 1)) 3)
;;                      (expt (factorial (* 2 k)) 3)
;;                      (expt (factorial k) 3)
;;                      (+ (* 126392 k k k k k)
;;                         (* 412708 k k k k)
;;                         (* 531578 k k k)
;;                         (* 336367 k k)
;;                         (* 104000 k)
;;                         12463))
;;                   (* (factorial (+ (* 3 k) 2))
;;                      (expt (factorial (+ 4 k) 3)) 3))))))
;;     (/ (big-sigma #'f 0 nterms) 24)))

;; (defparameter apéry (limit-stream->cf-stream (stream-map #'approximate-zeta-three (integers))))
