;;; -*- Lisp -*-

(in-package "LINEAR-FRACTIONAL-TRANSFORM")

(defgeneric add2 (left right)
  (:method ((left rational) (right rational))
    (cl:+ left right))
  (:method ((left float) right)
    (add2 (rational left) right))
  (:method (left (right float))
    (add2 left (rational right))))

(defgeneric mul2 (left right)
  (:method ((left rational) (right rational))
    (cl:* left right))
  (:method ((left float) right)
    (mul2 (rational left) right))
  (:method (left (right float))
    (mul2 left (rational right))))

(defgeneric negate (number)
  (:method ((number (eql 'infinity)))
    'infinity)
  (:method ((number rational))
    (cl:- number))
  (:method ((number float))
    (negate (rational number))))

(defgeneric sub2 (left right)
  (:method (left right)
    (add2 left (negate right))))

(defgeneric reciprocal (number)
  (:method ((number (eql 'infinity)))
    0)
  (:method ((number (eql 0)))
    'infinity)
  (:method ((number rational))
    (cl:/ 1 number))
  (:method ((number float))
    (reciprocal (rational number))))

(defgeneric div2 (left right)
  (:method (left right)
    (mul2 left (reciprocal right))))

(defun x+ (&rest args)
  (fold-left #'add2 0 args))

(defun x* (&rest args)
  (fold-left #'mul2 1 args))

(defun x- (leftmost &rest rights)
  (cond ((consp rights) (sub2 leftmost (fold-left #'add2 (car rights) (cdr rights))))
        ((null rights)  (negate leftmost))
        (t (error "Unexpected value ~s." rights))))

(defun x/ (leftmost &rest rights)
  (cond ((consp rights) (div2 leftmost (fold-left #'mul2 (car rights) (cdr rights))))
        ((null rights)  (reciprocal leftmost))
        (t (error "Unexpected value ~s." rights))))

(defgeneric x-exp (number)
  (:method ((number float))
    (x-exp (rational number))))

(defgeneric x-log (number)
  (:method ((number float))
    (x-log (rational number))))

(defgeneric x-sqrt (number)
  (:method ((number float))
    (x-sqrt (rational number))))
