;;; -*- Lisp -*-

(in-package "LINEAR-FRACTIONAL-TRANSFORMATION")

;;; Arithmetic utilities

(defun 2x1v< (a
              b

              c
              d)
  "Returns TRUE if a/b < c/d"
  (< (- (* a d) (* b c)) 0))

(defun 2x2m<2x1v (a b
                  c d

                  e
                  f)
  "Returns TRUE if both a/c and b/d are less than e/f"
  (and (2x1v< a
              c

              e
              f)
       (2x1v< b
              d

              e
              f)))

(defun 2x2m<2x2m (a b
                  c d

                  e f
                  g h)
  "Returns TRUE if both a/c and b/d are less than each of e/g and f/h."
  (and (2x2m<2x1v a b
                  c d

                  e
                  g)
       (2x2m<2x1v a b
                  c d

                  f
                  h)))

(defun disjoint? (a b
                  c d

                  e f
                  g h)
  "Returns TRUE if range [a/c, b/d] is disjoint from range [e/g, f/h]"
  (or (2x2m<2x2m a b
                 c d

                 e f
                 g h)
      (2x2m<2x2m e f
                 g h

                 a b
                 c d)))

(defun tensor-disjoint? (a b c d
                         e f g h)
  (disjoint? a c
             e g

             b d
             f h))

(defun 2x1-vector-vector-multiply (a b

                                   c
                                   d)
  (declare (type integer a b c d))
  (+ (* a c) (* b d)))

(defun 2x2-matrix-vector-multiply (a b
                                   c d

                                   e
                                   f

                                   receiver)
  "(funcall receiver i j) where
[i]   [a b]   [e]
[j] = [c d] x [f]
"
  (declare (type integer a b c d e f))
  (funcall receiver
           (2x1-vector-vector-multiply a b e f)
           (2x1-vector-vector-multiply c d e f)))

(defun 2x2-matrix-multiply (a b
                            c d

                            e f
                            g h

                            receiver)
  (declare (type integer a b c d e f g h))
  "(funcall receiver i j k l) where
   [i j]   [a b]   [e f]
   [k l] = [c d] x [g h]
"
  (2x2-matrix-vector-multiply a b
                              c d

                              e
                              g
                              (lambda (i k)
                                (declare (integer i k))
                                (2x2-matrix-vector-multiply a b
                                                            c d

                                                            f
                                                            h
                                                            (lambda (j l)
                                                              (declare (integer j l))
                                                              (funcall receiver i j k l))))))

  ;; (funcall receiver
  ;;          (+ (* a e) (* b g)) (+ (* a f) (* b h))
  ;;          (+ (* c e) (* d g)) (+ (* c f) (* d h))))

(defun 2x2-matrix-tame-inverse (a b
                                c d
                                receiver)
  ;;
  ;; Inverse is 1/(- (* a d) (* b c)) * [ d -b]
  ;;                                    [-c  a]
  ;; Tame inverse omits the reciprocal of the determinant.
  "(funcall receiver i j k l) where
   [i j]   [ d -b]
   [k l] = [-c  a]"
  (declare (integer a b c d))
  (funcall receiver
           d    (- b)
           (- c)   a))

(defun 2mx2t-multiply (a b
                       c d

                       e f g h
                       i j k l

                       receiver)
  (declare (integer a b c d e f g h i j k l))
  ;; Multiply a 2x2 matrix by a 2x2x2 tensor yielding a 2x2x2 tensor
  (funcall receiver
           (+ (* a e) (* b i)) (+ (* a f) (* b j)) (+ (* a g) (* b k)) (+ (* a h) (* b l))
           (+ (* c e) (* d i)) (+ (* c f) (* d j)) (+ (* c g) (* d k)) (+ (* c h) (* d l))))

(defun 2tx2m-multiply-x (a b c d
                         e f g h

                         i j
                         k l

                         receiver)
  (declare (integer a b c d e f g h i j k l))
  ;; Multiply a 2x2x2 tensor by a 2x2 matrix in the x direction yielding a 2x2x2 tensor
  (funcall receiver
           (+ (* a i) (* c k)) (+ (* b i) (* d k)) (+ (* a j) (* c l)) (+ (* b j) (* d l))
           (+ (* e i) (* g k)) (+ (* f i) (* h k)) (+ (* e j) (* g l)) (+ (* f j) (* h l))))

(defun 2tx2m-multiply-y (a b c d
                         e f g h

                         i j
                         k l

                         receiver)
  ;; Multiply a 2x2x2 tensor by a 2x2 matrix in the y direction yielding a 2x2x2 tensor
  (declare (integer a b c d e f g h i j k l))
  (funcall receiver
           (+ (* a i) (* b k)) (+ (* a j) (* b l)) (+ (* c i) (* d k)) (+ (* c j) (* d l))
           (+ (* e i) (* f k)) (+ (* e j) (* f l)) (+ (* g i) (* h k)) (+ (* g j) (* h l))))

;;;           A
;;;      <0  =0  >0
;;;   <0 -1  -1   0
;;; B =0 -1   0   1
;;;   >0  0   1   1

(defun s (a b)
  (cond ((minusp a) (if (plusp b)
                           0
                           -1))
        ((plusp a) (if (minusp b)
                           0
                           1))
        ((minusp b) -1)
        ((plusp b) 1)
        (t 0)))

;; (defun sign (a b)
;;   (cond ((< a 0) (if (<= b 0)
;;                      -1
;;                      0))
;;         ((= a 0) (cond ((< b 0) -1)
;;                        ((= b 0) 0)
;;                        (t 1)))
;;         (t (if (< b 0)
;;                0
;;                1))))

(defun 2x2-matrix-range? (a b
                          c d)
  (let ((l (s a c))
        (r (s b d)))
    (and (= l r)
         (not (= r 0)))))

(defun 2x2x2-tensor-range? (a b c d
                            e f g h)
  (let ((i (s a e))
        (j (s b f))
        (k (s c g))
        (l (s d h)))
    (and (= i j)
         (= j k)
         (= k l)
         (not (= l 0)))))

;;;
;;; Evaluators for homographic functions
;;;

(defun evaluate-lft (a b
                     c d numerator denominator)
  "Return lim x->n (ax + b)/(cx + d)"
  (check-type a integer)
  (check-type b integer)
  (check-type c integer)
  (check-type d integer)
  (check-type numerator integer)
  (check-type denominator (integer 0 *))
  (let ((result-numerator   (+ (* a numerator) (* b denominator)))
        (result-denominator (+ (* c numerator) (* d denominator))))
    (cond ((not (zerop result-denominator)) (/ result-numerator result-denominator))
          ((= (* a d) (* b c)) (cond ((not (zerop c)) (/ a c))
                                     ((not (zerop d)) (/ b d))
                                     ((or (not (zerop a)) (not (zerop b))) 'infinity)
                                     (t (error 'division-by-zero))))
          (t 'infinity))))

(defun evaluate-bilft (a b c d
                       e f g h m-num m-den n-num n-den)
  "Return lim x->m, y->n (axy + bx + cy + d)/(exy + fx + gy + h)"
  (check-type a integer)
  (check-type b integer)
  (check-type c integer)
  (check-type d integer)
  (check-type e integer)
  (check-type f integer)
  (check-type g integer)
  (check-type h integer)
  (check-type m-num integer)
  (check-type m-den (integer 0 *))
  (check-type n-num integer)
  (check-type n-den (integer 0 *))
  (let ((result-numerator   (+ (* a m-num n-num) (* b m-num n-den) (* c m-den n-num) (* d m-den n-den)))
        (result-denominator (+ (* e m-num n-num) (* f m-num n-den) (* g m-den n-num) (* h m-den n-den))))
    (cond ((not (zerop result-denominator)) (/ result-numerator result-denominator))
          ((not (zerop result-numerator)) 'infinity)
          ((= (* a f g h)
              (* b e g h)
              (* c e f h)
              (* d e f g)) (cond ((not (zerop h)) (/ d h))
                                 ((not (zerop g)) (/ c g))
                                 ((not (zerop f)) (/ b f))
                                 ((not (zerop e)) (/ a e))
                                 ((or (not (zerop a))
                                      (not (zerop b))
                                      (not (zerop c))
                                      (not (zerop d))) 'infinity)
                                 (t (error 'division-by-zero))))
          (t 'infinity))))

;;;
;;; Formatters for lft equations
;;;

(defun format-term (stream coefficient variables &key initial-term)
  "Format a term like ' + 3x' or ' - 2y'."
  (check-type coefficient integer)
  (check-type variables (or null string))
  ;; If initial-term is true, then omit spaces around the sign and suppress
  ;; positive signs.  Used for printing initial terms.
  (unless (zerop coefficient)
    (format stream "~:[ ~:[+~;-~] ~;~:[~;-~]~]~[~;~:;~:*~d~]~:[~;1~]~@[~a~]"
            initial-term
            (minusp coefficient)
            (abs coefficient)
            (and (null variables) (= (abs coefficient) 1))
            variables)))

(defun format-term-list (stream coefficients variables &key suppress-parens)
  "Format a list of terms given a list of coefficients and variables."
  (labels ((l1 (coefficients variables)
             (cond ((null coefficients) (format stream "0"))
                   ((zerop (car coefficients)) (l1 (cdr coefficients) (cdr variables)))
                   ((or suppress-parens
                        (every #'zerop (cdr coefficients)))
                    (l2 coefficients variables t))
                   (t (format stream "(")
                      (l2 coefficients variables t)
                      (format stream ")"))))

           (l2 (coefficients variables initial-term)
             (unless (null coefficients)
               (format-term stream (car coefficients) (car variables) :initial-term initial-term)
               (l2 (cdr coefficients) (cdr variables) nil))))

    (l1 coefficients variables)))

(defparameter *print-lambdas* nil
  "If T, print λx. in the algebraic printout of a LFT.")

(defun format-lft-equation (stream a b
                                   c d)
  (when *print-lambdas* (format stream "λx."))
  (cond ((and (zerop c)
              (= d 1))
         (format-term-list stream (list a b) '("x" nil) :suppress-parens t))
        ((and (= c 1)
              (zerop d))
         (if (zerop a)
             (format stream "~d/x" b)
             (format stream "~d ~:[+~;-~] ~d/x" a (minusp b) (abs b))))
        (t
         (format-term-list stream (list a b) '("x" nil))
         (format stream "/")
         (format-term-list stream (list c d) '("x" nil)))))

(defun format-bilft-equation (stream a b c d
                                     e f g h)
  (when *print-lambdas* (format stream "λxy."))
  (if (and (zerop e)
           (zerop f)
           (zerop g)
           (= h 1))
      (format-term-list stream (list a b c d) '("xy" "x" "y" nil) :suppress-parens t)
      (progn
        (format-term-list stream (list a b c d) '("xy" "x" "y" nil))
        (format stream "/")
        (format-term-list stream (list e f g h) '("xy" "x" "y" nil)))))

;;;
;;; Instances of lft functions
;;;

(eval-when (:compile-toplevel :load-toplevel :execute)
(defclass lft ()
  ((a :initarg :a
      :initform (error "Required initarg :a omitted")
      :type integer)
   (b :initarg :b
      :initform (error "Required initarg :b omitted")
      :type integer)
   (c :initarg :c
      :initform (error "Required initarg :c omitted")
      :type integer)
   (d :initarg :d
      :initform (error "Required initarg :d omitted")
      :type integer))
  (:metaclass sb-mop:funcallable-standard-class))

(defgeneric funcall-lft (lft arg)
  (:method (lft (arg (eql 'infinity))) (evaluate-lft
                                        (slot-value lft 'a)
                                        (slot-value lft 'b)
                                        (slot-value lft 'c)
                                        (slot-value lft 'd)
                                        1 0))
  (:method (lft (arg integer)) (evaluate-lft
                                (slot-value lft 'a)
                                (slot-value lft 'b)
                                (slot-value lft 'c)
                                (slot-value lft 'd)
                                arg 1))
  (:method (lft (arg rational)) (evaluate-lft
                                 (slot-value lft 'a)
                                 (slot-value lft 'b)
                                 (slot-value lft 'c)
                                 (slot-value lft 'd)
                                 (numerator arg) (denominator arg)))
  (:method (lft (arg float)) (funcall-lft lft (rational arg)))
  )

(defmethod initialize-instance :after ((instance lft) &rest initargs)
  (declare (ignore initargs))
  (sb-mop:set-funcallable-instance-function instance (lambda (arg) (funcall-lft instance arg))))
)

(defmethod print-object ((instance lft) stream)
  (print-unreadable-object (instance stream :type t)
    (format-lft-equation
     stream
     (slot-value instance 'a)
     (slot-value instance 'b)
     (slot-value instance 'c)
     (slot-value instance 'd))))

(eval-when (:compile-toplevel :load-toplevel :execute)
(defclass bilft ()
  ((a :initarg :a
      :initform (error "Required initarg :a omitted")
      :type integer)
   (b :initarg :b
      :initform (error "Required initarg :b omitted")
      :type integer)
   (c :initarg :c
      :initform (error "Required initarg :c omitted")
      :type integer)
   (d :initarg :d
      :initform (error "Required initarg :d omitted")
      :type integer)
   (e :initarg :e
      :initform (error "Required initarg :e omitted")
      :type integer)
   (f :initarg :f
      :initform (error "Required initarg :f omitted")
      :type integer)
   (g :initarg :g
      :initform (error "Required initarg :g omitted")
      :type integer)
   (h :initarg :h
      :initform (error "Required initarg :h omitted")
      :type integer))
  (:metaclass sb-mop:funcallable-standard-class))

(defgeneric funcall-bilft (bilft x y)
  (:method (bilft (x (eql 'infinity)) (y (eql 'infinity)))
    (evaluate-bilft
     (slot-value bilft 'a) (slot-value bilft 'b) (slot-value bilft 'c) (slot-value bilft 'd)
     (slot-value bilft 'e) (slot-value bilft 'f) (slot-value bilft 'g) (slot-value bilft 'h)
     1 0 1 0))
  (:method (bilft (x (eql 'infinity)) (y integer))
    (evaluate-bilft
     (slot-value bilft 'a) (slot-value bilft 'b) (slot-value bilft 'c) (slot-value bilft 'd)
     (slot-value bilft 'e) (slot-value bilft 'f) (slot-value bilft 'g) (slot-value bilft 'h)
     1 0 y 1))
  (:method (bilft (x (eql 'infinity)) (y rational))
    (evaluate-bilft
     (slot-value bilft 'a) (slot-value bilft 'b) (slot-value bilft 'c) (slot-value bilft 'd)
     (slot-value bilft 'e) (slot-value bilft 'f) (slot-value bilft 'g) (slot-value bilft 'h)
     1 0 (numerator y) (denominator y)))
  (:method (bilft (x integer) (y (eql 'infinity)))
    (evaluate-bilft
     (slot-value bilft 'a) (slot-value bilft 'b) (slot-value bilft 'c) (slot-value bilft 'd)
     (slot-value bilft 'e) (slot-value bilft 'f) (slot-value bilft 'g) (slot-value bilft 'h)
     x 1 1 0))
  (:method (bilft (x integer) (y integer))
    (evaluate-bilft
     (slot-value bilft 'a) (slot-value bilft 'b) (slot-value bilft 'c) (slot-value bilft 'd)
     (slot-value bilft 'e) (slot-value bilft 'f) (slot-value bilft 'g) (slot-value bilft 'h)
     x 1 y 1))
  (:method (bilft (x integer) (y rational))
    (evaluate-bilft
     (slot-value bilft 'a) (slot-value bilft 'b) (slot-value bilft 'c) (slot-value bilft 'd)
     (slot-value bilft 'e) (slot-value bilft 'f) (slot-value bilft 'g) (slot-value bilft 'h)
     x 1 (numerator y) (denominator y)))
  (:method (bilft (x rational) (y (eql 'infinity)))
    (evaluate-bilft
     (slot-value bilft 'a) (slot-value bilft 'b) (slot-value bilft 'c) (slot-value bilft 'd)
     (slot-value bilft 'e) (slot-value bilft 'f) (slot-value bilft 'g) (slot-value bilft 'h)
     (numerator x) (denominator x) 1 0))
  (:method (bilft (x rational) (y integer))
    (evaluate-bilft
     (slot-value bilft 'a) (slot-value bilft 'b) (slot-value bilft 'c) (slot-value bilft 'd)
     (slot-value bilft 'e) (slot-value bilft 'f) (slot-value bilft 'g) (slot-value bilft 'h)
     (numerator x) (denominator x) y 1))
  (:method (bilft (x rational) (y rational))
    (evaluate-bilft
     (slot-value bilft 'a) (slot-value bilft 'b) (slot-value bilft 'c) (slot-value bilft 'd)
     (slot-value bilft 'e) (slot-value bilft 'f) (slot-value bilft 'g) (slot-value bilft 'h)
     (numerator x) (denominator x) (numerator y) (denominator y)))
  (:method (bilft (x float) y)
    (funcall bilft (rational x) y))
  (:method (bilft x (y float))
    (funcall bilft x (rational y))))

(defmethod initialize-instance :after ((instance bilft) &rest initargs)
  (declare (ignore initargs))
  (sb-mop:set-funcallable-instance-function instance (lambda (x y) (funcall-bilft instance x y))))
)

(defmethod print-object ((instance bilft) stream)
  (print-unreadable-object (instance stream :type t)
    (format-bilft-equation
     stream
     (slot-value instance 'a)
     (slot-value instance 'b)
     (slot-value instance 'c)
     (slot-value instance 'd)
     (slot-value instance 'e)
     (slot-value instance 'f)
     (slot-value instance 'g)
     (slot-value instance 'h))))

;;;
;;; Constructors
;;;

(eval-when (:compile-toplevel :load-toplevel :execute)

(defparameter *reduce-to-lowest* nil
  "If T, use the GCD to reduce the coefficients of a LFT to lowest terms.")

(defun canonicalize-lft-coefficients (a b
                                      c d receiver)
  (cond ((or (minusp c)
             (and (zerop c)
                  (minusp d)))
         (canonicalize-lft-coefficients (- a) (- b)
                                        (- c) (- d) receiver))
        (*reduce-to-lowest*
         (let ((gcd (gcd a b c d)))
           (if (> gcd 1)
               (funcall receiver
                        (/ a gcd) (/ b gcd)
                        (/ c gcd) (/ d gcd))
               (funcall receiver
                        a b
                        c d))))
        ((and (evenp a)
              (evenp b)
              (evenp c)
              (evenp d))
         (canonicalize-lft-coefficients (/ a 2) (/ b 2)
                                        (/ c 2) (/ d 2) receiver))
        (t (funcall receiver
                    a b
                    c d))))

(defun %make-lft (a b c d)
  (canonicalize-lft-coefficients
   a b
   c d
   (lambda (a* b*
            c* d*)
     (make-instance 'lft
                    :a a* :b b*
                    :c c* :d d*))))

(defun make-lft (a b c d)
  (etypecase a
    (float (make-lft (rational a) b c d))
    (integer
     (etypecase b
       (float (make-lft a (rational b) c d))
       (integer
        (etypecase c
          (float (make-lft a b (rational c) d))
          (integer
           (etypecase d
             (float (make-lft a b c (rational d)))
             (integer (%make-lft a b
                                 c d))
             (rational (make-lft (* a (denominator d)) (* b (denominator d))
                                 (* c (denominator d)) (numerator d)))))
          (rational (make-lft (* a (denominator c)) (* b (denominator c))
                              (numerator c)         (* d (denominator c))))))
       (rational (make-lft (* a (denominator b)) (numerator b)
                           (* c (denominator b)) (* d (denominator b))))))
    (rational (make-lft (numerator a)         (* b (denominator a))
                        (* c (denominator a)) (* d (denominator a))))))

) ;; eval-when

(defparameter lft-identity
  (if (boundp 'lft-identity)
      lft-identity
      (make-lft 1 0
                0 1)))

(defparameter lft-negate
  (if (boundp 'lft-negate)
      lft-negate
      (make-lft -1 0
                0  1)))

(defmethod negate ((lft lft))
  (funcall lft-negate lft))

(defparameter lft-reciprocal
  (if (boundp 'lft-reciprocal)
      lft-reciprocal
      (make-lft 0 1
                1 0)))

(defmethod reciprocal ((lft lft))
  (funcall lft-reciprocal lft))

(eval-when (:compile-toplevel :load-toplevel :execute)

(defun canonicalize-bilft-coefficients (a b c d e f g h receiver)
  (cond ((or (minusp e)
             (and (zerop e)
                  (or (minusp f)
                      (and (zerop f)
                           (or (minusp g)
                               (and (zerop g)
                                    (minusp h)))))))
         (canonicalize-bilft-coefficients (- a) (- b) (- c) (- d)
                                          (- e) (- f) (- g) (- h) receiver))
        (*reduce-to-lowest*
         (let ((gcd (gcd a b c d e f g h)))
           (if (> gcd 1)
               (funcall receiver
                        (/ a gcd) (/ b gcd) (/ c gcd) (/ d gcd)
                        (/ e gcd) (/ f gcd) (/ g gcd) (/ h gcd))
               (funcall receiver
                        a b c d
                        e f g h))))
        ((and (evenp a)
              (evenp b)
              (evenp c)
              (evenp d)
              (evenp e)
              (evenp f)
              (evenp g)
              (evenp h))
         (canonicalize-bilft-coefficients (/ a 2) (/ b 2) (/ c 2) (/ d 2)
                                          (/ e 2) (/ f 2) (/ g 2) (/ h 2) receiver))
        (t (funcall receiver
                    a b c d
                    e f g h))))

(defun %make-bilft (a b c d
                    e f g h)
  (canonicalize-bilft-coefficients
   a b c d
   e f g h
   (lambda (a* b* c* d*
            e* f* g* h*)
     (make-instance 'bilft
                    :a a* :b b* :c c* :d d*
                    :e e* :f f* :g g* :h h*))))

(defun make-bilft (a b c d
                   e f g h)
  (etypecase a
    (float (make-bilft (rational a) b c d e f g h))
    (integer
     (etypecase b
       (float (make-bilft a (rational b) c d e f g h))
       (integer
        (etypecase c
          (float (make-bilft a b (rational c) d e f g h))
          (integer
           (etypecase d
             (float (make-bilft a b c (rational d) e f g h))
             (integer
              (etypecase e
                (float (make-bilft a b c d (rational e) f g h))
                (integer
                 (etypecase f
                   (float (make-bilft a b c d e (rational f) g h))
                   (integer
                    (etypecase g
                      (float (make-bilft a b c d e f (rational g) h))
                      (integer
                       (etypecase h
                         (float (make-bilft a b c d e f g (rational h)))
                         (integer
                          (%make-bilft a b c d
                                       e f g h))
                         (rational
                          (make-bilft
                           (* a (denominator h)) (* b (denominator h)) (* c (denominator h)) (* d (denominator h))
                           (* e (denominator h)) (* f (denominator h)) (* g (denominator h)) (numerator h)))))
                      (rational
                       (make-bilft
                        (* a (denominator g)) (* b (denominator g)) (* c (denominator g)) (* d (denominator g))
                        (* e (denominator g)) (* f (denominator g)) (numerator g)         (* h (denominator g))))))
                   (rational
                    (make-bilft
                     (* a (denominator f)) (* b (denominator f)) (* c (denominator f)) (* d (denominator f))
                     (* e (denominator f)) (numerator f)         (* g (denominator f)) (* h (denominator f))))))
                (rational
                 (make-bilft
                  (* a (denominator e)) (* b (denominator e)) (* c (denominator e)) (* d (denominator e))
                  (numerator e)         (* f (denominator e)) (* g (denominator e)) (* h (denominator e))))))
             (rational
              (make-bilft
               (* a (denominator d)) (* b (denominator d)) (* c (denominator d)) (numerator d)
               (* e (denominator d)) (* f (denominator d)) (* g (denominator d)) (* h (denominator d))))))
          (rational
           (make-bilft
            (* a (denominator c)) (* b (denominator c)) (numerator c) (* d (denominator c))
            (* e (denominator c)) (* f (denominator c)) (* g (denominator c)) (* h (denominator c))))))
       (rational
        (make-bilft
         (* a (denominator b)) (numerator b)         (* c (denominator b)) (* d (denominator b))
         (* e (denominator b)) (* f (denominator b)) (* g (denominator b)) (* h (denominator b))))))
    (rational
     (make-bilft
      (numerator a)         (* b (denominator a)) (* c (denominator a)) (* d (denominator a))
      (* e (denominator a)) (* f (denominator a)) (* g (denominator a)) (* h (denominator a))))))
) ;; End EVAL-WHEN

;;;
;;; Functions that need access to the coefficients
;;;

(defun inverse-lft (hf)
  (check-type hf lft)
  (let ((a (slot-value hf 'a))
        (b (slot-value hf 'b))
        (c (slot-value hf 'c))
        (d (slot-value hf 'd)))
    (2x2-matrix-tame-inverse a b
                             c d
                             #'%make-lft)))

(defmethod inverse ((function lft))
  (inverse-lft function))

(defun compose-lft-lft (left right)
  (check-type left lft)
  (check-type right lft)
  (let ((a (slot-value left 'a))
        (b (slot-value left 'b))
        (c (slot-value left 'c))
        (d (slot-value left 'd))

        (e (slot-value right 'a))
        (f (slot-value right 'b))
        (g (slot-value right 'c))
        (h (slot-value right 'd)))
    (2x2-matrix-multiply a b
                         c d

                         e f
                         g h
                         #'%make-lft)))

(defmethod compose2 ((left lft) (right lft))
  (compose-lft-lft left right))

(defmethod funcall-lft (lft (arg lft))
  (compose-lft-lft lft arg))

(defun transpose-bilft (bihf)
  (check-type bihf bilft)
  (let ((a (slot-value bihf 'a))
        (b (slot-value bihf 'b))
        (c (slot-value bihf 'c))
        (d (slot-value bihf 'd))
        (e (slot-value bihf 'e))
        (f (slot-value bihf 'f))
        (g (slot-value bihf 'g))
        (h (slot-value bihf 'h)))
    (make-instance 'bilft
                   :a a :b c :c b :d d
                   :e e :f g :g f :h h)))

(defun compose-lft-bilft (left right)
  (check-type left lft)
  (check-type right bilft)
  (let ((a (slot-value left 'a))
        (b (slot-value left 'b))
        (c (slot-value left 'c))
        (d (slot-value left 'd))

        (i (slot-value right 'a))
        (j (slot-value right 'b))
        (k (slot-value right 'c))
        (l (slot-value right 'd))
        (m (slot-value right 'e))
        (n (slot-value right 'f))
        (o (slot-value right 'g))
        (p (slot-value right 'h)))
    (2mx2t-multiply a b
                    c d

                    i j k l
                    m n o p
                    #'%make-bilft)))

(defmethod funcall-lft (lft (arg bilft))
  (compose-lft-bilft lft arg))

(defun compose-bilft-lft-x (left right)
  (check-type left bilft)
  (check-type right lft)
  (let ((a (slot-value left 'a))
        (b (slot-value left 'b))
        (c (slot-value left 'c))
        (d (slot-value left 'd))
        (e (slot-value left 'e))
        (f (slot-value left 'f))
        (g (slot-value left 'g))
        (h (slot-value left 'h))

        (i (slot-value right 'a))
        (j (slot-value right 'b))
        (k (slot-value right 'c))
        (l (slot-value right 'd)))
    (2tx2m-multiply-x a b c d
                      e f g h

                      i j
                      k l

                      #'%make-bilft)))

(defun compose-bilft-lft-y (left right)
  (check-type left bilft)
  (check-type right lft)
  (let ((a (slot-value left 'a))
        (b (slot-value left 'b))
        (c (slot-value left 'c))
        (d (slot-value left 'd))
        (e (slot-value left 'e))
        (f (slot-value left 'f))
        (g (slot-value left 'g))
        (h (slot-value left 'h))

        (i (slot-value right 'a))
        (j (slot-value right 'b))
        (k (slot-value right 'c))
        (l (slot-value right 'd)))
    (2tx2m-multiply-y a b c d
                      e f g h

                      i j
                      k l

                      #'%make-bilft)))

(defgeneric range? (object)
  (:method ((object lft))
    (2x2-matrix-range? (slot-value object 'a) (slot-value object 'b)
                       (slot-value object 'c) (slot-value object 'd)))
  (:method ((object bilft))
    (2x2x2-tensor-range?
     (slot-value object 'a) (slot-value object 'b) (slot-value object 'c) (slot-value object 'd)
     (slot-value object 'e) (slot-value object 'f) (slot-value object 'g) (slot-value object 'h))))

(defun pole-negative-p (lft)
  (check-type lft lft)
  (if (plusp (slot-value lft 'c))
      (not (minusp (slot-value lft 'd)))
      ;; by canonicalization, c cannot be negative, so c must be zero
      (plusp (slot-value lft 'd))))

(defun poles-negative-p (bilft)
  (check-type bilft bilft)
  (and (plusp (slot-value bilft 'e))
       (not (minusp (slot-value bilft 'f)))
       (not (minusp (slot-value bilft 'g)))
       (not (minusp (slot-value bilft 'h)))))

(defun bilft-disjoint? (bilft)
  (tensor-disjoint?
   (slot-value bilft 'a) (slot-value bilft 'b) (slot-value bilft 'c) (slot-value bilft 'd)
   (slot-value bilft 'e) (slot-value bilft 'f) (slot-value bilft 'g) (slot-value bilft 'h)))

;;;
;;; Workhorse functions
;;;

(defun lft-zero-p (lft if-always-zero if-never-zero if-sometimes-zero)
  (check-type lft lft)
  (cond ((and (zerop (slot-value lft 'a))
              (zerop (slot-value lft 'b))
              (or (plusp (slot-value lft 'c))
                  (plusp (slot-value lft 'd))))
         (funcall if-always-zero))
        ((and (not (minusp (slot-value lft 'd)))
              (or (and (plusp (slot-value lft 'a))
                       (plusp (slot-value lft 'b)))
                  (and (minusp (slot-value lft 'a))
                       (minusp (slot-value lft 'b)))))
         (funcall if-never-zero))
        (t (funcall if-sometimes-zero))))

(defun lft-minus-p (lft if-always-minus if-never-minus if-unknown)
  (check-type lft lft)
  (cond ((minusp (slot-value lft 'd)) (funcall if-unknown))
        ((and (minusp (slot-value lft 'a))
              (minusp (slot-value lft 'b)))
         (funcall if-always-minus))
        ((and (not (minusp (slot-value lft 'a)))
              (not (minusp (slot-value lft 'b))))
         (funcall if-never-minus))
        (t (funcall if-unknown))))

(defun lft-plus-p (lft if-always-plus if-never-plus if-unknown)
  (check-type lft lft)
  (cond ((minusp (slot-value lft 'd)) (funcall if-unknown))
        ((and (plusp (slot-value lft 'a))
              (plusp (slot-value lft 'b)))
         (funcall if-always-plus))
        ((and (not (plusp (slot-value lft 'a)))
              (not (plusp (slot-value lft 'b))))
         (funcall if-never-plus))
        (t (funcall if-unknown))))

(defun lft-non-negative-p (lft if-non-negative if-never-non-negative if-unknown)
  (check-type lft lft)
  (cond ((minusp (slot-value lft 'd)) (funcall if-unknown))
        ((and (not (minusp (slot-value lft 'a)))
              (not (minusp (slot-value lft 'b))))
         (funcall if-non-negative))
        ((and (minusp (slot-value lft 'a))
              (minusp (slot-value lft 'b)))
         (funcall if-never-non-negative))
        (t (funcall if-unknown))))

(defun lft-less-than-rat (lft rat if-less-than if-not-less-than if-unknown)
  (check-type lft lft)
  (cond ((minusp (slot-value lft 'd)) (funcall if-unknown))
        ((and (plusp (slot-value lft 'c))
              ;; a/c < n/d , ad < nc
              (< (* (slot-value lft 'a) (denominator rat))
                 (* (slot-value lft 'c) (numerator rat)))
              (plusp (slot-value lft 'd))
              (< (* (slot-value lft 'b) (denominator rat))
                 (* (slot-value lft 'd) (numerator rat))))
         (funcall if-less-than))
        ((and (not (< (* (slot-value lft 'a) (denominator rat))
                      (* (slot-value lft 'c) (numerator rat))))
              (not (< (* (slot-value lft 'b) (denominator rat))
                      (* (slot-value lft 'd) (numerator rat)))))
         (funcall if-not-less-than))
        (t (funcall if-unknown))))

(defun lft-greater-than-rat (lft rat if-greater-than if-not-greater-than if-unknown)
  (check-type lft lft)
  (cond ((minusp (slot-value lft 'd)) (funcall if-unknown))
        ((and (> (* (slot-value lft 'a) (denominator rat))
                 (* (slot-value lft 'c) (numerator rat)))
              (> (* (slot-value lft 'b) (denominator rat))
                 (* (slot-value lft 'd) (numerator rat))))
         (funcall if-greater-than))
        ((and (<= (* (slot-value lft 'a) (denominator rat))
                 (* (slot-value lft 'c) (numerator rat)))
              (<= (* (slot-value lft 'b) (denominator rat))
                 (* (slot-value lft 'd) (numerator rat))))
         (funcall if-not-greater-than))
        (t (funcall if-unknown))))

(defun lft-add-rat (rat)
  (make-lft 1 rat
            0 1))

(defun lft-subtract-rat (rat)
  (make-lft 1 (- rat)
            0 1))

(defun lft-multiply-by-rat (rat)
  (make-lft rat 0
            0   1))

(defun lft-divide-by-rat (rat)
  (make-lft 1 0
            0 rat))

(defmethod add2 ((left lft) (right rational))
  (funcall (lft-add-rat right) left))

(defmethod add2 ((left rational) (right lft))
  (funcall (lft-add-rat left) right))

(defmethod mul2 ((left lft) (right rational))
  (funcall (lft-multiply-by-rat right) left))

(defmethod mul2 ((left rational) (right lft))
  (funcall (lft-multiply-by-rat left) right))

(defun lft-truncate (lft if-success if-failure)
  (check-type lft lft)
  (if (and (plusp (slot-value lft 'c))
           (plusp (slot-value lft 'd)))
      (let ((t1 (truncate (slot-value lft 'a) (slot-value lft 'c)))
            (t2 (truncate (slot-value lft 'b) (slot-value lft 'd))))
        (if (= t1 t2)
            (funcall if-success t1 (x- lft t1))
            (funcall if-failure)))
      (funcall if-failure)))

(defun lft->single (lft if-success if-failure)
  (check-type lft lft)
  (if (and (plusp (slot-value lft 'c))
           (plusp (slot-value lft 'd)))
      (let ((s1 (coerce (/ (slot-value lft 'a) (slot-value lft 'c)) 'single-float))
            (s2 (coerce (/ (slot-value lft 'b) (slot-value lft 'd)) 'single-float)))
        (if (= s1 s2)
            (funcall if-success s1)
            (funcall if-failure)))
      (funcall if-failure)))

(defun lft->double (lft if-success if-failure)
  (check-type lft lft)
  (if (and (plusp (slot-value lft 'c))
           (plusp (slot-value lft 'd)))
      (let ((d1 (coerce (/ (slot-value lft 'a) (slot-value lft 'c)) 'double-float))
            (d2 (coerce (/ (slot-value lft 'b) (slot-value lft 'd)) 'double-float)))
        (if (= d1 d2)
            (funcall if-success d1)
            (funcall if-failure)))
      (funcall if-failure)))

(defun bilft-minus-p (bilft if-minus if-plus if-unknown)
  (check-type bilft bilft)
  (if (poles-negative-p bilft)
      (let ((bound-a (funcall bilft 0 0))
            (bound-b (funcall bilft 0 'infinity))
            (bound-c (funcall bilft 'infinity 0))
            (bound-d (funcall bilft 'infinity 'infinity)))
        (if (and (numberp bound-a)
                 (numberp bound-b)
                 (numberp bound-c)
                 (numberp bound-d))
            (cond ((and (or (minusp bound-a) (minusp bound-b) (minusp bound-c) (minusp bound-d))
                        (not (plusp bound-a))
                        (not (plusp bound-b))
                        (not (plusp bound-c))
                        (not (plusp bound-d)))
                   (funcall if-minus))
                  ((and (or (plusp bound-a) (plusp bound-b) (plusp bound-c) (plusp bound-d))
                        (not (minusp bound-a))
                        (not (minusp bound-b))
                        (not (minusp bound-c))
                        (not (minusp bound-d)))
                   (funcall if-plus))
                  (t (funcall if-unknown)))
            (funcall if-unknown)))
      (funcall if-unknown)))

(defun bilft-zero-p (bilft if-zero if-not-zero if-unknown)
  (check-type bilft bilft)
  (if (poles-negative-p bilft)
      (let ((bound-a (funcall bilft 0 0))
            (bound-b (funcall bilft 0 'infinity))
            (bound-c (funcall bilft 'infinity 0))
            (bound-d (funcall bilft 'infinity 'infinity)))
        (cond ((and (numberp bound-a) (zerop bound-a)
                    (numberp bound-b) (zerop bound-b)
                    (numberp bound-c) (zerop bound-c)
                    (numberp bound-d) (zerop bound-d))
               (funcall if-zero))
              ((or (and (numberp bound-a) (plusp bound-a)
                        (numberp bound-b) (plusp bound-b)
                        (numberp bound-c) (plusp bound-c)
                        (numberp bound-d) (plusp bound-d))
                   (and (numberp bound-a) (minusp bound-a)
                        (numberp bound-b) (minusp bound-b)
                        (numberp bound-c) (minusp bound-c)
                        (numberp bound-d) (minusp bound-d)))
               (funcall if-not-zero))
              (t (funcall if-unknown))))
      (funcall if-unknown)))

(defparameter bilft-add
  (if (boundp 'bilft-add)
      bilft-add
      (make-bilft 0 1 1 0
                  0 0 0 1)))

(defun bilft-add (left right)
  (funcall bilft-add left right))

(defparameter bilft-multiply
  (if (boundp 'bilft-multiply)
      bilft-multiply
      (make-bilft 1 0 0 0
                  0 0 0 1)))

(defun bilft-multiply (left right)
  (funcall bilft-multiply left right))

(defmethod add2 ((left lft) (right lft))
  (compose-bilft-lft-x
   (compose-bilft-lft-y
    bilft-add
    right)
   left))

(defmethod mul2 ((left lft) (right lft))
  (compose-bilft-lft-x
   (compose-bilft-lft-y
    bilft-multiply
    right)
   left))

(defmethod add2 ((left rational) (right bilft))
  (compose-lft-bilft (lft-add-rat left) right))

(defmethod add2 ((left bilft) (right rational))
  (compose-lft-bilft (lft-add-rat right) left))

(defmethod mul2 ((left rational) (right bilft))
  (compose-lft-bilft (lft-multiply-by-rat left) right))

(defmethod mul2 ((left bilft) (right rational))
  (compose-lft-bilft (lft-multiply-by-rat right) left))

(defun bilft-truncate (bilft if-success if-failure)
  (check-type bilft bilft)
  (if (and (plusp (slot-value bilft 'e))
           (plusp (slot-value bilft 'f))
           (plusp (slot-value bilft 'g))
           (plusp (slot-value bilft 'h)))
      (let ((t1 (truncate (slot-value bilft 'a) (slot-value bilft 'e)))
            (t2 (truncate (slot-value bilft 'b) (slot-value bilft 'f)))
            (t3 (truncate (slot-value bilft 'c) (slot-value bilft 'g)))
            (t4 (truncate (slot-value bilft 'd) (slot-value bilft 'h))))
        (if (= t1 t2 t3 t4)
            (funcall if-success t1 (x- bilft t1))
            (funcall if-failure)))
      (funcall if-failure)))
