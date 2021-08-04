;;; -*- Lisp -*-

(in-package "LINEAR-FRACTIONAL-TRANSFORM")

;;; Arithmetic utilities

(defun 2x2-matrix-multiply (a b
                            c d

                            e f
                            g h

                            receiver)
  "(funcall receiver i j k l) where
   [i j]   [a b]   [e f]   
   [k l] = [c d] x [g h]
"
  (funcall receiver
           (+ (* a e) (* b g)) (+ (* a f) (* b h))
           (+ (* c e) (* d g)) (+ (* c f) (* d h))))

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
  (funcall receiver
           d    (- b)
           (- c)   a))

(defun 2mx2t-multiply (a b
                       c d

                       e f g h
                       i j k l

                       receiver)
  ;; Multiply a 2x2 matrix by a 2x2x2 tensor yielding a 2x2x2 tensor
  (funcall receiver
           (+ (* a e) (* b i)) (+ (* a f) (* b j)) (+ (* a g) (* b k)) (+ (* a h) (* b l))
           (+ (* c e) (* d i)) (+ (* c f) (* d j)) (+ (* c g) (* d k)) (+ (* c h) (* d l))))

(defun 2tx2m-multiply-x (a b c d
                         e f g h

                         i j
                         k l

                         receiver)
  ;; Multiply a 2x2x2 tensor by a 2x2x2 matrix in the x direction yielding a 2x2x2 tensor
  (funcall receiver
           (+ (* a i) (* c k)) (+ (* b i) (* d k)) (+ (* a j) (* c l)) (+ (* b j) (* d l))
           (+ (* e i) (* g k)) (+ (* f i) (* h k)) (+ (* e j) (* g l)) (+ (* f j) (* h l))))

(defun 2tx2m-multiply-y (a b c d
                         e f g h

                         i j
                         k l

                         receiver)
  ;; Multiply a 2x2x2 tensor by a 2x2x2 matrix in the y direction yielding a 2x2x2 tensor
  (funcall receiver
           (+ (* a i) (* b k)) (+ (* a j) (* b l)) (+ (* c i) (* d k)) (+ (* c j) (* d l))
           (+ (* e i) (* f k)) (+ (* e j) (* f l)) (+ (* g i) (* h k)) (+ (* g j) (* h l))))

(defun sign (a b)
  (cond ((< a 0) (if (<= b 0)
                     -1
                     0))
        ((= a 0) (cond ((< b 0) -1)
                       ((= b 0) 0)
                       (t 1)))
        (t (if (< b 0)
               0
               1))))

(defun 2x2-matrix-refine? (a b
                           c d)
  (let ((a (sign a c))
        (b (sign b d)))
    (and (= a b)
         (not (= b 0)))))

(defun 2x2x2-tensor-refine? (a b c d
                             e f g h)
  (let ((a (sign a e))
        (b (sign b f))
        (c (sign c g))
        (d (sign d h)))
    (and (= a b)
         (= b c)
         (= c d)
         (not (= d 0)))))

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
                                     ((not (zerop result-numerator)) 'infinity)
                                     (t (error 'division-by-zero))))
          (t 'infinity))))

(defun new-evaluate-bilft (a b c d
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
                                 (t (error 'division-by-zero))))
          (t 'infinity))))

(defun old-evaluate-bilft (a b c d
                           e f g h m n)
  "Return lim x->m, y->n (axy + bx + cy + d)/(exy + fx + gy + h)"
  (check-type a integer)
  (check-type b integer)
  (check-type c integer)
  (check-type d integer)
  (check-type e integer)
  (check-type f integer)
  (check-type g integer)
  (check-type h integer)
  (check-type m (or number (eql infinity)))
  (check-type n (or number (eql infinity)))
  (cond ((eq m 'infinity) (if (eq n 'infinity)
                              (cond ((not (zerop e)) (/ a e))
                                    ((not (zerop a)) 'infinity)
                                    ((not (zerop (+ f g))) (/ (+ b c) (+ f g)))
                                    ((not (zerop (+ b c))) 'infinity)
                                    ((not (zerop h)) (/ d h))
                                    (t (error 'division-by-zero)))
                              (evaluate-lft (+ (* a n) b) (+ (* c n) d)
                                            (+ (* e n) f) (+ (* g n) h) 1 0)))
        ((eq n 'infinity) (evaluate-lft (+ (* a m) c) (+ (* b m) d)
                                        (+ (* e m) g) (+ (* f m) h) 1 0))
        (t (let ((num (+ (* a m n) (* b m) (* c n) d))
                 (den (+ (* e m n) (* f m) (* g n) h)))
             (cond ((not (zerop den)) (/ num den))
                   ((not (zerop num)) 'infinity)
                   ((= (* a f g h)
                       (* b e g h)
                       (* c e f h)
                       (* d e f g)) (cond ((not (zerop h)) (/ d h))
                                          ((not (zerop g)) (/ c g))
                                          ((not (zerop f)) (/ b f))
                                          ((not (zerop e)) (/ a e))
                       (t (error 'division-by-zero))))
                   (t 'infinity))))))

(defun evaluate-bilft (a b c d
                       e f g h m-num m-den n-num n-den)
  (let ((answer1 (new-evaluate-bilft a b c d
                                     e f g h m-num m-den n-num n-den))
        (answer2 (old-evaluate-bilft a b c d
                                     e f g h (if (zerop m-den)
                                                 'infinity
                                                 (/ m-num m-den))
                                     (if (zerop n-den)
                                         'infinity
                                         (/ n-num n-den)))))
    (if (or (and (numberp answer1)
                 (numberp answer2)
                 (= answer1 answer2))
            (and (eq answer1 'infinity)
                 (eq answer2 'infinity)))
        answer1
        (progn (format t "Answers differ for ~d ~d ~d ~d, ~d ~d ~d ~d, ~d ~d, ~d ~d, old: ~s, new: ~s.~%"
                       a b c d e f g h m-num m-den n-num n-den answer2 answer1)
               answer2))))

;;;
;;; Printers for lft equations
;;;

(defun print-term (coefficient variables stream &key initial-term)
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

(defun print-terms (coefficients variables stream &key suppress-parens)
  "Format a list of terms given a list of coefficients and variables."
  (labels ((l1 (coefficients variables)
             (cond ((null coefficients) (format stream "0"))
                   ((zerop (car coefficients)) (l1 (cdr coefficients) (cdr variables)))
                   ((or suppress-parens
                        (every #'zerop (cdr coefficients))) (l2 coefficients variables t))
                   (t (format stream "(")
                      (l2 coefficients variables t)
                      (format stream ")"))))

           (l2 (coefficients variables initial-term)
             (unless (null coefficients)
               (print-term (car coefficients) (car variables) stream :initial-term initial-term)
               (l2 (cdr coefficients) (cdr variables) nil))))

    (l1 coefficients variables)))

(defun print-lft-equation (a b
                           c d stream)
  (cond ((and (zerop c)
              (= d 1))
         (print-terms (list a b) '("x" nil) stream :suppress-parens t))
        ((and (= c 1)
              (zerop d))
         (if (zerop a)
             (format stream "~d/x" b)
             (format stream "~d ~:[+~;-~] ~d/x" a (minusp b) (abs b))))
        (t
         (print-terms (list a b) '("x" nil) stream)
         (format stream "/")
         (print-terms (list c d) '("x" nil) stream))))

(defun print-bilft-equation (a b c d
                             e f g h stream)
  (if (and (zerop e)
           (zerop f)
           (zerop g)
           (= h 1))
      (print-terms (list a b c d) '("xy" "x" "y" nil) stream :suppress-parens t)
      (progn
        (print-terms (list a b c d) '("xy" "x" "y" nil) stream)
        (format stream "/")
        (print-terms (list e f g h) '("xy" "x" "y" nil) stream))))

;;;
;;; Instances of lft functions
;;;

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

(defmethod print-object ((instance lft) stream)
  (print-unreadable-object (instance stream :type t)
    (print-lft-equation
     (slot-value instance 'a)
     (slot-value instance 'b)
     (slot-value instance 'c)
     (slot-value instance 'd)
     stream)))

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

(defmethod print-object ((instance bilft) stream)
  (print-unreadable-object (instance stream :type t)
    (print-bilft-equation
     (slot-value instance 'a)
     (slot-value instance 'b)
     (slot-value instance 'c)
     (slot-value instance 'd)
     (slot-value instance 'e)
     (slot-value instance 'f)
     (slot-value instance 'g)
     (slot-value instance 'h)
     stream)))

;;;
;;; Constructors
;;;

(defun canonicalize-lft-coefficients (a b
                                      c d receiver)
  (if (or (minusp c)
          (and (zerop c)
               (minusp d)))
      (canonicalize-lft-coefficients (- a) (- b)
                                     (- c) (- d) receiver)
      (let ((gcd (gcd a b c d)))
        (if (> gcd 1)
            (canonicalize-lft-coefficients (/ a gcd) (/ b gcd)
                                           (/ c gcd) (/ d gcd) receiver)
            (funcall receiver
                     a b
                     c d)))))

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
  (cond
    ((typep a 'float) (make-lft (rational a) b c d))
    ((typep b 'float) (make-lft a (rational b) c d))
    ((typep c 'float) (make-lft a b (rational c) d))
    ((typep d 'float) (make-lft a b c (rational d)))
    (t (etypecase a
         (integer
          (etypecase b
            (integer
             (etypecase c
               (integer
                (etypecase d
                  (integer (%make-lft a b
                                      c d))
                  (rational (make-lft (* a (denominator d)) (* b (denominator d))
                                      (* c (denominator d)) (numerator d)))))
               (rational (make-lft (* a (denominator c)) (* b (denominator c))
                                   (numerator c)         (* d (denominator c))))))
            (rational (make-lft (* a (denominator b)) (numerator b)
                                (* c (denominator b)) (* d (denominator b))))))
         (rational (make-lft (numerator a)         (* b (denominator a))
                             (* c (denominator a)) (* d (denominator a))))))))
(defconstant lft-identity
  (if (boundp 'lft-identity)
      lft-identity
      (make-lft 1 0
                0 1)))

(defconstant lft-negate
  (if (boundp 'lft-negate)
      lft-negate
      (make-lft -1 0
                0  1)))

(defmethod negate ((lft lft))
  (funcall lft-negate lft))

(defconstant lft-reciprocal
  (if (boundp 'lft-reciprocal)
      lft-reciprocal
      (make-lft 0 1
                1 0)))

(defmethod reciprocal ((lft lft))
  (funcall lft-reciprocal lft))

(defun canonicalize-bilft-coefficients (a b c d e f g h receiver)
  (if (or (minusp e)
          (and (zerop e)
               (or (minusp f)
                   (and (zerop f)
                        (or (minusp g)
                            (and (zerop g)
                                 (minusp h)))))))
      (canonicalize-bilft-coefficients (- a) (- b) (- c) (- d)
                                       (- e) (- f) (- g) (- h) receiver)
      (let ((gcd (gcd a b c d e f g h)))
        (if (> gcd 1)
            (canonicalize-bilft-coefficients (/ a gcd) (/ b gcd) (/ c gcd) (/ d gcd)
                                             (/ e gcd) (/ f gcd) (/ g gcd) (/ h gcd)
                                             receiver)
            (funcall receiver
                     a b c d
                     e f g h)))))

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
  (cond
    ((typep a 'float) (make-bilft (rational a) b c d e f g h))
    ((typep b 'float) (make-bilft a (rational b) c d e f g h))
    ((typep c 'float) (make-bilft a b (rational c) d e f g h))
    ((typep d 'float) (make-bilft a b c (rational d) e f g h))
    ((typep e 'float) (make-bilft a b c d (rational e) f g h))
    ((typep f 'float) (make-bilft a b c d e (rational f) g h))
    ((typep g 'float) (make-bilft a b c d e f (rational g) h))
    ((typep h 'float) (make-bilft a b c d e f g (rational h)))
    (t (etypecase a
         (integer
          (etypecase b
            (integer
             (etypecase c
               (integer
                (etypecase d
                  (integer
                   (etypecase e
                     (integer
                      (etypecase f
                        (integer
                         (etypecase g
                           (integer
                            (etypecase h
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
           (* e (denominator a)) (* f (denominator a)) (* g (denominator a)) (* h (denominator a))))))))

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

(defgeneric refine? (object)
  (:method ((object lft))
    (2x2-matrix-refine? (slot-value object 'a) (slot-value object 'b)
                        (slot-value object 'c) (slot-value object 'd)))
  (:method ((object bilft))
    (2x2x2-tensor-refine?
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

;;;
;;; Workhorse functions
;;;

(defun compose-lfts (lft &rest lfts)
  (fold-left #'compose-lft-lft
             lft
             lfts))

(defun lft-zero-p (lft if-zero if-not-zero if-unknown)
  (if (pole-negative-p lft)
      (let ((bound-a (funcall lft 0))
            (bound-b (funcall lft 'infinity)))
        (cond ((eq bound-a 'infinity) (funcall if-unknown))
              ((eq bound-b 'infinity) (funcall if-unknown))
              ((and (zerop bound-a)
                    (zerop bound-b))
               (funcall if-zero))
              ((or (and (plusp bound-a)
                        (plusp bound-b))
                   (and (minusp bound-a)
                        (minusp bound-b)))
               (funcall if-not-zero))
              (t (funcall if-unknown))))
      (funcall if-unknown)))

(defun lft-minus-p (lft if-minus if-plus if-unknown)
  (check-type lft lft)
  (if (pole-negative-p lft)
      (let ((bound-a (funcall lft 0))
            (bound-b (funcall lft 'infinity)))
        (cond ((eq bound-a 'infinity) (funcall if-unknown))
              ((eq bound-b 'infinity) (funcall if-unknown))
              ((and (or (minusp bound-a) (minusp bound-b))
                    (not (plusp bound-a))
                    (not (plusp bound-b)))
               (funcall if-minus))
              ((and (or (plusp bound-a) (plusp bound-b))
                    (not (minusp bound-a))
                    (not (minusp bound-b)))
               (funcall if-plus))
              (t (funcall if-unknown))))
      (funcall if-unknown)))

(defun lft-less-than-rat (lft rat if-less-than if-not-less-than if-unknown)
  (check-type lft lft)
  (if (pole-negative-p lft)
      (let ((bound-a (funcall lft 0))
            (bound-b (funcall lft 'infinity)))
        (cond ((eq bound-a 'infinity) (funcall if-unknown))
              ((eq bound-b 'infinity) (funcall if-unknown))
              ((and (< bound-a rat)
                    (< bound-b rat))
               (funcall if-less-than))
              ((and (> bound-a rat)
                    (> bound-b rat))
               (funcall if-not-less-than))
              (t (funcall if-unknown))))
      (funcall if-unknown)))

(defun lft-greater-than-rat (lft rat if-greater-than if-not-greater-than if-unknown)
  (check-type lft lft)
  (if (pole-negative-p lft)
      (let ((bound-a (funcall lft 0))
            (bound-b (funcall lft 'infinity)))
        (cond ((eq bound-a 'infinity) (funcall if-unknown))
              ((eq bound-b 'infinity) (funcall if-unknown))
              ((and (> bound-a rat)
                    (> bound-b rat))
               (funcall if-greater-than))
              ((and (< bound-a rat)
                    (< bound-b rat))
               (funcall if-not-greater-than))
              (t (funcall if-unknown))))
      (funcall if-unknown)))

(defun lft-add-rat (rat)
  (make-lft 1 rat
            0 1))

(defun lft-multiply-by-rat (rat)
  (make-lft rat 0
            0   1))

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
  (if (pole-negative-p lft)
      (let ((lower-bound (funcall lft 0))
            (upper-bound (funcall lft 'infinity)))
        (if (and (numberp lower-bound)
                 (numberp upper-bound))
            (let ((lower-floor (truncate lower-bound))
                  (upper-floor (truncate upper-bound)))
              (if (= lower-floor upper-floor)
                  (funcall if-success lower-floor (x- lft lower-floor))
                  (funcall if-failure)))
            (funcall if-failure)))
      (funcall if-failure)))

(defun lft->single (lft if-success if-failure)
  (check-type lft lft)
  (if (pole-negative-p lft)
      (let ((lower-bound (funcall lft 0))
            (upper-bound (funcall lft 'infinity)))
        (if (and (numberp lower-bound)
                 (numberp upper-bound))
            (let ((lower-single (coerce lower-bound 'single-float))
                  (upper-single (coerce upper-bound 'single-float)))
              (if (= lower-single upper-single)
                  (funcall if-success lower-single)
                  (funcall if-failure)))
            (funcall if-failure)))
      (funcall if-failure)))

(defun lft->double (lft if-success if-failure)
  (check-type lft lft)
  (if (pole-negative-p lft)
      (let ((lower-bound (funcall lft 0))
            (upper-bound (funcall lft 'infinity)))
        (if (and (numberp lower-bound)
                 (numberp upper-bound))
            (let ((lower-single (coerce lower-bound 'double-float))
                  (upper-single (coerce upper-bound 'double-float)))
              (if (= lower-single upper-single)
                  (funcall if-success lower-single)
                  (funcall if-failure)))
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

(defconstant bilft-add
  (if (boundp 'bilft-add)
      bilft-add
      (make-bilft 0 1 1 0
                  0 0 0 1)))

(defun bilft-add (left right)
  (funcall bilft-add left right))

(defconstant bilft-multiply
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
  (if (poles-negative-p bilft)
      (let ((bound-a (funcall bilft 0 0))
            (bound-b (funcall bilft 0 'infinity))
            (bound-c (funcall bilft 'infinity 0))
            (bound-d (funcall bilft 'infinity 'infinity)))
        (if (and (numberp bound-a)
                 (numberp bound-b)
                 (numberp bound-c)
                 (numberp bound-d))
            (let ((floor-a (truncate bound-a))
                  (floor-b (truncate bound-b))
                  (floor-c (truncate bound-c))
                  (floor-d (truncate bound-d)))
              (if (= floor-a floor-b floor-c floor-d)
                  (funcall if-success floor-a (x- bilft floor-a))
                  (funcall if-failure)))
            (funcall if-failure)))
      (funcall if-failure)))
