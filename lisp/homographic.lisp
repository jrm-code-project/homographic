;;; -*- Lisp -*-

(in-package "CL-USER")

;;; Basic arithmetic

(defun 2x2-matrix-multiply (a b
                            c d

                            e f
                            g h

                            receiver)
  (funcall receiver
           (+ (* a e) (* b g)) (+ (* a f) (* b h))
           (+ (* c e) (* d g)) (+ (* c f) (* d h))))

;;;
;;; Evaluators for homographic functions
;;;

(defun evaluate-homographic-function (a b
                                      c d n)
  "Return lim x->n (ax + b)/(cx + d)"
  (check-type a integer)
  (check-type b integer)
  (check-type c integer)
  (check-type d integer)
  (check-type n (or number (eql infinity)))
  (if (eq n 'infinity)
      (cond ((not (zerop c)) (/ a c))
            ((not (zerop a)) 'infinity)
            ((not (zerop d)) (/ b d))
            ((not (zerop b)) 'infinity)
            (t (error 'division-by-zero)))
      (let ((num (+ (* a n) b))
            (den (+ (* c n) d)))
        (cond ((not (zerop den)) (/ num den))
              ((not (zerop num)) 'infinity)
              ((= (* a d) (* b c)) (cond ((not (zerop c)) (/ a c))
                                         ((not (zerop d)) (/ b d))
                                         (t (error 'division-by-zero))))
              (t 'infinity)))))

(defun evaluate-bihomographic-function (a b c d
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
                              (evaluate-homographic-function (+ (* a n) b) (+ (* c n) d)
                                                             (+ (* e n) f) (+ (* g n) h) m)))
        ((eq n 'infinity) (evaluate-homographic-function (+ (* a m) c) (+ (* b m) d)
                                                         (+ (* e m) g) (+ (* f m) h) n))
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

;;;
;;; Printers for homographic equations
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

(defun print-homographic-equation (a b
                                   c d stream)
  (if (and (zerop c)
           (= d 1))
      (print-terms (list a b) '("x" nil) stream :suppress-parens t)
      (progn
        (print-terms (list a b) '("x" nil) stream)
        (format stream "/")
        (print-terms (list c d) '("x" nil) stream))))

(defun print-bihomographic-equation (a b c d
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
;;; Instances of homographic functions
;;;

(defclass homographic-function ()
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

(defmethod initialize-instance :after ((instance homographic-function) &rest initargs)
  (declare (ignore initargs))
  (sb-mop:set-funcallable-instance-function instance
                                            #'(lambda (x)
                                                (evaluate-homographic-function
                                                 (slot-value instance 'a)
                                                 (slot-value instance 'b)
                                                 (slot-value instance 'c)
                                                 (slot-value instance 'd)
                                                 x))))

(defmethod print-object ((instance homographic-function) stream)
  (print-unreadable-object (instance stream :type t :identity t)
    (print-homographic-equation
     (slot-value instance 'a)
     (slot-value instance 'b)
     (slot-value instance 'c)
     (slot-value instance 'd)
     stream)))

(defclass bihomographic-function ()
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

(defmethod initialize-instance :after ((instance bihomographic-function) &rest initargs)
  (declare (ignore initargs))
  (sb-mop:set-funcallable-instance-function instance
                                            #'(lambda (x y)
                                                (evaluate-bihomographic-function
                                                 (slot-value instance 'a)
                                                 (slot-value instance 'b)
                                                 (slot-value instance 'c)
                                                 (slot-value instance 'd)
                                                 (slot-value instance 'e)
                                                 (slot-value instance 'f)
                                                 (slot-value instance 'g)
                                                 (slot-value instance 'h)
                                                 x y))))

(defmethod print-object ((instance bihomographic-function) stream)
  (print-unreadable-object (instance stream :type t :identity t)
    (print-bihomographic-equation
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

(defun canonicalize-homographic-coefficients (a b
                                              c d receiver)
  (if (or (minusp c)
          (and (zerop c)
               (minusp d)))
      (canonicalize-homographic-coefficients (- a) (- b)
                                             (- c) (- d) receiver)
      (let ((gcd (gcd a b c d)))
        (if (> gcd 1)
            (canonicalize-homographic-coefficients (/ a gcd) (/ b gcd)
                                                   (/ c gcd) (/ d gcd) receiver)
            (funcall receiver
                     a b
                     c d)))))

(defun make-homographic-function (a b c d)
  (canonicalize-homographic-coefficients
   a b
   c d
   (lambda (a b
            c d)
     (make-instance 'homographic-function
                    :a a :b b
                    :c c :d d))))

(defun canonicalize-bihomographic-coefficients (a b c d e f g h receiver)
  (if (or (minusp e)
          (and (zerop e)
               (or (minusp f)
                   (and (zerop f)
                        (or (minusp g)
                            (and (zerop g)
                                 (minusp h)))))))
      (canonicalize-bihomographic-coefficients (- a) (- b) (- c) (- d)
                                               (- e) (- f) (- g) (- h) receiver)
      (let ((gcd (gcd a b c d e f g h)))
        (if (> gcd 1)
            (canonicalize-bihomographic-coefficients (/ a gcd) (/ b gcd) (/ c gcd) (/ d gcd)
                                                     (/ e gcd) (/ f gcd) (/ g gcd) (/ f gcd)
                                                     receiver)
            (funcall receiver
                     a b c d
                     e f g h)))))

(defun make-bihomographic-function (a b c d
                                    e f g h)
  (canonicalize-bihomographic-coefficients
   a b c d
   e f g h
   (lambda (a b c d
            e f g h)
     (make-instance 'bihomographic-function
                    :a a :b b :c c :d d
                    :e e :f f :g g :h h))))
     
;;;
;;; Functions that need access to the coefficients
;;;

(defun compose-two-homographic-functions (left right)
  (check-type left homographic-function)
  (check-type right homographic-function)
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
                         #'make-homographic-function)))
   

(defun pole-negative-p (homographic-function)
  (check-type homographic-function homographic-function)
  (or (and (plusp (slot-value homographic-function 'c))
           (plusp (slot-value homographic-function 'd)))
      (and (minusp (slot-value homographic-function 'c))
           (minusp (slot-value homographic-function 'd)))))

(defun poles-negative-p (bihomographic-function)
  (check-type bihomographic-function bihomographic-function)
  (or (and (plusp (slot-value bihomographic-function 'e))
           (plusp (slot-value bihomographic-function 'f))
           (plusp (slot-value bihomographic-function 'g))
           (plusp (slot-value bihomographic-function 'h)))
      (and (minusp (slot-value bihomographic-function 'e))
           (minusp (slot-value bihomographic-function 'f))
           (minusp (slot-value bihomographic-function 'g))
           (minusp (slot-value bihomographic-function 'h)))))

;;;
;;; Workhorse functions
;;;

(defun compose-homographic-functions (&rest homographic-functions)
  (fold-left #'compose-two-homographic-functions
             (make-homographic-function 1 0 0 1)
             homographic-functions))

(defun homographic-function-integer-value (homographic-function if-success if-failure)
  (check-type homographic-function homographic-function)
  (if (pole-negative-p homographic-function)
      (let ((lower-bound (funcall homographic-function 0))
            (upper-bound (funcall homographic-function 'infinity)))
        (if (and (numberp lower-bound)
                 (numberp upper-bound))
            (let ((lower-floor (floor lower-bound))
                  (upper-floor (floor upper-bound)))
              (if (= lower-floor upper-floor)
                  (funcall if-success lower-floor)
                  (funcall if-failure)))
            (funcall if-failure)))
      (funcall if-failure)))

(defun bihomographic-function-integer-value (bihomographic-function if-success if-failure)
  (check-type bihomographic-function bihomographic-function)
  (if (poles-negative-p bihomographic-function)
      (let ((bound-a (funcall bihomographic-function 0 0))
            (bound-b (funcall bihomographic-function 0 'infinity))
            (bound-c (funcall bihomographic-function 'infinity 0))
            (bound-d (funcall bihomographic-function 'infinity 'infinity)))
        (if (and (numberp bound-a)
                 (numberp bound-b)
                 (numberp bound-c)
                 (numberp bound-d))
            (let ((floor-a (floor bound-a))
                  (floor-b (floor bound-b))
                  (floor-c (floor bound-c))
                  (floor-d (floor bound-d)))
              (if (= floor-a floor-b floor-c floor-d)
                  (funcall if-success floor-a)
                  (funcall if-failure)))
            (funcall if-failure)))
      (funcall if-failure)))



