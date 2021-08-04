;;; -*- Lisp -*-

(in-package "LINEAR-FRACTIONAL-TRANSFORM")

(defclass binary-expression ()
  ((generation :initarg :generation
               :initform 0)
   (bilft :initarg :bilft
          :initform (error "Required initarg :bilft was omitted."))
   (delayed-left :initarg :delayed-left
                 :initform (error "Required initarg :delayed-left was omitted."))
   (delayed-right :initarg :delayed-right
                  :initform (error "Required initarg :delayed-right was omitted."))))

(defmethod funcall-lft (lft (arg binary-expression))
  (make-instance 'binary-expression
                 :generation (slot-value arg 'generation)
                 :bilft (funcall lft (slot-value arg 'bilft))
                 :delayed-left (slot-value arg 'delayed-left)
                 :delayed-right (slot-value arg 'delayed-right)))

(defun refine-left (binary-expression)
  (let ((left (force (slot-value binary-expression 'delayed-left))))
    (make-instance 'binary-expression
                   :generation (+ (slot-value binary-expression 'generation) 1)
                   :bilft (compose-bilft-lft-x
                           (slot-value binary-expression 'bilft)
                           (if (empty-stream? left)
                               (make-lft 1 0
                                         0 0)
                               (stream-car left)))
                   :delayed-left (if (empty-stream? left)
                                     (delay the-empty-stream)
                                     (stream-delayed-cdr left))
                   :delayed-right (slot-value binary-expression 'delayed-right))))

(defun refine-right (binary-expression)
  (let ((right (force (slot-value binary-expression 'delayed-right))))
    (make-instance 'binary-expression
                   :generation (+ (slot-value binary-expression 'generation) 1)
                   :bilft (compose-bilft-lft-y
                           (slot-value binary-expression 'bilft)
                           (if (empty-stream? right)
                               (make-lft 1 0
                                         0 0)
                               (stream-car right)))
                   :delayed-left (slot-value binary-expression 'delayed-left)
                   :delayed-right (if (empty-stream? right)
                                      (delay the-empty-stream)
                                      (stream-delayed-cdr right)))))

(defun refine-fairly (binary-expression)
  (if (zerop (mod (slot-value binary-expression 'generation) 2))
      (refine-left binary-expression)
      (refine-right binary-expression)))

(defun refine-widest (binary-expression)
  (let* ((bilft (slot-value binary-expression 'bilft))
         (atzero (funcall bilft 0 0))
         (atxinf (funcall bilft 'infinity 0))
         (atyinf (funcall bilft 0 'infinity)))
    (cond ;((and (not (numberp atxinf)) (numberp atyinf)) (refine-right binary-expression))
          ;((and (numberp atxinf) (not (numberp atyinf))) (refine-left binary-expression))
          ((or (not (numberp atzero))
               (not (numberp atxinf))
               (not (numberp atyinf)))
           (refine-fairly binary-expression))
          (t
           (let ((dx (abs (- atxinf atzero)))
                 (dy (abs (- atyinf atzero))))
             (cond ((= dy dx) (refine-fairly binary-expression))
                   ((> dx dy) (refine-left binary-expression))
                   (t (refine-right binary-expression))))))))

(defun binary-expression-minus-p (binary-expression)
  (bilft-minus-p (slot-value binary-expression 'bilft)
                 (constantly t)
                 (constantly nil)
                 (lambda ()
                   (binary-expression-minus-p (refine-widest binary-expression)))))

(defun binary-expression-zero-p (binary-expression)
  (bilft-zero-p (slot-value binary-expression 'bilft)
                (constantly t)
                (constantly nil)
                (lambda ()
                  (binary-expression-zero-p (refine-widest binary-expression)))))

(defun binary-expression-truncate (binary-expression)
  (bilft-truncate
   (slot-value binary-expression 'bilft)
   (lambda (integer-value fractional-part)
     (values integer-value
             (make-instance 'binary-expression
                            :generation (slot-value binary-expression 'generation)
                            :bilft fractional-part
                            :delayed-left (slot-value binary-expression 'delayed-left)
                            :delayed-right (slot-value binary-expression 'delayed-right))))
   (lambda ()
     (binary-expression-truncate (refine-widest binary-expression)))))

(defun bilft-compose-refine? (lft bilft if-success if-failure)
  (let ((attempt (funcall lft bilft)))
    (if (refine? attempt)
        (funcall if-success attempt)
        (funcall if-failure))))

(defun bilft-try-factor-sign (bilft if-success if-failure)
  (bilft-compose-refine?
   inverse-lft-sign-positive bilft
   (lambda (new-bilft)
     (funcall if-success lft-sign-positive new-bilft))
   (lambda ()
     (bilft-compose-refine?
      inverse-lft-sign-negative bilft
      (lambda (new-bilft)
        (funcall if-success lft-sign-negative new-bilft))
      (lambda ()
        (bilft-compose-refine?
         inverse-lft-sign-zero bilft
         (lambda (new-bilft)
           (funcall if-success lft-sign-zero new-bilft))
         (lambda ()
           (bilft-compose-refine?
            inverse-lft-sign-infinity bilft
            (lambda (new-bilft)
              (funcall if-success lft-sign-infinity new-bilft))
            if-failure))))))))

(defun bilft-try-factor-digit (bilft if-success if-failure)
  (bilft-compose-refine?
   inverse-lft-digit-one bilft
   (lambda (new-bilft)
     (funcall if-success lft-digit-one new-bilft))
   (lambda ()
     (bilft-compose-refine?
      inverse-lft-digit-minus-one bilft
      (lambda (new-bilft)
        (funcall if-success lft-digit-minus-one new-bilft))
      (lambda ()
        (bilft-compose-refine?
         inverse-lft-digit-zero bilft
         (lambda (new-bilft)
           (funcall if-success lft-digit-zero new-bilft))
         if-failure))))))

(defun binary-expression-factor-sign (binary-expression)
  (bilft-try-factor-sign
   (slot-value binary-expression 'bilft)
   (lambda (sign factor)
     (values sign (make-instance 'binary-expression
                                  :generation (slot-value binary-expression 'generation)
                                  :bilft factor
                                  :delayed-left (slot-value binary-expression 'delayed-left)
                                  :delayed-right (slot-value binary-expression 'delayed-right))))
   (lambda ()
     (binary-expression-factor-sign (refine-widest binary-expression)))))

(defun binary-expression-factor-digit (binary-expression)
  (bilft-try-factor-digit
   (slot-value binary-expression 'bilft)
   (lambda (digit factor)
     (values digit (make-instance 'binary-expression
                                  :generation (slot-value binary-expression 'generation)
                                  :bilft factor
                                  :delayed-left (slot-value binary-expression 'delayed-left)
                                  :delayed-right (slot-value binary-expression 'delayed-right))))
   (lambda ()
     (binary-expression-factor-digit (refine-widest binary-expression)))))

(defun binary-expression->lft-digit-stream (binary-expression)
  (multiple-value-bind (digit remainder) (binary-expression-factor-digit binary-expression)
    (cons-lft-stream digit 
                     (binary-expression->lft-digit-stream remainder))))

(defun binary-expression->lft-stream (binary-expression)
  (multiple-value-bind (sign remainder) (binary-expression-factor-sign binary-expression)
    (cons-lft-stream sign (binary-expression->lft-digit-stream remainder))))

(defmethod print-object ((obj binary-expression) cl:stream)
  (print-unreadable-object (obj cl:stream :type t)
    (print-lft-stream-as-decimal
     (binary-expression->lft-stream obj)
     cl:stream)))

(defmethod funcall-bilft (bilft (x lft-stream) (y lft-stream))
  (binary-expression->lft-stream
   (make-instance 'binary-expression
                  :generation 0
                  :bilft (compose-bilft-lft-x
                          (compose-bilft-lft-y
                           bilft
                           (stream-car y))
                          (stream-car x))
                  :delayed-left (stream-delayed-cdr x)
                  :delayed-right (stream-delayed-cdr y))))

(defmethod funcall-bilft (bilft (x lft-stream) (y promise))
  (binary-expression->lft-stream
   (make-instance 'binary-expression
                  :generation 0
                  :bilft (compose-bilft-lft-x
                           bilft
                           (stream-car x))
                  :delayed-left (stream-delayed-cdr x)
                  :delayed-right y)))

(defun unfold-expression-tree (counter generate left)
  (binary-expression->lft-stream
   (make-instance 'binary-expression
                  :generation 0
                  :bilft (compose-bilft-lft-x (funcall generate counter) (stream-car left))
                  :delayed-left (stream-delayed-cdr left)
                  :delayed-right (delay (unfold-expression-tree (+ counter 1) generate left)))))

(defmethod add2 ((left lft-stream) (right lft-stream))
  (bilft-add left right))

(defmethod mul2 ((left lft-stream) (right lft-stream))
  (bilft-multiply left right))

(defun %atan-lft-stream (lft-stream)
  (let ((lft-stream* (funcall inverse-lft-sign-zero lft-stream)))
    (binary-expression->lft-stream
     (make-instance 'binary-expression
                    :generation 0
                    :bilft (compose-bilft-lft-x
                            (make-bilft 1 1 -1 -1
                                        2 0  0  2)
                            (stream-car lft-stream*))
                    :delayed-left (stream-delayed-cdr lft-stream*)
                    :delayed-right (delay (unfold-expression-tree
                                           1 (lambda (n)
                                               (make-bilft (+ (* 2 n) 1) n 0 (+ n 1)
                                                           (+ n 1)       0 n (+ (* 2 n) 1)))
                                           lft-stream*))))))

(defun %exp-lft-stream (lft-stream)
  (unfold-expression-tree
   0 (lambda (n)
       (make-bilft (+ (* 2 n) 2) (+ (* 2 n) 1) (* 2 n) (+ (* 2 n) 1)
                   (+ (* 2 n) 1) (* 2 n) (+ (* 2 n) 1) (+ (* 2 n) 2)))
   (inverse-lft-sign-zero lft-stream)))

(defun %log-lft-stream (lft-stream)
  (binary-expression->lft-stream
   (make-instance 'binary-expression
                  :generation 0
                  :bilft (compose-bilft-lft-x
                          (make-bilft 1 1 -1 -1
                                      0 1  1  0)
                          (stream-car lft-stream))
                  :delayed-left (stream-delayed-cdr lft-stream)
                  :delayed-right (delay (unfold-expression-tree
                                         1 (lambda (n)
                                             (make-bilft n (+ (* n 2) 1)       (+ n 1) 0
                                                         0       (+ n 1) (+ (* n 2) 1) n))
                                         lft-stream)))))

(defun %sqrt-lft-stream (lft-stream)
  (unfold-expression-tree
   0 (constantly (make-bilft 1 2 1 0
                             0 1 2 1))
   lft-stream))

(defun %square-lft-stream (lft-stream)
  (bilft-multiply lft-stream lft-stream))

(defun %tan-lft-stream (lft-stream)
  (let ((lft-stream* (inverse-lft-sign-zero lft-stream)))
    (binary-expression->lft-stream
     (make-instance 'binary-expression
                    :generation 0
                    :bilft (compose-bilft-lft-x
                            (make-bilft 1 1 -1 -1
                                        2 0  0  2)
                            (stream-car lft-stream*))
                    :delayed-left (stream-delayed-cdr lft-stream*)
                    :delayed-right (delay
                                    (unfold-expression-tree
                                     1 (lambda (n)
                                         (make-bilft (+ (* 2 n) 1) (- (* 2 n) 1) (+ (* 2 n) 1) (+ (* 2 n) 3)
                                                     (+ (* 2 n) 3) (+ (* 2 n) 1) (- (* 2 n) 1) (+ (* 2 n) 1)))
                                     lft-stream*))))))
