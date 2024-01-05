;;; -*- Lisp -*-

(in-package "LINEAR-FRACTIONAL-TRANSFORMATION")

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
  (let ((delayed-left (slot-value binary-expression 'delayed-left)))
    (if (null delayed-left)
        (error "Attempt to refine-left on empty stream.")
        (let ((left (force delayed-left)))
          (make-instance 'binary-expression
                         :generation (+ (slot-value binary-expression 'generation) 1)
                         :bilft (compose-bilft-lft-x
                                 (slot-value binary-expression 'bilft)
                                 (if (empty-stream? left)
                                     (make-lft 1 1
                                               0 0)
                                     (stream-car left)))
                         :delayed-left (if (empty-stream? left)
                                           nil
                                           (stream-delayed-cdr left))
                         :delayed-right (slot-value binary-expression 'delayed-right))))))

(defun refine-right (binary-expression)
  (let ((delayed-right (slot-value binary-expression 'delayed-right)))
    (if (null delayed-right)
        (error "Attempt to refine-right on empty stream.")
        (let ((right (force delayed-right)))
          (make-instance 'binary-expression
                         :generation (+ (slot-value binary-expression 'generation) 1)
                         :bilft (compose-bilft-lft-y
                                 (slot-value binary-expression 'bilft)
                                 (if (empty-stream? right)
                                     (make-lft 1 1
                                               0 0)
                                     (stream-car right)))
                         :delayed-left (slot-value binary-expression 'delayed-left)
                         :delayed-right (if (empty-stream? right)
                                            nil
                                            (stream-delayed-cdr right)))))))

(defun refine-fairly (binary-expression)
  (if (zerop (mod (slot-value binary-expression 'generation) 2))
      (if (null (slot-value binary-expression 'delayed-left))
          (refine-right binary-expression)
          (refine-left binary-expression))
      (if (null (slot-value binary-expression 'delayed-right))
          (refine-left binary-expression)
          (refine-right binary-expression))))

(defun refine-disjoint (binary-expression)
  (if (bilft-disjoint? (slot-value binary-expression 'bilft))
      (refine-right binary-expression)
      (refine-left binary-expression)))

(defun binary-expression-minus-p (binary-expression)
  (bilft-minus-p (slot-value binary-expression 'bilft)
                 (constantly t)
                 (constantly nil)
                 (lambda ()
                   (binary-expression-minus-p (refine-disjoint binary-expression)))))

(defun binary-expression-zero-p (binary-expression)
  (bilft-zero-p (slot-value binary-expression 'bilft)
                (constantly t)
                (constantly nil)
                (lambda ()
                  (binary-expression-zero-p (refine-disjoint binary-expression)))))

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
     (binary-expression-truncate (refine-disjoint binary-expression)))))

(defun binary-expression-decompose-sign (binary-expression)
  (try-decompose-sign
   (slot-value binary-expression 'bilft)
   (lambda (sign factor)
     (values sign (make-instance 'binary-expression
                                  :generation (slot-value binary-expression 'generation)
                                  :bilft factor
                                  :delayed-left (slot-value binary-expression 'delayed-left)
                                  :delayed-right (slot-value binary-expression 'delayed-right))))
   (lambda ()
     (binary-expression-decompose-sign (refine-disjoint binary-expression)))))

(defun binary-expression-decompose-digit (binary-expression)
  (try-decompose-digit
   (slot-value binary-expression 'bilft)
   (lambda (digit factor)
     (values digit (make-instance 'binary-expression
                                  :generation (slot-value binary-expression 'generation)
                                  :bilft factor
                                  :delayed-left (slot-value binary-expression 'delayed-left)
                                  :delayed-right (slot-value binary-expression 'delayed-right))))
   (lambda ()
     (binary-expression-decompose-digit (refine-disjoint binary-expression)))))

(defun binary-expression->lft-digit-stream (binary-expression)
  (multiple-value-bind (digit remainder) (binary-expression-decompose-digit binary-expression)
    (cons-lft-stream digit 
                     (binary-expression->lft-digit-stream remainder))))

(defun binary-expression->lft-stream (binary-expression)
  (multiple-value-bind (sign remainder) (binary-expression-decompose-sign binary-expression)
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

(defun unfold-expression-tree (counter generate left)
  (funcall (funcall generate counter)
           left
           (delay-lft-stream (unfold-expression-tree (+ counter 1) generate left))))

(defun unfold-expression-tree-1 (root-bilft generate left)
  (funcall root-bilft
           left
           (delay-lft-stream (unfold-expression-tree 1 generate left))))

(defmethod add2 ((left lft-stream) (right lft-stream))
  (bilft-add left right))

(defmethod mul2 ((left lft-stream) (right lft-stream))
  (bilft-multiply left right))

(defun %atan-lft-stream (lft-stream)
  (unfold-expression-tree-1
   (make-bilft 1 1 -1 -1
               2 0  0  2)
   (lambda (n)
     (make-bilft (+ (* 2 n) 1) n 0 (+ n 1)
                 (+ n 1)       0 n (+ (* 2 n) 1)))
   (funcall (inverse-lft lft-sign-zero) lft-stream)))

(defun %exp-lft-stream (lft-stream)
  (unfold-expression-tree
   0 (lambda (n)
       (make-bilft (+ (* 2 n) 2) (+ (* 2 n) 1) (* 2 n) (+ (* 2 n) 1)
                   (+ (* 2 n) 1) (* 2 n) (+ (* 2 n) 1) (+ (* 2 n) 2)))
   (funcall (inverse-lft lft-sign-zero) lft-stream)))

(defun %log-lft-stream (lft-stream)
  (unfold-expression-tree-1
   (make-bilft 1 1 -1 -1
               0 1  1  0)
   (lambda (n)
     (make-bilft n (+ (* n 2) 1)       (+ n 1) 0
                 0       (+ n 1) (+ (* n 2) 1) n))
   lft-stream))

(defun %sqrt-lft-stream (lft-stream)
  ;; (unfold-expression-tree
  ;;  0 (constantly (make-bilft 1 2 1 0
  ;;                            0 1 2 1))
  ;;  lft-stream)
  (let ((stream nil))
    (setq stream
          (funcall (make-bilft 1 2 1 0
                               0 1 2 1)
                   lft-stream
                   (delay-lft-stream stream)))
    stream))

(defun %square-lft-stream (lft-stream)
  (bilft-multiply lft-stream lft-stream))

(defun %tan-lft-stream (lft-stream)
  (unfold-expression-tree-1
   (make-bilft 1 1 -1 -1
               2 0  0  2)
   (lambda (n)
     (make-bilft (+ (* 2 n) 1) (- (* 2 n) 1) (+ (* 2 n) 1) (+ (* 2 n) 3)
                 (+ (* 2 n) 3) (+ (* 2 n) 1) (- (* 2 n) 1) (+ (* 2 n) 1)))
   (funcall (inverse-lft lft-sign-zero) lft-stream)))

(defun %rat-tan (x)
  (check-type x (rational -1 1))
  (%tan-lft-stream
   (cf-stream->lft-stream
    (rational->cf-stream x))))
