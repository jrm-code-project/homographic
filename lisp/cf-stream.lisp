;;; -*- Lisp -*-

(in-package "LINEAR-FRACTIONAL-TRANSFORMATION")

(eval-when (:compile-toplevel :load-toplevel :execute)

(defclass cf-stream (stream::stream)
  ())

(defmacro cons-cf-stream (stream-car stream-cdr)
  `(make-instance 'cf-stream
                  :car ,stream-car
                  :delayed-cdr (delay ,stream-cdr)))
)

(defun cf-stream->lft-stream (cf-stream)
  (cond ((empty-stream? cf-stream) the-empty-stream)
        ((or (minusp (stream-car cf-stream))
             (and (zerop (stream-car cf-stream))
                  (minusp (stream-car (stream-cdr cf-stream)))))
         (cons-lft-stream (make-lft -1 0
                                     0 1)
                          (cf-stream->lft-stream (negate cf-stream))))
        (t (cons-lft-stream (make-lft (stream-car cf-stream) 1
                                      1                      0)
                            (cf-stream->lft-stream (stream-cdr cf-stream))))))

(defmethod print-object ((object cf-stream) cl:stream)
  (print-unreadable-object (object cl:stream :type t)
    (print-lft-stream-as-decimal (cf-stream->lft-stream object) cl:stream)))

(defun lft-stream->cf-stream (lft-stream)
  (if (empty-stream? lft-stream)
      the-empty-stream
      (multiple-value-bind (integer-part fractional-part) (lft-stream-truncate lft-stream)
        (cons-cf-stream integer-part
                        (multiple-value-bind (zerop refined) (lft-stream-zero-p fractional-part)
                          (if zerop
                              the-empty-stream
                              (lft-stream->cf-stream (reciprocal refined))))))))

(defun canonicalize-cf-stream (cf-stream)
  (if (and (= 1 (stream-car (stream-cdr cf-stream)))
           (empty-stream? (stream-cdr (stream-cdr cf-stream))))
      (cons-stream (+ (stream-car cf-stream) 1) the-empty-stream)
      (cons-stream (stream-car cf-stream) (canonicalize-cf-stream (stream-cdr cf-stream)))))

(defun cf-stream-truncate-delayed (delayed-cf-stream n-elements)
  (if (zerop n-elements)
      the-empty-stream
      (let ((cf-stream (force delayed-cf-stream)))
        (cons-cf-stream (stream-car cf-stream)
                        (cf-stream-truncate-delayed (stream-delayed-cdr cf-stream) (- n-elements 1))))))

(defun cf-stream-truncate (cf-stream n-elements)
  (if (zerop n-elements)
      the-empty-stream
      (cons-cf-stream (stream-car cf-stream)
                      (cf-stream-truncate-delayed (stream-delayed-cdr cf-stream) (- n-elements 1)))))

(defun coefficients->cf-stream (coefficients)
  (if (null coefficients)
      the-empty-stream
      (cons-cf-stream (car coefficients) (coefficients->cf-stream (cdr coefficients)))))

(defun rational->cf-stream (rational)
  (multiple-value-bind (integer-part remainder) (truncate rational)
    (cons-cf-stream integer-part (if (zerop remainder)
                                     the-empty-stream
                                     (rational->cf-stream (reciprocal remainder))))))

(defgeneric ->cf-stream (object)
  (:method ((object lft-stream))
    (lft-stream->cf-stream object))
  (:method ((object ratio))
    (rational->cf-stream object))
  (:method ((object single-float))
    (rational->cf-stream (rational object)))
  (:method ((object double-float))
    (rational->cf-stream (rational object))))

(defun cf-stream->rational (cf-stream)
  (if (empty-stream? cf-stream)
      'infinity
      (+ (stream-car cf-stream) (reciprocal (cf-stream->rational (stream-cdr cf-stream))))))

(defmethod reciprocal ((object cf-stream))
  (if (zerop (stream-car object))
      (stream-cdr object)
      (cons-cf-stream 0 object)))

(defmethod negate ((object cf-stream))
  (if (empty-stream? object)
      the-empty-stream
      (cons-cf-stream (- (stream-car object)) (negate (stream-cdr object)))))

(defun cf-stream-prefix (left right)
  (if (and (typep left 'cf-stream)
           (typep right 'cf-stream)
           (= (stream-car left) (stream-car right)))
      (cons-cf-stream (stream-car left) (cf-stream-prefix (stream-cdr left) (stream-cdr right)))
      the-empty-stream))

(defun cf-stream-suffix (left right)
  (if (and (typep left 'cf-stream)
           (typep right 'cf-stream)
           (= (stream-car left) (stream-car right)))
      (cf-stream-suffix (stream-cdr left) (stream-cdr right))
      right))

(defun cf-stream-append2-delayed (left delayed-right)
  (if (empty-stream? left)
      (force delayed-right)
      (cons-cf-stream (stream-car left)
                      (cf-stream-append2-delayed (stream-cdr left) delayed-right))))

(defun cf-stream-flatten-append (stream-of-streams)
  (stream-fold-right-delayed #'cf-stream-append2-delayed '() stream-of-streams))

;; (defun limit-stream->cf-stream (limit-stream)
;;   (let* ((cfs (stream-map #'->cf-stream limit-stream))
;;          (prefixes (cons-stream the-empty-stream (stream-map #'cf-stream-prefix cfs (stream-cdr cfs)))))
;;     (cf-stream-flatten-append (stream-map #'cf-stream-suffix prefixes (stream-cdr prefixes)))))
