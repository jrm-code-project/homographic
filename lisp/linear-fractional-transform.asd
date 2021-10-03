;;; -*- Lisp -*-

(defsystem "linear-fractional-transform"
  :depends-on ("promise" "series" "stream" "utilities")
  :components ((:file "package")
               (:file "generics" :depends-on ("package"))
               (:file "lft" :depends-on ("generics" "package"))
               (:file "lft-stream" :depends-on ("generics" "lft" "package"))
               (:file "binary-expression" :depends-on ("generics" "lft" "lft-stream" "package"))
               (:file "cf-stream" :depends-on ("generics" "lft" "lft-stream" "package"))
               (:file "lft-arith" :depends-on ("binary-expression" "generics" "lft" "lft-stream" "package"))))
