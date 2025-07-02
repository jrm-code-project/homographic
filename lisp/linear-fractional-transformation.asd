;;; -*- Lisp -*-

(defsystem "linear-fractional-transformation"
  :depends-on ("fold" "generic-arithmetic" "named-let" "promise" "series" "stream" "utilities")
  :components ((:file "package")
               (:file "lft" :depends-on ("package"))
               (:file "lft-stream" :depends-on ("lft" "package"))
               (:file "binary-expression" :depends-on ("lft" "lft-stream" "package"))
               (:file "cf-stream" :depends-on ("lft" "lft-stream" "package"))
               (:file "lft-arith" :depends-on ("binary-expression" "lft" "lft-stream" "package"))))
