;;; -*- Lisp -*-

(defsystem "linear-fractional-transformation"
  :depends-on ("fold" "generic-arithmetic" "named-let" "promise" "series" "stream" "utilities")
  :in-order-to ((test-op (test-op "linear-fractional-transformation/tests")))
  :components ((:file "package")
               (:file "lft" :depends-on ("package"))
               (:file "lft-stream" :depends-on ("lft" "package"))
               (:file "binary-expression" :depends-on ("lft" "lft-stream" "package"))
               (:file "cf-stream" :depends-on ("lft" "lft-stream" "package"))
               (:file "lft-arith" :depends-on ("binary-expression" "lft" "lft-stream" "package"))))

(defsystem "linear-fractional-transformation/tests"
  :description "Test suite for linear-fractional-transformation"
  :depends-on ("linear-fractional-transformation" "fiveam")
  :components ((:module "tests"
                :components ((:file "package")
                             (:file "tests" :depends-on ("package"))
                             (:file "arith-tests" :depends-on ("package"))
                             (:file "lft-stream-tests" :depends-on ("package"))
                 (:file "binary-expression-tests" :depends-on ("package"))
                             (:file "cf-stream-tests" :depends-on ("package"))
                             (:file "integration-tests"
                              :depends-on ("package" "tests" "arith-tests"
                                           "lft-stream-tests" "binary-expression-tests"
                                           "cf-stream-tests")))))
  :perform (test-op (op c)
                     (uiop:symbol-call :fiveam :run!
                                       (uiop:find-symbol* :lft-suite :lft-tests))))
