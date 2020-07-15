(import (mcmc hmc))

(define main
  (let ((U (lambda (q) (* -0.5 (reduce-sum (elementwise-square q)))))
        (grad-U (lambda (q) (scalar-product -1 q)))
        (current-q (sample-p 20))
        (epsilon 0.1)
        (L 10))
    (HMC U grad-U epsilon L current-q)))