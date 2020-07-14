;; This module defines Hamiltonian Monte-Carlo (HMC).
; C.f. ... TODO: Add reference.


; vector utils

(define (size vector) (vector-length vector))

(define (reduce-sum vector)
  (vector-fold (lambda (index summation x) (+ summation x)) 0 vector))

(define (elementwise-square vector)
  (define (square x) (* x x))
  (vector-map square vector))

(define (elementwise-add vector-1 vector-2)
  (vector-map (lambda (x y) (+ x y)) vector-1 vector-2))

(define (scalar-product scalar vector)
  (define (prod-s x) (* scalar x))
    (vector-map prod-s vector))


; state utils

(define (make-state p q) (cons p q))
(define (get-p state) (car state))
(define (get-q state) (cdr state))

(define (sample-p size)
  (define (sample x) (- 1 (random 2.0)))
  (vector-map sample (make-vector size 1)))


; dynamics

(define (leap-frog-integrator p+dp/2 q+dq dt steps)
  (define (leap-frog-iter i state)
    (let ((p (get-p state))
          (q (get-q state)))
      (if (eq? i steps)
          (make-state p q)
          (let ((half-p (p+dp/2 p q dt)))
            (let ((next-q (q+dq half-p q dt)))
              (let ((next-p (p+dp/2 half-p next-q dt)))
                (leap-frog-iter (+ i 1) (make-state next-p next-q))))))))
  (lambda (state) (leap-frog-iter 1 state)))

(define (K p) (* 0.5 (reduce-sum (elementwise-square p))))


; main function

(define (HMC U grad-U epsilon L current-q)
  (define (q+dq p q dt)
    (elementwise-add q (scalar-product dt p)))
  (define (p+dp/2 p q dt)
    (elementwise-add p (scalar-product (* -0.5 dt) (grad-U q))))
  (define leap-frog (leap-frog-integrator p+dp/2 q+dq epsilon L))

  (define (energy-diff current-state proposed-state)
    (let ((current-K (K (get-p current-state)))
          (current-U (U (get-q current-state)))
          (proposed-K (K (get-p proposed-state)))
          (proposed-U (U (get-q proposed-state))))
      (- (+ current-K current-U) (+ proposed-K proposed-U))))

  (let ((current-state (make-state (sample-p (size current-q)) current-q)))
    (let ((proposed-state (leap-frog current-state)))
      (if (< (random 1.0)
             (energy-diff current-state proposed-state))
          (get-q proposed-state)
          (get-q current-state)))))


; test

(define main
  (let ((U (lambda (q) (* -0.5 (reduce-sum (elementwise-square q)))))
        (grad-U (lambda (q) (scalar-product -1 q)))
        (current-q (sample-p 20))
        (epsilon 0.1)
        (L 10))
    (HMC U grad-U epsilon L current-q)))