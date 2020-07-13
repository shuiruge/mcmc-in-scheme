; p, q, and the like, are vectors.
; TODO: implement `length`, `reduce-sum`, `elementwise-square`, `scalar-product`, and `element-add`.
; TODO: implement `random-normal` and `random-uniform`

(define (length vector) ...)
(define (reduce-sum vector) ...)
(define (elementwise-square vector) ...)
(define (elementwise-add vector-1 vector-2) ...)
(define (scalar-product scalar vector) ...)
(define (random-normal size mu sigma) ...)
(define (random-uniform min-value max-value) ...)

(define (make-state p q) (cons p q))
(define (get-p state) (car state))
(define (get-q state) (cdr state))

(define (K p) (* 0.5 (reduce-sum (elementwise-square p))))

(define (hamiltonian-monte-carlo U grad-U epsilon L current-q)
    (define current-p (random-normal (length current-q) 0 1))
    (define current-U (U current-q))
    (define current-K (K current-p))

    (define (p-dynamics p q dt)
        (scalar-product (* -1 dt) (grad-U q)))
    (define (q-dynamics p q dt)
        (scalar-product dt p))
    (define (p+dp p q)
        (elementwise-add p (p-dynamics p q epsilon)))
    (define (q+dq p q dt)
        (elementwise-add q (q-dynamics p q epsilon)))
    (define (p+dp/2 p q)
        (elementwise-add p (p-dynamics p q (* 0.5 epsilon))))
    (define (get-next-state i state)
        (define p (get-p state))
        (define q (get-q state))
        (if (eq? i L)
            (make-pair p (q+dq p q))
            (get-next-state (+ i 1) (make-state (p+dp p q) (q+dq p q)))))

    (let ((next-state
            (get-next-state 1 (make-state (p+dp/2 current-p current-q)
                                          current-q)))
          (p (scalar-product -1 (p+dp/2 (get-p next-state))))
          (q (get-q next-state))
          (proposed-U (U q))
          (proposed-K (K p))
          (alpha (random-uniform 0 1)))
         (if (< alpha (+ current-K current-U (* -1 proposed-K) (* -1 proposed-U)))
             q current-q)))