;;
;; LOH model estimation script.
;; Copyright 2009, R. V. Craiu and A. F. Di Narzo

;;
;; General settings:
;;
(define burnin (* 1000 40)) ;number of burnin batches
(define thin 1)      ;thinning interval
(define batches (* 19000 40)) ;number of batches for each chain
(define starting-values ;chains starting values
  (vector
   #(0.0 -0.666666666666667 -1.2 -14.2857142857143)
   #(-1.0 0.666666666666667 -0.4 -8.57142857142857)
   #(1.5 0.222222222222222 -0.24 -19.1836734693878)
   #(-1.75 1.55555555555556 0.56 -13.469387755102)
   ))

(load "simul_aux.scm") ;load auxiliary stuff
(define dim 4) ;state space dimension

;; Starting mixture estimates:
(define beta-hat-0 (v2gv #(0.8 0.2)))
(define mu-hat-0 (vector-map v2gv
                             #(#(2.2 -1.4 1.4 12.2)
                               #(-2.2 2.2 -1.15 -13.25))))
(define sigma-hat-0 (vector-ec (: i 2) (diag 10.0 dim)))
;;
;; End of settings
;;

;; Load LOH model distribution function
(load-extension
 "./libguiledistribLOH.so"
 "scm_init_distribLOH_module")

;;Target distribution: LOH posterior
(define pi (lambda (x) (distrib-LOH '() x)))

(define (run-RAPTOR startVal rng fout)
  "run a RAPTOR simulation, starting from 'startVal'.
   Outputs the chain to 'fout' in plain txt format.
   rng: random number generator.
   Returns a list of diagnostics"
  (let*
      ((x (v2gv startVal))
       (s (make-amh 'raptor
                    rng ;random number generator
                    pi ;target distribution object
                    x ;current chain value
                    (* burnin thin) ;burnin length
                    (diag 1000.0 dim) ;starting proposal var-cov. matrix
                    beta-hat-0
                    mu-hat-0
                    sigma-hat-0))
       (monitor-base '()))
    (mcmclib-raptor-set-alpha (slot-ref s 'c-ref) 0.3)
    (do-ec (: i (* burnin thin)) (update s)) ;burnin
    (set! monitor-base (make-monitor x)) ;monitor mean, AR, MSJD
    (do-ec (: i batches)
           (begin
             (do-ec (: ii thin)
                    (begin
                      (update s)
                      (update monitor-base)))
             (gsl-vector-fprintf fout x "%f")))
    (free s)
    (list
     (get-monitor-value monitor-base 'means)
     (get-monitor-value monitor-base 'ar)
     (get-monitor-value monitor-base 'msjd))))


(define (run-AM startVal rng fout)
  "run an AM simulation, starting from 'startVal'.
   Outputs the chain to 'fout' in plain txt format.
   rng: random number generator.
   Returns a list of diagnostics"
  (let*
      ((x (v2gv startVal))
       (s (make-amh
           'gauss-am
           rng ;random number generator
           pi ;target distribution object
           x ;current chain value
           (diag 1000.0 dim) ;starting proposal var-cov. matrix
           (* burnin thin)))
       (monitor-base '()))
    (do-ec (: i (* burnin thin)) (update s)) ;burnin
    (set! monitor-base (make-monitor x)) ;monitor mean, AR, MSJD
    (do-ec (: i batches)
           (begin
             (do-ec (: ii thin)
                    (begin
                      (update s)
                      (update monitor-base)))
             (gsl-vector-fprintf fout x "%f")))
    (free s)
    (list
     (get-monitor-value monitor-base 'means)
     (get-monitor-value monitor-base 'ar)
     (get-monitor-value monitor-base 'msjd))))

(define (simulate chain-fun filename)
  (let*
      ((rng (new-gsl-rng (gsl-rng-default)))
       (diags (list
               (make-vector dim 0.0)
               (make-vector dim 0.0)
               (make-vector dim 0.0)
               (make-array 0.0 2 2)))
       (gx (new-gsl-vector dim))
       (fout (fopen filename "w")))
    (do-ec (: i (vector-length starting-values))
           (set! diags (update-means
                        diags
                        (chain-fun (vector-ref starting-values i) rng fout) i)))
    (fclose fout)
    diags))

;about 16mins each on a 2.6Ghz CPU
(define results-RAPTOR (time (simulate run-RAPTOR "chains-RAPTOR.dat")))
(with-output-to-file "LOH-RAPTOR.out" (lambda () (pretty-print results-RAPTOR)))

(define results-AM (time (simulate run-AM "chains-AM.dat")))
(with-output-to-file "LOH-AM.out" (lambda () (pretty-print results-AM)))
