(TeX-add-style-hook
 "gaussian_variance_vs_time_2_standalone"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz")))
   (TeX-run-style-hooks
    "latex2e"
    "gaussian_variance_vs_time_2"
    "standalone"
    "standalone10"
    "pgfplots")
   (LaTeX-add-lengths
    "figurewidth"))
 :latex)

