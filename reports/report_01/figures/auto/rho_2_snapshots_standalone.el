(TeX-add-style-hook
 "rho_2_snapshots_standalone"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz")))
   (TeX-run-style-hooks
    "latex2e"
    "rho_2_snapshots"
    "standalone"
    "standalone10"
    "pgfplots")
   (LaTeX-add-lengths
    "figurewidth"))
 :latex)

