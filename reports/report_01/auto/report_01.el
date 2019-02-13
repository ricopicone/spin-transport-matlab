(TeX-add-style-hook
 "report_01"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("amsart" "onecolumn")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem") ("xcolor" "hyperref" "x11names")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "amsart"
    "amsart10"
    "inputenc"
    "fontenc"
    "graphicx"
    "grffile"
    "longtable"
    "wrapfig"
    "rotating"
    "ulem"
    "amsmath"
    "textcomp"
    "amssymb"
    "capt-of"
    "hyperref"
    "xcolor"
    "cleveref"
    "tikz"
    "subcaption")
   (TeX-add-symbols
    '("Diff" 1))
   (LaTeX-add-labels
    "sec:org50cf7cc"
    "sec:org96530b0"
    "eq:diffusion"
    "eq:solution"
    "eq:solution-dirac"
    "sec:org70efa3f"
    "fig:1a"
    "fig:1b"
    "fig:1"
    "fig:2a"
    "fig:2b"
    "fig:2"
    "sec:orgece42eb")
   (LaTeX-add-bibliographies))
 :latex)

