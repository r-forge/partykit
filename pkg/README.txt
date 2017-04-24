
There are two versions, actually three, of partykit in this folder.

(1) partykitR1: Based on CRAN version 1.1-1 with some minor patches.
    This is essentially just a fall-back option.

(2) partykit: This is the development version and contains two
    different flavours. The first one contains a completely new ctree()
    but the mob() code is equivalent to partykit 1.1-1.
    The second one also contains completely new mob() code.

    To compile the "new ctree / old mob" variant, do:

    - comment "prune.modelparty" in partykit/NAMESPACE
    - have the line ^.*R2.*$ in partykit/.Rbuildignore

    Run R CMD build partykit --no-build-vignettes 
    _FIRST_ and then R CMD build partykit on the 
    extracted tar-ball (R CMD build ignores .Rbuildignore)

    All files (in /R, /man, /vignettes and /tests) with *_R2.*
    corresponding to new mob() code are removed.

    To compile the "new ctree / new mob" variant, do

    - uncomment "prune.modelparty" in partykit/NAMESPACE
    - have the line ^.*R1.*$ in partykit/.Rbuildignore

    Run R CMD build partykit --no-build-vignettes 
    _FIRST_ and then R CMD build partykit on the 
    extracted tar-ball (R CMD build ignores .Rbuildignore)

    All files (in /R, /man, /vignettes and /tests) with *_R1.*
    corresponding to old mob() code are removed.

    By default, the new ctree() / new mob() variant is compiled.
