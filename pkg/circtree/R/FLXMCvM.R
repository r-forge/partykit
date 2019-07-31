## Flexmix driver von Bettina Gruen
#
# Der flexmix Driver nimmt an, dass die Daten als Winkel auf [0, 2pi] kommen und 
# wandelt diese dann um auf (x, y) Koordinaten am Einheitskreis. Es wird dann der 
# Parameter theta zurueckgegeben, der erwarteter Richtung mal kappa entspricht.
FLXMCvM <- function(formula = .~., solve_kappa = movMF:::solve_kappa_Newton_Fourier) {
    z <- methods::new("FLXMC", weighted = TRUE, formula = formula,
             dist = "vM", name = "model-based von Mises clustering")
    z@preproc.y <- function(x) {
        if (ncol(x) > 1)
            stop("for the von Mises distribution y must be univariate")
        if (any(abs(x) > pi))
            stop("values must be in [-pi, pi]")
        cbind(cos(x), sin(x))
    }
    
    z@defineComponent <- function(para) {
        predict <-  function(x, ...)
            matrix(para$theta, nrow = nrow(x), ncol = 2, byrow = TRUE)

        logLik <- function(x, y)
            movMF::dmovMF(y, theta = para$theta, log = TRUE)

        methods::new("FLXcomponent", parameters = list(theta = para$theta),
            predict = predict, logLik = logLik, df = para$df)
    }
    
    z@fit <- function(x, y, w, ...) {
        ybar <- colMeans(w * y) / mean(w)
        rbar <- sqrt(sum(ybar^2))
        kappa <- solve_kappa(rbar, 2)
        z@defineComponent(list(theta = ybar / rbar * kappa, df = 2))
    }
    z
}


