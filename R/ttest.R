ttest <-
function (x, y = NULL, alternative = c("two.sided", "less", "greater"), 
    mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95, 
    LG=c("NULL","LOG2","LOG","LOG10"),
    ...) 
{
    alternative <- match.arg(alternative)
    LG<-match.arg(LG)
    if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
        stop("'mu' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
        if (paired) 
            xok <- yok <- complete.cases(x, y)
        else {
            yok <- !is.na(y)
            xok <- !is.na(x)
        }
        y <- y[yok]
    }
    else {
        dname <- deparse(substitute(x))
        if (paired) 
            stop("'y' is missing for paired test")
        xok <- !is.na(x)
        yok <- NULL
    }
    x <- x[xok]
    if (paired) {
        x <- x - y
        y <- NULL
    }
    nx <- length(x)
    mx <- mean(x)
    vx <- var(x)
    if (is.null(y)) {
        if (nx < 2) 
            stop("not enough 'x' observations")
        df <- nx - 1
        
        if(LG=="LOG2"){
         stderr <- sqrt(vx/nx+log2(nx/nx-1)^2)
        }else if(LG=="LOG"){
         stderr <- sqrt(vx/nx+log(nx/nx-1)^2)
        }else if(LG=="LOG10"){
        stderr <- sqrt(vx/nx+log10(nx/nx-1)^2)
        }else{
        stderr <- sqrt(vx/nx+(nx/nx-1)^2)
         }

        if (stderr < 10 * .Machine$double.eps * abs(mx)) 
            stop("data are essentially constant")
        tstat <- (mx - mu)/stderr
        method <- if (paired) 
            "Paired t-test"
        else "One Sample t-test"
        estimate <- setNames(mx, if (paired) 
            "mean of the differences"
        else "mean of x")
    }
    else {
        ny <- length(y)
        if (nx < 1 || (!var.equal && nx < 2)) 
            stop("not enough 'x' observations")
        if (ny < 1 || (!var.equal && ny < 2)) 
            stop("not enough 'y' observations")
        if (var.equal && nx + ny < 3) 
            stop("not enough observations")
        my <- mean(y)
        vy <- var(y)
        method <- paste(if (!var.equal) 
            "Welch", "Two Sample t-test")
        estimate <- c(mx, my)
        names(estimate) <- c("mean of x", "mean of y")
        if (var.equal) {
            df <- nx + ny - 2
            v <- 0
            if (nx > 1) 
                v <- v + (nx - 1) * vx
            if (ny > 1) 
                v <- v + (ny - 1) * vy
            v <- v/df

            if(LG=="LOG2"){
            stderr <- sqrt(v * (1/nx + 1/ny)+log2((nx+ny)/(nx+ny-1))^2)
            }else if(LG=="LOG"){
             stderr <- sqrt(v * (1/nx + 1/ny)+log((nx+ny)/(nx+ny-1))^2)
             }else if(LG=="LOG10"){
             stderr <- sqrt(v * (1/nx + 1/ny)+log10((nx+ny)/(nx+ny-1))^2)
             }else{
            stderr <- sqrt(v * (1/nx + 1/ny)+((nx+ny)/(nx+ny-1))^2)
             }

        }
        else {
            stderrx <- sqrt(vx/nx)
            stderry <- sqrt(vy/ny)

            if(LG=="LOG2"){
             stderr <- sqrt(stderrx^2 + stderry^2+log2((nx+ny)/(nx+ny-1))^2)

           }else if(LG=="LOG"){
               stderr <- sqrt(stderrx^2 + stderry^2+log((nx+ny)/(nx+ny-1))^2)

           }else if(LG=="LOG10"){
             stderr <- sqrt(stderrx^2 + stderry^2+log10((nx+ny)/(nx+ny-1))^2)
           }else{
            stderr <- sqrt(stderrx^2 + stderry^2+((nx+ny)/(nx+ny-1))^2)
            }

            df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 
                1))
        }
        if (stderr < 10 * .Machine$double.eps * max(abs(mx), 
            abs(my))) 
            stop("data are essentially constant")
        tstat <- (mx - my - mu)/stderr
    }
    if (alternative == "less") {
        pval <- pt(tstat, df)
        cint <- c(-Inf, tstat + qt(conf.level, df))
    }
    else if (alternative == "greater") {
        pval <- pt(tstat, df, lower.tail = FALSE)
        cint <- c(tstat - qt(conf.level, df), Inf)
    }
    else {
        pval <- 2 * pt(-abs(tstat), df)
        alpha <- 1 - conf.level
        cint <- qt(1 - alpha/2, df)
        cint <- tstat + c(-cint, cint)
    }
    cint <- mu + cint * stderr
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if (paired || !is.null(y)) 
        "difference in means"
    else "mean"
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval, 
        conf.int = cint, estimate = estimate, null.value = mu, 
        alternative = alternative, method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
}
