#' @title Estimating Parameters via the Expected-Maximization algorithm.
#' @description
#' This function uses the Expected-Maximization (EM) algorithm to estimate the maximum likelihood estimates (MLE) of the parameters.
#' @param counts_data This is a gene-by-cell raw count matrix.
#' @param group This is a vector of factors with 2 levels that serves as a grouping variable that assigns each column (cell) of the count matrix to one of two predefined experimental conditions.
#' @param CellCluster This is a vector of factors with the optimal number of cluster levels, with each entry assigning a cell from the count matrix's columns to its computationally inferred cluster identity.
#' @param CellAuxil This is a vector of factors of cell-level metadata that provides additional information about each cell, such as its batch, donor, or quality control metrics, with each entry aligning to a column in the count matrix.
#' @param weights This is the observation-wise weight for the cells, with a default value of NULL.
#' @param muoffset This is the offset parameter for mean (mu), with a default value of NULL.
#' @param phioffset This is the offset parameter for the zero inflation (phi) parameter, with a default value of NULL.
#' @param maxit This is the maximum number of iterations for the Expected-Maximization (EM) algorithm.
#' @param eps This is the convergence criteria for the Expected-Maximization (EM) algorithm.
#' @importFrom stats model.matrix
#' @importFrom stats glm.fit
#' @importFrom MASS glm.nb
#' @returns
#' This function returns a list of results.
#' @importFrom stats model.matrix
#' @importFrom stats glm.fit
#' @importFrom MASS glm.nb
#' @export
#' @examples
#' # Do not run.
#' # Paste the example into your R console and run it.
#' # Load the test data.
#' data("TestData")
#' counts_data <- as.matrix(TestData$CountData)
#' group <- c(rep(1,200), rep(2,200))
#' CellCluster <- c(rep(1,60), rep(2,90), rep(3, 80), rep(4, 95), rep(5, 75))
#' group <- as.factor(group)
#' CellCluster <- as.factor(CellCluster)
#' results <- suppressWarnings(SwarnSeq::ZINBEM(counts_data = counts_data, group = group, CellCluster = CellCluster))
ZINBEM <- function(counts_data, group, CellCluster, CellAuxil = NULL, weights = NULL, muoffset = NULL, phioffset = NULL, maxit = NULL, eps = NULL) {
    # Calling all the functions required for this function to run..
    zinb_loglikfun <- function(par, Ycounts, X, Z, weights, offsetx, offsetz) {
        kx <- ncol(X)
        kz <- ncol(Z)
        beta <- par[seq_len(length.out = kx)]
        gamma <- par[(kx + 1):(kx + kz)]
        theta <- exp(par[(kx + kz + 1)])
        eta <- as.vector(X %*% beta + offsetx)
        mu <- exp(eta)
        etaz <- as.vector(Z %*% gamma + offsetz)
        pi0 <- plogis(etaz)
        Y0 <- which(Ycounts == 0)
        Y1 <- which(Ycounts > 0)
        loglik0 <- log(pi0 + exp(log(1 - pi0) + dnbinom(0, mu = mu, size = theta, log = TRUE)))
        loglik1 <- log(1 - pi0) + dnbinom(Ycounts, mu = mu, size = theta, log = TRUE)
        loglik <- -sum(weights[Y0] * loglik0[Y0]) - sum(weights[Y1] * loglik1[Y1])
        return(loglik)
    }
    gradNegBin <- function(par, Ycounts, X, Z, weights, offsetx, offsetz) {
        kx <- ncol(X)
        kz <- ncol(Z)
        beta <- par[seq_len(length.out = kx)]
        gamma <- par[(kx + 1):(kx + kz)]
        theta <- exp(par[(kx + kz + 1)])
        eta <- as.vector(X %*% beta + offsetx)
        mu <- exp(eta)
        etaz <- as.vector(Z %*% gamma + offsetz)
        pi0 <- plogis(etaz)
        Y1 <- as.numeric(Ycounts > 0)
        clogdense0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
        dense0 <- pi0 * (1 - Y1) + exp(log(1 - pi0) + clogdense0)
        wres_count <- ifelse(Y1, Ycounts - mu * (Ycounts + theta) / (mu + theta), -exp(-log(dense0) + log(1 - pi0) + clogdense0 + log(theta) - log(mu + theta) + log(mu)))
        logit_deriv <- make.link("logit")$mu.eta(etaz)
        wres_zero <- ifelse(Y1, -1 / ((1 - pi0) * logit_deriv), (logit_deriv - exp(clogdense0) * logit_deriv) / dense0)
        wres_theta <- theta * ifelse(Y1, digamma(Ycounts + theta) - digamma(theta) + log(theta) - log(mu + theta) + 1 - (Ycounts + theta) / (mu + theta), exp(-log(dense0) + log(1 - pi0) + clogdense0) * (log(theta) - log(mu + theta) + 1 - theta / (mu + theta)))
        gradients_beta <- colSums((wres_count * weights) * X)
        gradients_gamma <- colSums((wres_zero * weights) * Z)
        gradients_theta <- sum(wres_theta * weights)
        return(c(gradients_beta, gradients_gamma, gradients_theta))
    }
    ################### Functions are to be used in the later stages of this process #######################
    show.custom.metod <- function(x) message(x);show.custom.warning <- function(x) warning(x)
    ################### Parameters setup ###################
    group <- as.factor(group);CellCluster <- as.factor(CellCluster)
    if (is.null(CellAuxil)) {
        X <- model.matrix(~ group + CellCluster)
        Z <- model.matrix(~ group + CellCluster)
    }else {
        CellType <- as.factor(CellAuxil)
        X <- model.matrix(~ group + CellCluster + CellType)
        Z <- model.matrix(~ group + CellCluster + CellType)
    }
    n <- ncol(counts_data);m <- nrow(counts_data);kx <- ncol(X);kz <- ncol(Z)
    if (is.null(weights)) {weights <- 1}
    if (length(weights) == 1) {weights <- rep.int(weights, n)}
    weights <- as.vector(weights)
    if (is.null(maxit)) {maxit <- 100}
    if (is.null(eps)) {eps <- 1e-4}
    if (is.null(muoffset)) {muoffset <- 0} # OffsetX
    offsetx <- muoffset
    if (length(offsetx) == 1) {offsetx <- rep.int(offsetx, n)}
    if (is.null(phioffset)) {phioffset <- 0} # OffsetZ
    offsetz <- phioffset
    if (length(offsetz) == 1) {offsetz <- rep.int(offsetz, n)}
    rm(muoffset, phioffset)
    r_names <- rownames(counts_data);c_names <- colnames(counts_data)
    c_clust_length <- length(levels(CellCluster));c_auxil_length <- length(levels(CellAuxil))
    if (!is.null(CellAuxil)) {
        res <- list(parms = matrix(nrow = m, ncol = 3, dimnames = list(r_names, c("Mean", "ZeroInflation", "Dispersion"))), count_group = matrix(nrow = m, ncol = 2, dimnames = list(r_names, c("(Intercept)", "groupInfected"))), count_clust = matrix(nrow = m, ncol = c_clust_length, dimnames = list(r_names, levels(CellCluster))), count_cell_auxil = matrix(nrow = m, ncol = c_auxil_length, dimnames = list(r_names, levels(CellAuxil))), zero_group = matrix(nrow = m, ncol = 2, dimnames = list(r_names, c("(Intercept.Zero)", "groupInfected.Zero"))), zero_clust = matrix(nrow = m, ncol = c_clust_length, dimnames = list(r_names, levels(CellCluster))), zero_cell_auxil = matrix(nrow = m, ncol = c_auxil_length, dimnames = list(r_names, levels(CellAuxil))), `No. of iterations` = vector(length = m), DE.stat = matrix(nrow = m, ncol = 4, dimnames = list(r_names, c("Stat.DE", "Pval.DE", "Stat.DZI", "Pval.DZI"))));names(res$`No. of iterations`) <- r_names
    }else {
        res <- list(parms = matrix(nrow = m, ncol = 3, dimnames = list(r_names, c("Mean", "ZeroInflation", "Dispersion"))), count_group = matrix(nrow = m, ncol = 2, dimnames = list(r_names, c("(Intercept)", "groupInfected"))), count_clust = matrix(nrow = m, ncol = c_clust_length, dimnames = list(r_names, levels(CellCluster))), zero_group = matrix(nrow = m, ncol = 2, dimnames = list(r_names, c("(Intercept.Zero)", "groupInfected.Zero"))), zero_clust = matrix(nrow = m, ncol = c_clust_length, dimnames = list(r_names, levels(CellCluster))), `No. of iterations` = vector(length = m), DE.stat = matrix(nrow = m, ncol = 4, dimnames = list(r_names, c("Stat.DE", "Pval.DE", "Stat.DZI", "Pval.DZI"))));names(res$`No. of iterations`) <- r_names
    }
    glm.nb.error <- vector();glm.nb.limit.reach.error <- vector()
    for (i in seq_len(length.out = m)) {
        show.custom.metod(c("SwarnSeq is analyzing the gene at index ",i,"."))
        Y <- counts_data[i,];Y0 <- Y == 0;Y1 <- Y > 0
        # EM function setup..
        model_count <- glm.fit(x = X, y = Y, family = poisson(), weights = weights, offset = offsetz)
        model_zero <- glm.fit(x = Z, y = Y0, family = binomial(link = "logit"), weights = weights, offset = offsetz)
        model_zero2 <- model_zero
        start <- list(count = model_count$coefficients, zero = model_zero$coefficients, theta = 1)
        start$count[is.na(start$count)] <- 0;start$zero[is.na(start$zero)] <- 0
        # Fitted values
        mui <- model_count$fitted.values;probi <- model_zero$fitted.values
        probi <- probi / (probi + (1 - probi) * dnbinom(0, size = start$theta, mu = mui))
        probi[Y1] <- 0;rownames(probi) <- rownames(Y)
        ll_new <- zinb_loglikfun(par = c(model_count$coefficients, model_zero$coefficients, log(start$theta)), Ycounts = Y, X = X, Z = Z, weights = weights, offsetx = offsetx, offsetz = offsetz)
        ll_old <- 2 * ll_new;iter <- 0;Convergence <- TRUE;start2 <- start
        tryCatch(
            expr = {
                while (abs(ll_old - ll_new) > eps) {
                    ll_old <- ll_new
                    model_count <- MASS::glm.nb(Y ~ 0 + X + offset(offsetx), weights = weights * (1 - probi), start = rep(0, kx), init.theta = start$theta)
                    model_zero <- glm.fit(x = Z, y = probi, family = binomial(link = "logit"), weights = weights, offset = offsetz)
                    start <- list(count = model_count$coefficients, zero = model_zero$coefficients, theta = model_count$theta)
                    mui <- model_count$fitted.values;probi <- model_zero$fitted.values
                    probi <- probi / (probi + (1 - probi) * dnbinom(0, size = start$theta, mu = mui))
                    probi[Y1] <- 0
                    ll_new <- zinb_loglikfun(par = c(model_count$coefficients, model_zero$coefficients, log(start$theta)), Ycounts = Y, X = X, Z = Z, weights = weights, offsetx = offsetx, offsetz = offsetz)
                    iter <- iter + 1;Convergence <- TRUE
                    if (iter > maxit) {
                        Convergence <- FALSE
                        show.custom.metod(c("Convergence was not achieved in maxit."))
                        break
                    }
                }
                start$iter <- iter;start$Convergence <- Convergence;start$loglik_diff <- ll_old - ll_new
            },
            error = function(e) {
                show.custom.warning("Exiting MASS::glm.nb method..")
                show.custom.metod(c("An error occured: ", e$message))
            }
        )
        if (!is.na(ll_new) && !is.na(ll_old)) {
            if (abs(ll_old - ll_new) <= eps & is.null(start$iter)) {
                start <- start
                start$iter <- iter
            }else if (abs(ll_old - ll_new) <= eps) {
                start <- start
                start$iter <- start$iter
            }else if (is.null(start$loglik_diff)) {
                model_count <- glm.fit(x = X, y = Y, family = poisson(), weights = weights, offset = offsetx)
                c <- c(model_count$coefficients, model_zero2$coefficients)
                fit <- optim(par = c(c,1), fn = function(par) zinb_loglikfun(par, Y, X, Z, weights, offsetx, offsetz), gr = function(par) gradNegBin(par, Y, X, Z, weights, offsetx, offsetz), method = "BFGS", control = list(), hessian = FALSE)
                start <- list()
                start$count <- fit$par[seq_len(length.out = kx)];rownames(start$count) <- rownames(Y)
                start$zero <- fit$par[(kx + 1):(kx + kz)];rownames(start$zero) <- rownames(Y)
                start$theta <- 1;start$iter <- 0;glm.nb.error <- append(glm.nb.error, i)
            }else if (!start$Convergence) {
                glm.nb.limit.reach.error <- append(glm.nb.limit.reach.error, i)
                model_count <- glm.fit(x = X, y = Y, family = poisson(), weights = weights, offset = offsetx)
                c <- c(model_count$coefficients, model_zero2$coefficients)
                fit <- optim(par = c(c,1), fn = function(par) zinb_loglikfun(par, Y, X, Z, weights, offsetx, offsetz), gr = function(par) gradNegBin(par, Y, X, Z, weights, offsetx, offsetz), method = "BFGS", control = list(), hessian = FALSE)
                start <- list();start$count <- fit$par[seq_len(length.out = kx)];rownames(start$count) <- rownames(Y)
                start$zero <- fit$par[(kx + 1):(kx + kz)];rownames(start$zero) <- rownames(Y)
                start$theta <- 1;start$iter <- 0
            }
        }else {
            start <- start2;start$iter <- 0;glm.nb.error <- append(glm.nb.error, i)
        }
        # Count p-values
        mu <- exp(X %*% start$count + offsetx)
        mean.finite <- function(x) mean(x[is.finite(x)])
        pop_mean <- mean.finite(mu)
        phi <- plogis(Z %*% start$zero + offsetz)
        pop_phi <- mean.finite(phi)
        theta <- start$theta
        # Statistical testing for DE
        start$count0 <- start$count
        start$count0[2] <- 0
        mu0 <- exp(X %*% start$count0 + offsetx)
        likH0 <- sum(dzinb(Y, size = start$theta, mu = mu0, rho = phi, log = TRUE))
        likH <- sum(dzinb(Y, size = start$theta, mu = mu, rho = phi, log = TRUE))
        lrt.stat <- -2 * (likH0 - likH)
        lrt.stat <- ifelse(lrt.stat < 0, 0, lrt.stat)
        pval <- exp(pchisq(lrt.stat, df = 1, lower.tail = FALSE, log.p = TRUE))
        # Differential zero inflation
        start$zero0 <- start$zero
        # Under null hypothesis
        start$zero0[2] <- 0
        phi0 <- plogis(Z %*% start$zero0 + offsetz)
        phi.likH0 <- sum(dzinb(Y, size = start$theta, mu = mu, rho = phi0, log = TRUE))
        phi.likH <- sum(dzinb(Y, size = start$theta, mu = mu, rho = phi, log = TRUE))
        lrt.stat.phi <- -2 * (phi.likH0 - phi.likH)
        lrt.stat.phi <- ifelse(lrt.stat.phi < 0, 0, lrt.stat.phi)
        pval.phi <- exp(pchisq(lrt.stat.phi, df = 1, lower.tail = FALSE, log.p = TRUE))
        # DE stat
        DE.stat <- c(lrt.stat, pval, lrt.stat.phi, pval.phi)
        # Distribution characterization
        params <- c(pop_mean, pop_phi, start$theta)
        count_group <- start$count[seq_len(length.out = 2)];iter <- start$iter
        ############ Returning the results ##############
        if (!is.null(CellAuxil)) {
            nclust <- length(levels(CellCluster))
            count_clust <- start$count[-((nclust + 2):length(start$count))][-2]
            zero_group <- start$zero[seq_len(length.out = 2)]
            zero_clust <- start$zero[-((nclust + 2):length(start$zero))][-2]
            ncellauxil <- length(levels(CellAuxil))
            count_cell_auxil <- start$count[-(2:(nclust + 2))]
            zero_cell_auxil <- start$zero[-(2:(nclust + 2))]
            res$parms[i,] <- c(params)
            res$count_group[i,] <- count_group
            res$count_clust[i,] <- count_clust
            res$count_cell_auxil <- count_cell_auxil
            res$zero_group[i,] <- zero_group
            res$zero_clust[i,] <- zero_clust
            res$zero_cell_auxil <- zero_cell_auxil
            res$`No. of iterations`[i] <- iter
            res$DE.stat[i,] <- DE.stat
        }else {
            nclust <- length(levels(CellCluster))
            count_clust <- start$count[-2]
            zero_group <- start$zero[seq_len(length.out = 2)]
            zero_clust <- start$zero[-2]
            res$parms[i,] <- c(params)
            res$count_group[i,] <- count_group
            res$count_clust[i,] <- count_clust
            res$zero_group[i,] <- zero_group
            res$zero_clust[i,] <- zero_clust
            res$`No. of iterations`[i] <- iter
            res$DE.stat[i,] <- DE.stat
        }
    }
    rm(params, count_group, count_clust, zero_group, zero_clust, iter)
    res$glm.nb.error.gene.index <- glm.nb.error
    res$glm.nb.limit.reach.error.gene.index <- glm.nb.limit.reach.error
    show.custom.metod(c(" "))
    show.custom.metod(c("\nThis process has been completed successfully.\n"))
    return(res)
}
