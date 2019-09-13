
std_beta <- function(Z, N) {
    Z / sqrt(N + Z^2)
}

std_se <- function(Z, N) {
    1 / sqrt(N + Z^2)
}

V_MMR_delta <- function(G, g, C, Sigma) {
    C_inv <- solve(C)
    GCG_inv <- solve(t(G) %*% solve(C) %*% G)
    
    df_dg <- GCG_inv %*% t(G) %*% C_inv
    df_dG <- (GCG_inv %x% (t(g) %*% C_inv %*% (-G %*% GCG_inv %*% t(G) %*% C_inv + diag(nrow(G))))) +
        ((-t(g) %*% C_inv %*% G %*% GCG_inv) %x% (GCG_inv %*% t(G) %*% C_inv))
    J <- cbind(df_dG, df_dg)
    
    V <- J %*% Sigma %*% t(J)
    
    return(V)
}

V_MMR_delta2 <- function(G, g, C, Sigma) {
    require(matrixcalc)
    
    M <- nrow(G)
    K <- ncol(G)
    C_inv <- solve(C)
    GCG_inv <- solve(t(G) %*% solve(C) %*% G)
    
    df_dg <- GCG_inv %*% t(G) %*% C_inv
    df_dG <- (t(C_inv %*% g) %x% diag(K)) %*% ((G %x% diag(K)) %*% (-GCG_inv %x% GCG_inv) %*% ((t(C_inv %*% G) %x% diag(K)) %*% commutation.matrix(M, K) + (diag(K) %x% t(G)) %*% (diag(K) %x% C_inv)) + (diag(M) %x% GCG_inv) %*% commutation.matrix(M, K))
    J <- cbind(df_dG, df_dg)
    
    V <- J %*% Sigma %*% t(J)
    
    return(V)
}

#### TOY example with 3 SNPs and 2 genes ####
set.seed(1234)
n <- 1000
snp1 <- sample(0:2, n, replace = TRUE, prob = c(0.7^2, 2*0.3*0.7, 0.3^2))
snp2 <- sample(0:2, n, replace = TRUE, prob = c(0.7^2, 2*0.3*0.7, 0.3^2))
snp3 <- sample(0:2, n, replace = TRUE, prob = c(0.7^2, 2*0.3*0.7, 0.3^2))
expr1 <- 0.2*snp1 + rnorm(n)
expr2 <- 0.2*snp2 + 0.2*snp3 + rnorm(n)
trait <- 0.2*expr1 + 0.2*expr2 + rnorm(n)

Zscores <- c(summary(lm(expr1 ~ snp1))$coefficients[2, 3], summary(lm(expr1 ~ snp2))$coefficients[2, 3], summary(lm(expr1 ~ snp3))$coefficients[2, 3],
             summary(lm(expr2 ~ snp1))$coefficients[2, 3], summary(lm(expr2 ~ snp2))$coefficients[2, 3], summary(lm(expr2 ~ snp3))$coefficients[2, 3],
             summary(lm(trait ~ snp1))$coefficients[2, 3], summary(lm(trait ~ snp2))$coefficients[2, 3], summary(lm(trait ~ snp3))$coefficients[2, 3])
betas <- std_beta(Zscores, n)
SEs <- std_se(Zscores, n)

G <- matrix(betas[1:6], nrow = 3, ncol = 2, dimnames = list(paste0("snp", 1:3), paste0("expr", 1:2)))
g <- matrix(betas[7:9], nrow = 3, ncol = 1, dimnames = list(paste0("snp", 1:3), "trait"))
C <- diag(3) # We have no correlation between SNPs
R <- diag(2+1) # For simplicity, we assume no gene-gene and gene-trait correlations
Sigma <- (SEs %*% t(SEs)) * (C %x% R)

V_MMR_delta(G, g, C, Sigma)
