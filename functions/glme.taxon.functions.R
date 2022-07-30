
library("lme4")
library("lmerTest")
library(stats)
library(rstatix)
library(ggplot2)
library(reshape2)
library(export)
library(ggpubr)
library(ggsci)
library("pheatmap")

my_theme <- theme_bw(base_rect_size = 0.8)+
  theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank(), 
        text=element_text(size=14), axis.ticks.x = element_blank())

my_theme2 <- theme(text = element_text(color = "black",  size=14),
                   axis.text = element_text(color = "black",  size=12),
                   axis.line = element_line(size = 0.6),
                   axis.ticks.length  = unit(1, 'mm'),
                   axis.ticks = element_line(size = 0.6),
                   axis.text.x = element_text(vjust=0.5),
                   plot.title = element_text(size=12),
                   plot.subtitle = element_text(size=8))


######################### 从数据中subset提取数据
# dis <- "MM","ALL","LM"
# lev <- "p","c","o","f","g","s"
sub1 <- function(data0, dis, lev){
  if(lev=="p"){
    data1 <- data0[,c(1:8,grep("p__",colnames(data0)))]
    data <- data1[,c(grep("c__",colnames(data1), invert = T))]
  }else if(lev=="c"){
    data1 <- data0[,c(1:8,grep("c__",colnames(data0)))]
    data <- data1[,c(grep("o__",colnames(data1), invert = T))]
  }else if(lev=="o"){
    data1 <- data0[,c(1:8,grep("o__",colnames(data0)))]
    data <- data1[,c(grep("f__",colnames(data1), invert = T))]
  }else if(lev=="f"){
    data1 <- data0[,c(1:8,grep("f__",colnames(data0)))]
    data <- data1[,c(grep("g__",colnames(data1), invert = T))]
  }else if(lev=="g"){
    data1 <- data0[,c(1:8,grep("g__",colnames(data0)))]
    data <- data1[,c(grep("s__",colnames(data1), invert = T))]
  }else if(lev=="s"){
    data1 <- data0[,c(1:8,grep("s__",colnames(data0)))]
    data <- data1
  }
  data_MM <- data[which(data$disease1==dis & data$Stage2 != "CY"),]
  return(data_MM)
}




######################### glmer test

taxa.lmer <- function(data2){
  res <- list()
  for (k in colnames(data2)[-c(1:8)]){
    tmp <- data2[,c("subject","Stage2","Outcome1",k)]
    #tmp[,4] <- log2(tmp[,4]+1)
    fit <- try(re <- glmer(data=tmp, get(k) ~ Outcome1+ (1|Stage2) + (1|subject), family = "poisson"))
  # fit <- try(re <- glmer(data=tmp, get(k) ~ Outcome1+ (1|Stage2) + (1|subject), family=negative.binomial(theta = 10)))
    if("try-error" %in% class(fit)){
      print(k)
      }else{
      re2 <- summary(re)
      res[[eval(k)]] <- re2
    }
  }
  return(res)
}


taxa.lmer11 <- function(data2){
  res <- list()
  for (k in colnames(data2)[-c(1:8)]){
    tmp <- data2[,c("subject","Stage2","Outcome1",k)]
    #tmp[,4] <- log2(tmp[,4]+1)
    fit <- try(re <- lmer(data=tmp, get(k) ~ Outcome1 + (1|Stage2) + (1|subject)))
    # fit <- try(re <- glmer(data=tmp, get(k) ~ Outcome1+ (1|Stage2) + (1|subject), family=negative.binomial(theta = 10)))
    if("try-error" %in% class(fit)){
      print(k)
    }else{
      re2 <- summary(re)
      res[[eval(k)]] <- re2
    }
  }
  return(res)
}


######################### glmer ~ CRSgrade

taxa.lmer2 <- function(data2){
  res <- list()
  for (k in colnames(data2)[-c(1:8)]){
    tmp <- data2[,c("subject","Stage2","CRSgrade",k)]
   # tmp[,4] <- log2(tmp[,4]+1)
    fit <- try(re <- glmer(data=tmp, get(k) ~ CRSgrade+ (1|Stage2) + (1|subject), family = "poisson"))
    # fit <- try(re <- glmer(data=tmp, get(k) ~ CRSgrade+ (1|Stage2) + (1|subject), family=negative.binomial(theta = 10)))
    if("try-error" %in% class(fit)){
      print(k)
    }else{
      re2 <- summary(re)
      res[[eval(k)]] <- re2
    }
  }
  return(res)
}

taxa.lmer22 <- function(data2){
  res <- list()
  for (k in colnames(data2)[-c(1:8)]){
    tmp <- data2[,c("subject","Stage2","CRSgrade",k)]
    #tmp[,4] <- log2(tmp[,4]+1)
    fit <- try(re <- lmer(data=tmp, get(k) ~ CRSgrade + (1|Stage2) + (1|subject)))
    # fit <- try(re <- glmer(data=tmp, get(k) ~ Outcome1+ (1|Stage2) + (1|subject), family=negative.binomial(theta = 10)))
    if("try-error" %in% class(fit)){
      print(k)
    }else{
      re2 <- summary(re)
      res[[eval(k)]] <- re2
    }
  }
  return(res)
}


######################### glmer ~ CRSgrade

taxa.lmer3 <- function(data2){
  res <- list()
  for (k in colnames(data2)[-c(1:8)]){
    tmp <- data2[,c("subject","Stage2","CRSgrade",k)]
   # tmp[,4] <- log2(tmp[,4]+1)
    fit <- try(re <- glmer(data=tmp, get(k) ~ CRSgrade + (1|subject), family = "poisson"))
    # fit <- try(re <- glmer(data=tmp, get(k) ~ CRSgrade+ (1|Stage2) + (1|subject), family=negative.binomial(theta = 10)))
    if("try-error" %in% class(fit)){
      print(k)
    }else{
      re2 <- summary(re)
      res[[eval(k)]] <- re2
    }
  }
  return(res)
}


taxa.lmer33 <- function(data2){
  res <- list()
  for (k in colnames(data2)[-c(1:8)]){
    tmp <- data2[,c("subject","Stage2","CRSgrade",k)]
    #tmp[,4] <- log2(tmp[,4]+1)
    fit <- try(re <- lm(data=tmp, get(k) ~ CRSgrade ))
    # fit <- try(re <- glmer(data=tmp, get(k) ~ Outcome1+ (1|Stage2) + (1|subject), family=negative.binomial(theta = 10)))
    if("try-error" %in% class(fit)){
      print(k)
    }else{
      re2 <- summary(re)
      res[[eval(k)]] <- re2
    }
  }
  return(res)
}


######################### extract one level of taxon id 
exid <- function(x, lev="g"){
  if(lev=="g"){
    y <- lapply(x, function(x){strsplit(x, split="__")[[1]][7]}) %>% unlist()
  }else if(lev=="f"){
    y <- lapply(x, function(x){strsplit(x, split="__")[[1]][6]}) %>% lapply(.,function(x){substr(x, start = 1, stop = nchar(x)-2)}) %>% unlist()
  }else if(lev=="o"){
    y <- lapply(x, function(x){strsplit(x, split="__")[[1]][5]}) %>% lapply(.,function(x){substr(x, start = 1, stop = nchar(x)-2)}) %>% unlist()
  }else if(lev=="c"){
    y <- lapply(x, function(x){strsplit(x, split="__")[[1]][4]}) %>% lapply(.,function(x){substr(x, start = 1, stop = nchar(x)-2)}) %>% unlist()
  }else if(lev=="p"){
    y <- lapply(x, function(x){strsplit(x, split="__")[[1]][3]}) %>% lapply(.,function(x){substr(x, start = 1, stop = nchar(x)-2)}) %>% unlist()
  }
  return(y)
}







################# remove rare (>90% 0)
rare.remove <- function(x){
  k <- apply(x==0, 2, sum)
  n <- nrow(x)*0.8
  x2 <- x[,-which(k > n)]
  return(x2)
}



################# masigpro
T.fit2 <- function (data, design = data$dis, step.method = "backward", 
                    min.obs = data$min.obs, alfa = data$Q, nvar.correction = FALSE, 
                    family = gaussian(), epsilon = 1e-05, item = "gene") 
{
  if (is.list(data)) {
    dat <- as.matrix(data$SELEC)
    dat <- rbind(c(rep(1, ncol(dat))), dat)
    groups.vector <- data$groups.vector
    groups.vector <- c(groups.vector[nchar(groups.vector) == 
                                       min(nchar(groups.vector))][1], groups.vector)
    edesign <- data$edesign
    G <- data$g
    family <- data$family
  }
  else {
    G <- nrow(data)
    data <- rbind(c(rep(1, ncol(data))), data)
    dat <- as.matrix(data)
    count.na <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[apply(dat, 1, count.na) >= min.obs, ]
    groups.vector = NULL
    edesign = NULL
  }
  dis <- as.data.frame(design)
  dat <- dat[, as.character(rownames(dis))]
  g <- (dim(dat)[1] - 1)
  n <- dim(dat)[2]
  p <- dim(dis)[2]
  vars.in <- colnames(dis)
  sol <- coefficients <- group.coeffs <- t.score <- sig.profiles <- NULL
  influ.info <- matrix(NA, nrow = nrow(dis), ncol = 1)
  rownames(influ.info) <- rownames(dis)
  if (nvar.correction) 
    alfa <- alfa/ncol(dis)
  for (i in c(2:83,85:122,124:(g + 1))) {
    y <- as.numeric(dat[i, ])
    name <- rownames(dat)[i]
    if (step.method == "backward") {
      reg <- stepback(y = y, d = dis, alfa = alfa, family = family,  
                      epsilon = epsilon)
    }
    else if (step.method == "forward") {
      reg <- stepfor(y = y, d = dis, alfa = alfa, family = family, 
                     epsilon = epsilon)
    }
    else if (step.method == "two.ways.backward") {
      reg <- two.ways.stepback(y = y, d = dis, alfa = alfa, 
                               family = family, epsilon = epsilon)
    }
    else if (step.method == "two.ways.forward") {
      reg <- two.ways.stepfor(y = y, d = dis, alfa = alfa, 
                              family = family, epsilon = epsilon)
    }
    else stop("stepwise method must be one of backward, forward, two.ways.backward, two.ways.forward")
    div <- c(1:round(g/100)) * 100
    if (is.element(i, div)) 
      print(paste(c("fitting ", item, i, "out of", g), 
                  collapse = " "))
    lmf <- glm(y ~ ., data = as.data.frame(dis), family = family, 
               epsilon = epsilon)
    result <- summary(lmf)
    novar <- vars.in[!is.element(vars.in, names(result$coefficients[, 
                                                                    4]))]
    influ <- influence.measures(reg)$is.inf
    influ <- influ[, c(ncol(influ) - 3, ncol(influ) - 1)]
    influ1 <- which(apply(influ, 1, all))
    if (length(influ1) != 0) {
      paste.names <- function(a) {
        paste(names(a)[a], collapse = "/")
      }
      match <- match(rownames(dis), rownames(influ))
      influ <- as.data.frame(apply(influ, 1, paste.names))
      influ.info <- cbind(influ.info, influ[match, ])
      colnames(influ.info)[ncol(influ.info)] <- name
    }
    result <- summary(reg)
    if ((!(result$aic == -Inf) & !is.na(result$aic) & family$family == 
         "gaussian") | family$family != "gaussian") {
      k <- i
      model.glm.0 <- glm(y ~ 1, family = family, epsilon = epsilon)
      if (family$family == "gaussian") {
        test <- anova(model.glm.0, reg, test = "F")
        p.value = test[6][2, 1]
      }
      else {
        test <- anova(model.glm.0, reg, test = "Chisq")
        p.value = test[5][2, 1]
      }
      bondad <- (reg$null.deviance - reg$deviance)/reg$null.deviance
      if (bondad < 0) {
        bondad = 0
      }
      beta.coeff <- result$coefficients[, 1]
      beta.p.valor <- result$coefficients[, 4]
      coeff <- rep(0, (length(vars.in) + 1))
      if (length(novar) != 0) {
        for (m in 1:length(novar)) {
          coeff[position(dis, novar[m]) + 1] <- NA
        }
      }
      p.valor <- t <- as.numeric(rep(NA, (length(vars.in) + 
                                            1)))
      if (result$coefficients[, 4][rownames(result$coefficients) == 
                                   "(Intercept)"] < alfa) {
        coeff[1] <- result$coefficients[, 1][rownames(result$coefficients) == 
                                               "(Intercept)"]
        p.valor[1] <- result$coefficients[, 4][rownames(result$coefficients) == 
                                                 "(Intercept)"]
        t[1] <- result$coefficients[, 3][rownames(result$coefficients) == 
                                           "(Intercept)"]
      }
      for (j in 2:length(coeff)) {
        if (is.element(vars.in[j - 1], rownames(result$coefficients))) {
          coeff[j] <- result$coefficients[, 1][rownames(result$coefficients) == 
                                                 vars.in[j - 1]]
          p.valor[j] <- result$coefficients[, 4][rownames(result$coefficients) == 
                                                   vars.in[j - 1]]
          t[j] <- result$coefficients[, 3][rownames(result$coefficients) == 
                                             vars.in[j - 1]]
        }
      }
      if (!all(is.na(p.valor))) {
        sol <- rbind(sol, as.numeric(c(p.value, bondad, 
                                       p.valor)))
        coefficients <- rbind(coefficients, coeff)
        t.score <- rbind(t.score, t)
        sig.profiles <- rbind(sig.profiles, y)
        h <- nrow(sol)
        rownames(sol)[h] <- name
        rownames(coefficients)[h] <- name
        rownames(t.score)[h] <- name
        rownames(sig.profiles)[h] <- name
      }
    }
  }
  if (!is.null(sol)) {
    sol <- as.data.frame(sol)
    coefficients <- as.data.frame(coefficients)
    coeffic <- coefficients
    t.score <- as.data.frame(t.score)
    sig.profiles <- as.data.frame(sig.profiles)
    colnames(sol) <- c("p-value", "R-squared", "p.valor_beta0", 
                       paste("p.valor_", vars.in, sep = ""))
    colnames(coefficients) <- c("beta0", paste("beta", vars.in, 
                                               sep = ""))
    colnames(t.score) <- c("t.score_beta0", paste("t.score_", 
                                                  vars.in, sep = ""))
    colnames(sig.profiles) <- colnames(dat)
    if (!is.null(groups.vector) & !is.null(edesign)) {
      groups <- colnames(edesign)[3:ncol(edesign)]
      degree <- (length(groups.vector)/length(groups)) - 
        1
      for (w in 1:nrow(coefficients)) {
        A <- NULL
        col.names <- NULL
        for (l in 1:length(groups)) {
          B <- reg.coeffs(coefficients = coefficients[w, 
                                                      ], groups.vector = groups.vector, group = groups[l])
          cols <- paste(rep(groups[l], each = length(B)), 
                        paste("beta", c(0:(length(B) - 1)), sep = ""), 
                        sep = "_")
          A <- c(A, B)
          col.names <- c(col.names, cols)
        }
        group.coeffs <- (rbind(group.coeffs, A))
      }
      colnames(group.coeffs) <- col.names
      rownames(group.coeffs) <- rownames(coefficients)
    }
  }
  if (ncol(influ.info) > 2) {
    print(paste("Influence:", ncol(influ.info) - 1, "genes with influential data at slot influ.info. Model validation for these genes is recommended"))
  }
  influ.info <- influ.info[, -1]
  output <- list(sol, sig.profiles, coefficients, as.data.frame(group.coeffs), 
                 t.score, vars.in, G, g, dat, dis, step.method, groups.vector, 
                 edesign, influ.info)
  names(output) <- c("sol", "sig.profiles", "coefficients", 
                     "group.coeffs", "t.score", "variables", "G", "g", "dat", 
                     "dis", "step.method", "groups.vector", "edesign", "influ.info")
  output
}

