PlotGroups2 <- function (data, edesign = NULL, time = edesign[, 1], 
                         groups = edesign[,c(3:ncol(edesign))], repvect = edesign[, 2], 
                         show.fit = FALSE, dis = NULL, step.method = "backward", 
                         min.obs = 2, alfa = 0.05, 
                         nvar.correction = FALSE, summary.mode = "median", show.lines = TRUE, 
                         groups.vector = NULL, xlab = "Time", ylab = "Expression value", 
                         cex.xaxis = 1, ylim = NULL, main = NULL, cexlab = 0.8, legend = TRUE, 
                         sub = NULL, item = NULL) {
  tmp=list()
  if (!is.vector(data)) {
    if (summary.mode == "representative") {
      distances <- apply(as.matrix(dist(data, diag = TRUE, 
                                        upper = TRUE)), 1, sum)
      representative <- names(distances)[distances == 
                                           min(distances)]
      yy <- as.numeric(data[rownames(data) == representative, 
                            ])
      sub <- paste("Representative:", representative)
    }
    else if (summary.mode == "median") {
      yy <- apply(as.matrix(data), 2, median, na.rm = TRUE)
      if (is.null(sub)) {
        sub <- paste("Median profile of", nrow(data), 
                     "features", sep = " ")
      }
    }
    else stop("not valid summary.mode")
    if (dim(data)[1] == 1) {
      sub <- rownames(data)
    }
  }
  else if (length(data) != 0) {
    yy <- as.numeric(data)
    sub <- rownames(data)
  }
  else stop("empty data")
  if (is.null(ncol(groups))) {
    ncol = 1
    legend = FALSE
    codeg = "group"
  }
  else {
    ncol = ncol(groups)
    codeg <- as.character(colnames(groups))
  }
  reps <- i.rank(repvect)
  y <- vector(mode = "numeric", length = length(unique(reps)))
  x <- vector(mode = "numeric", length = length(unique(reps)))
  g <- matrix(nrow = length(unique(reps)), ncol = ncol)
  for (k in 1:length(y)) {
    y[k] <- mean(yy[reps == k], na.rm = TRUE)
    x[k] <- mean(time[reps == k])
    for (j in 1:ncol) {
      g[k, j] <- mean(groups[reps == k, j])
    }
  }
  if (is.null(ylim)) 
    ylim = c(min(as.numeric(yy), na.rm = TRUE), max(as.numeric(yy), 
                                                    na.rm = TRUE))
  abcissa <- x
  xlim = c(min(abcissa, na.rm = TRUE), max(abcissa, na.rm = TRUE))
  color1 <- as.numeric(sort(factor(colnames(groups)))) + 1
  color1_LJJ <- color1
  color1_LJJ <- unlist(lapply(color1_LJJ, function(x){gsub(color1[1],pal_lancet()(9)[1],x)}))
  color1_LJJ <- unlist(lapply(color1_LJJ, function(x){gsub(color1[2],pal_lancet()(9)[7],x)}))
  color2 <- groups
  for (j in 1:ncol) {
    color2[, j] <- color2[, j] * j
  }
  color2 <- as.vector(apply(color2, 1, sum) + 1)
  color2_LJJ <- color2
  color2_LJJ <- unlist(lapply(color2_LJJ, function(x){gsub(color1[1],pal_lancet()(9)[1],x)}))
  color2_LJJ <- unlist(lapply(color2_LJJ, function(x){gsub(color1[2],pal_lancet()(9)[7],x)}))
  
  #plot(x = time, y = yy, pch = 21, col = color2_LJJ, xlab = xlab, ylab = ylab, 
  #     xaxt = "n", main = main, sub = sub, ylim = ylim, xlim = xlim, 
  #     cex = cexlab)
  #***
  plot(x = time, y = yy, pch = 21, col = NA, xlab = xlab, ylab = ylab, 
       xaxt = "n", main = main, sub = sub, ylim = ylim, xlim = xlim, 
       cex = cexlab)
  
  axis(1, at = unique(abcissa), labels = unique(abcissa), 
       cex.axis = cex.xaxis)
  if (show.fit) {
    rm <- matrix(yy, nrow = 1, ncol = length(yy))
    rownames(rm) <- c("ratio medio")
    colnames(rm) <- rownames(dis)
    fit.y <- T.fit(rm, design = dis, step.method = step.method, 
                   min.obs = min.obs, alfa = alfa, nvar.correction = nvar.correction)
    betas <- fit.y$coefficients
  }
  for (i in 1:ncol(groups)) {
    group <- g[, i]
    if ((show.fit) && !is.null(betas)) {
      li <- c(2:6)
      a <- reg.coeffs(coefficients = betas, groups.vector = groups.vector, 
                      group = colnames(groups)[i])
      a <- c(a, rep(0, (7 - length(a))))
      curve(a[1] + a[2] * x + a[3] * (x^2) + a[4] * (x^3) + 
              a[5] * (x^4) + a[6] * (x^5) + a[7] * (x^5), 
            from = min(time), to = max(time), col = color1_LJJ[i], 
            add = TRUE, lty = 2, lwd=2)
    }
    if (show.lines) {
      lx <- abcissa[group != 0]
      ly <- y[group != 0]
      ord <- order(lx)
      lxo <- lx[ord]
      lyo <- ly[ord]
      #***
      points(lxo, lyo, col = color1_LJJ[i])
      lines(lxo, lyo, col = color1_LJJ[i], lwd=2)
      tmp <- c(tmp,list(lxo,lyo))
    }
  }
  op <- par(bg = "white")
  if (legend) 
    legend("topright", legend = codeg, 
           text.col = color1_LJJ, col = color1_LJJ, cex = cexlab, lty = 1, 
           yjust = 0, bty="n")
  par(op)
  
  
  OUTPUT <- list(lxo, lyo, time, yy, tmp)
  names(OUTPUT) <- c("x", "y","point.x","point.y","coordinate")
  OUTPUT
}
see.genes2 <- function (data, edesign = data$edesign, time.col = 1, repl.col = 2, 
                        group.cols = c(3:ncol(edesign)), names.groups = colnames(edesign)[3:ncol(edesign)], 
                        cluster.data = 1, groups.vector = data$groups.vector, k = 9, 
                        k.mclust = FALSE, cluster.method = "hclust", distance = "cor", 
                        agglo.method = "ward.D", show.fit = FALSE, dis = NULL, step.method = "backward", 
                        min.obs = 3, alfa = 0.05, nvar.correction = FALSE, show.lines = TRUE, 
                        iter.max = 500, summary.mode = "median", color.mode = "rainbow", 
                        cexlab = 1, legend = TRUE, newX11 = TRUE, ylim = NULL, main = NULL, 
                        item = "genes", ...) {
  time = edesign[, time.col]
  repvect = edesign[, repl.col]
  groups = edesign[, group.cols]
  narrays <- length(time)
  if (!is.null(dim(data))) {
    dat <- log10(as.data.frame(data)+1)
    clusterdata <- data
  }
  else {
    clusterdata <- data[[cluster.data]]
    dat <- log10(as.data.frame(data$sig.profiles)+1)
  }
  if (nrow(dat) > 1) {
    dat <- as.data.frame(dat[, (ncol(dat) - length(time) + 
                                  1):ncol(dat)])
    count.noNa <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[which(apply(as.matrix(dat), 1, count.noNa) >= 
                       length(unique(repvect))), ]
    #clusterdata <- dat
    if (any(is.na(clusterdata))) {
      if (cluster.method == "kmeans" || cluster.method == 
          "Mclust") {
        if (all(cluster.data != 1, cluster.data != "sig.profiles")) {
          clusterdata[is.na(clusterdata)] <- 0
        }
        else {
          mean.replic <- function(x) {
            tapply(as.numeric(x), repvect, mean, na.rm = TRUE)
          }
          MR <- t(apply(clusterdata, 1, mean.replic))
          if (any(is.na(MR))) {
            row.mean <- t(apply(MR, 1, mean, na.rm = TRUE))
            MRR <- matrix(row.mean, nrow(MR), ncol(MR))
            MR[is.na(MR)] <- MRR[is.na(MR)]
          }
          data.noNA <- matrix(NA, nrow(clusterdata), 
                              ncol(clusterdata))
          u.repvect <- unique(repvect)
          for (i in 1:nrow(clusterdata)) {
            for (j in 1:length(u.repvect)) {
              data.noNA[i, repvect == u.repvect[j]] = MR[i, 
                                                         u.repvect[j]]
            }
          }
          clusterdata <- data.noNA
        }
      }
    }
    if (!is.null(clusterdata)) {
      k <- min(k, nrow(dat), na.rm = TRUE)
      if (cluster.method == "hclust") {
        if (distance == "cor") {
          dcorrel <- matrix(rep(1, nrow(clusterdata)^2), 
                            nrow(clusterdata), nrow(clusterdata)) - 
            cor(t(clusterdata), use = "pairwise.complete.obs")
          clust <- hclust(as.dist(dcorrel), method = agglo.method)
          c.algo.used = paste(cluster.method, "cor", 
                              agglo.method, sep = "_")
        }
        else {
          clust <- hclust(dist(clusterdata, method = distance), 
                          method = agglo.method)
          c.algo.used = paste(cluster.method, distance, 
                              agglo.method, sep = "_")
        }
        cut <- cutree(clust, k = k)
      }
      else if (cluster.method == "kmeans") {
        cut <- kmeans(clusterdata, k, iter.max)$cluster
        c.algo.used = paste("kmeans", k, iter.max, sep = "_")
      }
      else if (cluster.method == "Mclust") {
        if (k.mclust) {
          my.mclust <- Mclust(clusterdata)
          k = my.mclust$G
        }
        else {
          my.mclust <- Mclust(clusterdata, k)
        }
        cut <- my.mclust$class
        c.algo.used = paste("Mclust", k, sep = "_")
      }
      else stop("Invalid cluster algorithm")
      #   if (newX11) 
      #     X11()
      #  groups <- as.matrix(groups)
      #  colnames(groups) <- names.groups
      #   if (k <= 4) 
      #     par(mfrow = c(2, 2))
      #  else if (k <= 6) 
      #    par(mfrow = c(3, 2))
      #  else if (k > 6) 
      #    par(mfrow = c(3, 3))
      #  for (i in 1:(k)) {
      #    PlotProfiles(data = dat[cut == i, ], repvect = repvect, 
      #                 main = i, ylim = ylim, color.mode = color.mode, 
      #                cond = rownames(edesign), item = item, ...)
      # }
      if (newX11) 
        X11()
      if (k <= 4) {
        par(mfrow = c(2, 2))
        cexlab = 0.6
      }
      else if (k <= 6) {
        par(mfrow = c(3, 2))
        cexlab = 0.6
      }
      else if (k > 6) {
        par(mfrow = c(3, 3))
        cexlab = 0.35
      }
      for (j in 1:(k)) {
        PlotGroups2(data = dat[cut == j, ], show.fit = show.fit, 
                    dis = dis, step.method = step.method, min.obs = min.obs, 
                    alfa = alfa, nvar.correction = nvar.correction, 
                    show.lines = show.lines, time = time, groups = groups, 
                    repvect = repvect, summary.mode = summary.mode, 
                    xlab = "time", main = paste("Cluster", j, 
                                                sep = " "), ylim = ylim, cexlab = cexlab, 
                    legend = legend, groups.vector = groups.vector, 
                    item = item, ...)
      }
    }
    else {
      print("warning: impossible to compute hierarchical clustering")
      c.algo.used <- NULL
      cut <- 1
    }
  }
  else if (nrow(dat) == 1) {
    # if (newX11) 
    #   X11()
    # PlotProfiles(data = dat, repvect = repvect, main = NULL, 
    #             ylim = ylim, color.mode = color.mode, cond = rownames(edesign), 
    #             ...)
    if (newX11) 
      X11()
    PlotGroups2(data = dat, show.fit = show.fit, dis = dis, 
                step.method = step.method, min.obs = min.obs, alfa = alfa, 
                nvar.correction = nvar.correction, show.lines = show.lines, 
                time = time, groups = groups, repvect = repvect, 
                summary.mode = summary.mode, xlab = "time", main = main, 
                ylim = ylim, cexlab = cexlab, legend = legend, groups.vector = groups.vector, 
                ...)
    c.algo.used <- NULL
    cut <- 1
  }
  else {
    print("warning: NULL data. No visualization possible")
    c.algo.used <- NULL
    cut <- NULL
  }
  OUTPUT <- list(cut, c.algo.used, groups)
  names(OUTPUT) <- c("cut", "cluster.algorithm.used", "groups")
  OUTPUT
}


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