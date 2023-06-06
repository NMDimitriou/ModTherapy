kdeggpairs <- function(x, n_1d=20, n_2d=25, labels, density=TRUE,
                             contour=TRUE, ...) {
  require(MASS)
  require(sp)
  require(ggplot2)
  require(bayestestR)
  require(dplyr)
  
  par(cex.axis=0.95) # magnify the axes font
  
  
  fun.lower <- function(x1, x2, ...) {
    if (is.factor(x1)) x1 <- as.integer(x1)
    if (is.factor(x2)) x1 <- as.integer(x2)
    OK <- length(unique(x1)) > 2 && length(unique(x2)) > 2
    
    if (!density && !contour) n_2d <- 0
    
    if (n_2d > 0 && OK) {
      if (density || contour) {
        d <- tryCatch(kde2d(x1, x2, n=n_2d), error=function(e) {
          if (e$message == "bandwidths must be strictly positive") {
            cat("Warning: resetting the bandwidths")
            hx <- 0.1*(max(x1)-min(x1)) #had 0.1*
            hy <- 0.1*(max(x2)-min(x2))
            d <- kde2d(x1, x2, n=n_2d, h=c(hx, hy))
          } else cat("Error caught! Handling not implemented!")
        })
      }
      
      if (density) {
        # library("colorspace"); vv<-c("white", heat_hcl(12))
        # library(RColorBrewer); vv<-c("white", rev(brewer.pal(9, "YlOrRd")))
        #vv<-c(rev(bpy.colors(n=10, cutoff.tails=0.3, alpha=0.7))) #"white", 
        vv<-c( bpy.colors(n=100, cutoff.tails=0.3, alpha=0.8)) # 
        
        image(d, col=vv, add=TRUE)
      }
      
      if (contour) graphics:::contour(d, add=TRUE, nlevels=8) #nlevels=8
    } else points(x1, x2)
  }
  
  fun.upper <- function(x1, x2, ...) {
    if (is.factor(x1)) x1 <- as.integer(x1)
    if (is.factor(x2)) x1 <- as.integer(x2)
    
    vv<-bpy.colors(n=100, cutoff.tails=0.3, alpha=0.2)
    
    # This adds a column of color values based on the last column values
    #datc <- vv[as.numeric(cut(as.numeric(loglh), breaks=100))]
    datc <- vv[cut(loglh, breaks=100)]
    points(x1, x2, col=datc, pch=16, cex=0.5)
  }
  
  fun.diag.kde <- function(x, ...) {
    
    posterior <- x
    ci_eti <- ci(posterior,ci=0.95, method = "ETI",verbose = TRUE)
    # Plot the distribution and add the limits of the two CIs
    posterior %>% 
      estimate_density(extend=TRUE) %>% 
      ggplot(aes(x = x, y = y)) +
      geom_area(fill = "orange") +
      theme_classic() +
      # HDI in blue
      #geom_vline(xintercept = ci_hdi$CI_low, color = "royalblue", size = 3) +
      #geom_vline(xintercept = ci_hdi$CI_high, color = "royalblue", size = 3) +
      # Quantile in red
      geom_vline(xintercept = ci_eti$CI_low, color = "red", size = 1) +
      geom_vline(xintercept = ci_eti$CI_high, color = "red", size = 1)
    
  }  
  
    
    
    nd <- dim(x)[2]
    loglh <- x[, nd]
    pairs(x[, 1:(nd-2)], labels=labels, lower.panel=fun.lower,
          upper.panel=fun.upper, diag.panel=fun.diag.kde)
    # was x[,1:(nd-1)]
    
    pm <- ggpairs(
      x[, 1:(nd-2)], columnLabels = labels,
      lower = list(continuous = fun.lower),
      diag  = list(continuous = fun.diag.kde),
      upper = list(continuous = fun.upper)
    )
    pm
    invisible(NULL)
  }