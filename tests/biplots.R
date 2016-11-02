biplot.clr <- function(x, ...)
{

  #Calcula biplot composicional

  #Con el argumento groups se indica la variable que contiene las etiquetas
  #de grupo.

  biplot2 <-
    function(x, y, groups=NULL, cpv=NULL, var.axes = TRUE, col, cex = rep(par("cex"), 2),
             xlabs = NULL, ylabs = NULL, expand=1, xlim = NULL, ylim = NULL,
             arrow.len = 0.1,
             main = NULL, sub = NULL, xlab = "", ylab = "", ...)
    {

      #Vectores de colores y caracteres

      if(is.null(groups)){
        color <- "dodgerblue1"
        char <- 1
      }
      else {
        cod.vector <- rep(0,length(groups))
        lev <- levels(groups)
        cod.label <- rep(0,length(lev))
        for(i in 1:length(lev)){
          cod.vector[which(groups==as.character(lev[i]))] <- i+14 #n? color y char
          cod.label[i] <- i+14
        }
        color <- cod.vector-13
        char <- cod.vector
      }

      n <- nrow(x)
      p <- nrow(y)
      if(missing(xlabs)) {
        xlabs <- dimnames(x)[[1L]]
        if(is.null(xlabs)) xlabs <- 1L:n
      }
      xlabs <- as.character(xlabs)
      dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
      if(missing(ylabs)) {
        ylabs <- dimnames(y)[[1L]]
        if(is.null(ylabs)) ylabs <- paste("Var", 1L:p)
      }
      ylabs <- as.character(ylabs)
      dimnames(y) <- list(ylabs, dimnames(y)[[2L]])

      if(length(cex) == 1L) cex <- c(cex, cex)
      if(missing(col)) {
        col <- par("col")
        if (!is.numeric(col)) col <- match(col, palette(), nomatch=1L)
        col <- c(col, col + 1L)
      }
      else if(length(col) == 1L) col <- c(col, col)


      #En la siguiente funci?n es donde se especifica todo lo relacionado con colores, tama?os
      #s?mbolos, etc.


      unsigned.range <- function(x)
        c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
      rangx1 <- unsigned.range(x[, 1L])
      rangx2 <- unsigned.range(x[, 2L])
      rangy1 <- unsigned.range(y[, 1L])
      rangy2 <- unsigned.range(y[, 2L])

      if(missing(xlim) && missing(ylim))
        xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
      else if(missing(xlim)) xlim <- rangx1
      else if(missing(ylim)) ylim <- rangx2
      ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
      on.exit(par(op))
      op <- par(pty = "s")
      if(!is.null(main))
        op <- c(op, par(mar = par("mar")+c(0,0,1,0)))
      plot(x, type = "p", xlim = xlim, ylim = ylim,
           col = color,pch=char,cex=0.75, #color, s?mbolo y tama?os de los puntos
           xaxt = "n",
           yaxt = "n",  #Quita ejes
           xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
      #text(x[1:90,], xlabs, cex = cex[1L], col = col[1L], ...)
      par(new = TRUE)
      plot(y, axes = FALSE, type = "n", xlim = xlim*ratio, ylim = ylim*ratio,
           xlab = "", ylab = "", col = col[1L], ...)
      #    axis(3, col = col[2L], ...) #Segundo Eje
      #    axis(4, col = col[2L], ...)
      #    box(col = col[1L])
      text(y*0.95, labels=ylabs, cex = 0.8, col = "black", ...)
      if(var.axes)
        arrows(0, 0, y[,1L] * 0.8, y[,2L] * 0.8,
               col = "black", #Arrow color
               length=0.07)   #Tama?o punta flecha

      if(!is.null(groups))  #Leyenda
        legend("bottomleft", cex=0.85,bty="n",leg=lev, col=cod.label-13,
               pch=cod.label)

      # Añadir Prop. Var. Explicada

      #legend("bottomright",cex=0.85,bty="n",
      #leg=paste("Prop. Var. Explained: ",round(cpv[2],2)))

      box(col="black") #Quitar la caja (ejecutarlo varias veces para dejar borde blanco)

      invisible()
    }


  # Transformaci?n clr y PCA

  x<-log(x)-rowMeans(log(x))
  x<-princomp(x)

  #Calcula prop. var. explicada

  cpv <- cumsum(x$sdev^2 / sum(x$sdev^2))

  choices <- 1L:2
  pc.biplot <- FALSE
  scale <- 1
  if(length(choices) != 2) stop("length of choices must be 2")
  if(!length(scores <- x$scores))
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))),domain = NA)
  lam <- x$sdev[choices]
  if(is.null(n <- x$n.obs)) n <- 1
  lam <- lam * sqrt(n)
  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) lam <- lam^scale else lam <- 1
  if(pc.biplot) lam <- lam / sqrt(n)
  biplot2(t(t(scores[, choices]) / lam), t(t(x$loadings[, choices]) * lam),
          cpv=cpv, ...)
  invisible()
}

biplot <- function(x, ...)
{

  # Computes PCA biplot with an option to distinguise groups (argument groups)

  biplot2 <-
    function(x, y, groups=NULL, cpv=NULL, var.axes = TRUE, col, cex = rep(par("cex"), 2),
             xlabs = NULL, ylabs = NULL, expand=1, xlim = NULL, ylim = NULL,
             arrow.len = 0.1,
             main = NULL, sub = NULL, xlab = "", ylab = "", ...)
    {

      # Colours and character vectors

      if(is.null(groups)){
        color <- "blue"
        char <- 1
      }
      else {
        cod.vector <- rep(0,length(groups))
        lev <- levels(groups)
        cod.label <- rep(0,length(lev))
        for(i in 1:length(lev)){
          cod.vector[which(groups==as.character(lev[i]))] <- i #no. colour and char
          cod.label[i] <- i
        }
        color <- cod.vector
        char <- cod.vector
      }

      n <- nrow(x)
      p <- nrow(y)
      if(missing(xlabs)) {
        xlabs <- dimnames(x)[[1L]]
        if(is.null(xlabs)) xlabs <- 1L:n
      }
      xlabs <- as.character(xlabs)
      dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
      if(missing(ylabs)) {
        ylabs <- dimnames(y)[[1L]]
        if(is.null(ylabs)) ylabs <- paste("Var", 1L:p)
      }
      ylabs <- as.character(ylabs)
      dimnames(y) <- list(ylabs, dimnames(y)[[2L]])

      if(length(cex) == 1L) cex <- c(cex, cex)
      if(missing(col)) {
        col <- par("col")
        if (!is.numeric(col)) col <- match(col, palette(), nomatch=1L)
        col <- c(col, col + 1L)
      }
      else if(length(col) == 1L) col <- c(col, col)

      # The next function manages all related to colours, sizes, symbols, etc.

      unsigned.range <- function(x)
        c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
      rangx1 <- unsigned.range(x[, 1L])
      rangx2 <- unsigned.range(x[, 2L])
      rangy1 <- unsigned.range(y[, 1L])
      rangy2 <- unsigned.range(y[, 2L])

      if(missing(xlim) && missing(ylim))
        xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
      else if(missing(xlim)) xlim <- rangx1
      else if(missing(ylim)) ylim <- rangx2
      ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
      on.exit(par(op))
      op <- par(pty = "s")
      if(!is.null(main))
        op <- c(op, par(mar = par("mar")+c(0,0,1,0)))
      plot(x, type = "p", xlim = xlim, ylim = ylim,
           col = color,pch=char,cex=0.75, # point colours, sizes and characters
           xaxt = "n",
           yaxt = "n",  # remove axes
           xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
      #text(x[1:90,], xlabs, cex = cex[1L], col = col[1L], ...)
      par(new = TRUE)
      plot(y, axes = FALSE, type = "n", xlim = xlim*ratio, ylim = ylim*ratio,
           xlab = "", ylab = "", col = col[1L], ...)
      #    axis(3, col = col[2L], ...) # second axis
      #    axis(4, col = col[2L], ...)
      #    box(col = col[1L])
      text(y, labels=ylabs, cex = 0.9, col = "black", ...)
      if(var.axes)
        arrows(0, 0, y[,1L] * 0.8, y[,2L] * 0.8,
               col = "black", # Arrow colour
               length=0.07)   # Arrow end size

      if(!is.null(groups))  # Legend
        legend("bottomleft", cex=0.85,bty="n",leg=lev, col=cod.label,
               pch=cod.label)

      # Add prop. var. explained

      legend("topleft",cex=0.85,bty="n",
             leg=paste("Prop. Var. Explained: ",round(cpv[2],2)))

      box(col="white") #Quitar la caja

      invisible()
    }

  # PCA

  x<-princomp(x)

  # Compute prop. var explained

  cpv <- cumsum(x$sdev^2 / sum(x$sdev^2))

  choices <- 1L:2
  pc.biplot <- FALSE
  scale <- 1
  if(length(choices) != 2) stop("length of choices must be 2")
  if(!length(scores <- x$scores))
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))),domain = NA)
  lam <- x$sdev[choices]
  if(is.null(n <- x$n.obs)) n <- 1
  lam <- lam * sqrt(n)
  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) lam <- lam^scale else lam <- 1
  if(pc.biplot) lam <- lam / sqrt(n)
  biplot2(t(t(scores[, choices]) / lam), t(t(x$loadings[, choices]) * lam),
          cpv=cpv, ...)
  invisible()

}
