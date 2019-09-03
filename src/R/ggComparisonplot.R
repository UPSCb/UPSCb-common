library(MASS)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
ggComparisonPlot <- function(x,y){
  DF <- data.frame("x" = x, "y" = y)
  dens <- kde2d(x,y)
  gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
  names(gr) <- c("xgr", "ygr", "zgr")
  mod <- loess(zgr~xgr*ygr, data=gr,span = 0.1)
  DF$pointdens <- predict(mod, newdata=data.frame(xgr=x, ygr=y))
  pl1 <- ggplot(DF, aes(x=x,y=y, color=pointdens)) + geom_point() + 
    scale_color_gradientn(colours = rev(brewer.pal(11,"RdYlBu"))) +
    coord_cartesian(xlim = c(-0.1,1.1), ylim = c(-0.1,1.1)) +
    theme(legend.position = 'none', axis.title = element_blank(), plot.margin = unit(c(1,1,1,1),"mm"))
  plxh <- ggplot(DF, aes(x = x)) + geom_histogram(binwidth=0.05) +  
    coord_cartesian(xlim = c(-0.1,1.1)) +
    scale_y_reverse() +
    theme(legend.position = 'none', axis.title = element_blank(), plot.margin = unit(c(1,1,1,1),"mm"))
  plyh <- ggplot(DF, aes(x = y)) + geom_histogram(binwidth=0.05) +  
    coord_cartesian(xlim = c(-0.1,1.1)) +
    coord_flip() +
    theme(legend.position = 'none', axis.title = element_blank(), plot.margin = unit(c(1,1,1,1),"mm"))
  plyb <- ggplot(DF, aes(x = 1, y=y)) + geom_boxplot() + 
    coord_cartesian(ylim = c(-0.1,1.1)) +
    theme(legend.position = 'none', axis.title = element_blank(), plot.margin = unit(c(1,1,1,1),"mm"))
  plxb <- ggplot(DF, aes(x = 1, y=x)) + geom_boxplot() + 
    coord_cartesian(xlim = c(-0.1,1.1)) +
    coord_flip() +
    theme(legend.position = 'none', axis.title = element_blank(), plot.margin = unit(c(1,1,1,1),"mm"))
  layout <- matrix(data = c(0,1,0,2,3,4,0,5,0), nrow = 3, byrow = TRUE)
  plots <- list(plxb,plyb,pl1,plyh,plxh)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout),
                                             widths = c(0.3,1,0.3), 
                                             heights = c(0.3,1,0.3))))
  
  for (i in 1:5) {
    # Get the i,j matrix positions of the regions that contain this subplot
    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
    
    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                    layout.pos.col = matchidx$col))
  }
}