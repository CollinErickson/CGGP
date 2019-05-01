# Show possible blocks

uo <- SG$uo[1:SG$uoCOUNT,]
po <- SG$po[1:SG$poCOUNT,]
bdf <- rbind(cbind(data.frame(uo), Indexes="Selected"),
             cbind(data.frame(po), Indexes="Valid"))

bdf$xmin <- bdf$X1-1
bdf$xmax <- bdf$X1
bdf$ymin <- bdf$X2-1
bdf$ymax <- bdf$X2
ggplot() + geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Indexes), bdf, color="black") + scale_fill_manual(values=c('yellow', 'cyan')) + coord_fixed()


# Plot ancestors
plot_ancestors <- function(..., SG, uo=SG$uo[1:SG$uoCOUNT,]) {
  
  bdf <- rbind(cbind(data.frame(X1=uo[,1], X2=uo[,2]), Indexes="Selected"))
  
  bdf$xmin <- bdf$X1-1
  bdf$xmax <- bdf$X1
  bdf$ymin <- bdf$X2-1
  bdf$ymax <- bdf$X2
  arrowdf <- NULL
  for (irow in 1:nrow(uo)) {
    if (uo[irow,1] > 1) {
      arrowdf <- rbind(arrowdf,
                       data.frame(xhead=uo[irow,1]-.75, xtail=uo[irow,1]-1.25, yhead=uo[irow,2]-.5, ytail=uo[irow,2]-.5))
    }
    if (uo[irow,2] > 1) {
      arrowdf <- rbind(arrowdf,
                       data.frame(xhead=uo[irow,1]-.5, xtail=uo[irow,1]-.5, yhead=uo[irow,2]-.75, ytail=uo[irow,2]-1.25)
      )
    }
  }
  ggplot() + geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), bdf, fill='yellow', color="black") + scale_fill_manual(values=c('yellow', 'cyan')) + coord_fixed() +
    geom_segment(aes(x=xhead, y=yhead, xend=xtail, yend=ytail), arrowdf, arrow=arrow(length=unit(.03, "npc")), size=2) + xlab(expression(i[1])) + ylab(expression(i[2])) + ggtitle("Ancestor blocks")
}
eg4 <- expand.grid(1:4,1:4)
plot_ancestors(uo=eg4[c(1,2,3,4,5,6,9, 7),])
plot_ancestors(uo=eg4[c(1,2,3,4,8,11,5,9),])
