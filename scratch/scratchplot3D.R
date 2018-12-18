skinny <- NULL
heatmatrix <- matrix(NaN, SG$d, SG$d)
for (i in 1:SG$d) {
  heatmatrix[i,i] <- max(SG$designindex[,i])
  skinny <- rbind(skinny, c(i, i, max(SG$designindex[,i])))
}
for (i in 1:(SG$d-1)) {
  for (j in (i+1):SG$d) {
    heatmatrix[i,j] <- heatmatrix[j,i] <- max(apply(SG$designindex[,c(i,j)], 1, min))
    skinny <- rbind(skinny,
                    c(i, j, max(apply(SG$designindex[,c(i,j)], 1, min))),
                    c(j, i, max(apply(SG$designindex[,c(i,j)], 1, min)))
    )
  }
}
# Useless, not what I wanted
plot3D::hist3D(z = heatmatrix, scale = FALSE, expand = 0.01, bty = "g", phi = 20,
       col = "#0072B2", border = "black", shade = 0.2, ltheta = 90,
       space = 0.3, ticktype = "detailed", d = 2)



m2 <- matrix(0, max(SG$designindex), max(SG$designindex))
i <- 1
j <- 2
m2ij <- SG$designindex[,c(i,j)]
dfij <- data.frame(m2ij)
pij <- plyr::ddply(dfij, c('X1','X2'), function(x) nrow(x))
# heatmatrix <- matrix(0, max(SG$designindex), max(SG$designindex))
for (k in 1:nrow(pij)) {
  m2[pij[k,1], pij[k,2]] <- pij[k,3]
}

plot3D::hist3D(x=1:max(SG$designindex), y=1:max(SG$designindex),
               z = m2, scale = FALSE, expand = 0.01, bty = "g", phi = 20,
               col = "#0072B2", border = "black", shade = 0.2, ltheta = 90,
               space = 0.3, ticktype = "detailed", d = 2)
# Can't get it to leave off zeros. Try using colvar, or col=ifelse(m2>0,'blue','white)
