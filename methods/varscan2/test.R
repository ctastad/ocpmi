data(iris)
pdf('testplot.pdf')
plot(iris$Petal.Length, iris$Petal.Width)
dev.off()

ds <- data(iris)
pdf(paste0(ds, ".pdf"))
plot(iris$Petal.Length, iris$Petal.Width)
dev.off()
