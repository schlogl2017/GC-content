library("ggplot2")
data(iris)

g <- ggplot(data=iris, aes(Sepal.Length)) +
     geom_histogram()

ggsave('iris_sepal_length_histogram.png')
