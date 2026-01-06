library(ggplot2)
library(ggtext)
library(glue)

df <- data.frame(
  x = c("A", "B"),
  y = c(1, 2)
)

colors <- c("A" = "red", "B" = "blue")

p <- ggplot(df, aes(x, y)) +
  geom_point() +
  scale_x_discrete(labels = function(x) {
    glue("<span style='color:{colors[x]}'>{x}</span>")
  }) +
  theme(
    axis.text.x = element_markdown()
  )

ggsave("debug_ggtext.png", p)
print(p)
