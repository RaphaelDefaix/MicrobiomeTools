#############################################################
## Add statistical annotations to alpha diversity plot
#############################################################

add_alpha_annotations <- function(
    plot,
    annotations,
    text_size = 3
){
  
  plot +
    
    geom_text(
      data = annotations,
      aes(
        x = Day,
        y = y.position,
        label = p.adj.signif
      ),
      inherit.aes = FALSE,
      size = text_size
    )
  
}
