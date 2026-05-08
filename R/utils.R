# z_ig or "responsibility" or
e_step <- function(density){
  total_density <- rowSums(density)
  responsibilities <- density / total_density
  list(responsibilities = responsibilities, total_density = total_density)
}
