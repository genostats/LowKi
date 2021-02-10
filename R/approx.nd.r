
# x = point où on veut l'approx (longueur n)
# L = liste de vecteurs (longueur n) donnant les coordonnées de la grille où les valeurs sont connues
# F = array des valeurs  
# rule = comment gérer valeurs hors grille (1 = NA, 2 = constant, 3 = prolongement linéaire...)
#        même longueur que x
# donne l'approx en ( L[[1]][ I[1] ], L[[2]][ I[2] ], ... , x) 
approx.nd <- function(x, L, F, rule.left = 1, rule.right = rule.left) {
  n <- length(x)
  if(length(L) != n) stop("Dimensions mismatch")
  if(n > 1) {
    if(length(rule.left) == 1)
      rule.left <- rep(rule.left, n)
    if(length(rule.right) == 1)
      rule.right <- rep(rule.right, n)
  }
  aprx(x, L, F, rule.left, rule.right, integer(n), 0)
}