# this function doesn't monitor press perturbations to validate. 
# Instead just simulate and get stable matrices, 
# record outcome (+ve, -ve, neutral) for landward and seaward mangroves
# get weights 

system.sim_press <- function (n.sims, edges, required.groups = c(0), 
                              sampler = community.sampler(edges, required.groups),  
                              perturb) {
  stableout <- list()
  stablews <- list()
  stable <- 0
  unstable <- 0
  
  labels <- node.labels(edges)
  index <- function(name) {
    k <- match(name, labels)
    if (any(is.na(k))) 
      warning("Unknown nodes:", paste(name[is.na(k)], collapse = " "))
    k
  }
  
  k.perturb <- index(names(perturb))
  S.press <- double(length(labels))
  S.press[k.perturb] <- -perturb
  
  while (stable < n.sims) {
    z <- sampler$select(runif(1))
    W <- sampler$community()
    if (!stable.community(W) & !is.null(solve(W, S.press))){
      unstable <- unstable + 1
      next
    } else{
    stable <- stable + 1
    stableout[[stable]] <- data.frame(nsim = stable, var = labels, outcome = solve(W, S.press))
    stablews[[stable]] <- data.frame(nsim = stable, param = sampler$weight.labels, weight = sampler$weights(W))
    }
    }
  
  stableout <- do.call(rbind, stableout)
  stablews <- do.call(rbind, stablews)
  stability <- data.frame(Num_unstable = unstable, Num_stable = stable, Pot_stability =  stable/(unstable+stable))
  list(edges = edges, stability.df = stability, stableoutcome = stableout, stableweights = stablews)
}

# this function can only do one press scenario at a time
# unlike original 'system.simulate' which can accept models that are valid under
# multiple scenarios

system.sim_press_val <- function (n.sims, edges, required.groups = c(0), 
          sampler = community.sampler(edges, required.groups),  
          perturb, monitor, epsilon = 1e-05) {
  
  allout <- list()
  valout <- matrix(0, n.sims, length(node.labels(edges)))
  ws <- matrix(0, n.sims, nrow(edges))
  total <- 0
  stable <- 0
  accepted <- 0
  
  labels <- node.labels(edges)
  index <- function(name) {
    k <- match(name, labels)
    if (any(is.na(k))) 
      warning("Unknown nodes:", paste(name[is.na(k)], collapse = " "))
    k
  }
  
  k.perturb <- index(names(perturb))
  k.monitor <- index(names(monitor))
  S.press <- double(length(labels))
  S.press[k.perturb] <- -perturb
  monitor <- sign(monitor)
  
  while (accepted < n.sims) {
    total <- total + 1
    z <- sampler$select(runif(1))
    W <- sampler$community()
    if (!stable.community(W)) 
      next
    stable <- stable + 1
    allout[[total]] <- tryCatch(solve(W, S.press), error = function(e) NULL)
      s <- tryCatch(solve(W, S.press), error = function(e) NULL)
      val <- !is.null(s) && all(signum(s[k.monitor], epsilon) == monitor)
    if(val == FALSE)
      next
    accepted <- accepted + 1
    valout[accepted, ] <- s
    ws[accepted, ] <- sampler$weights(W)
  }
  allout <- do.call(rbind, allout)
  colnames(allout) <- node.labels(edges)
  colnames(valout) <- node.labels(edges)
  colnames(ws) <- sampler$weight.labels
  list(edges = edges, allout = allout, valout = valout, valweights = ws, total = total, stable = stable, 
       accepted = accepted)
}

# save DiagrammR plot

save_png <- function(plot, path){
  DiagrammeRsvg::export_svg(plot) %>%
    charToRaw() %>%
    rsvg::rsvg() %>%
    png::writePNG(path)
}
