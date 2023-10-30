# below are a set of functions either taken from, or adapted from, the QPress R package: https://github.com/SWotherspoon/QPress
# devtools::install_github("SWotherspoon/QPress",ref="Constrain")

# simulate stable matrices, no perturbations

system.sim <- function (n.sims, constrainedigraph, required.groups = c(0), from, to, class, spatial, arid, prob, cdev,
                              sampler = community.sampler_con1(constrainedigraph, required.groups, from, to, class, spatial = 'Y', arid, cdev)) {
  
  edges1 <- constrainedigraph$edges
  labels <- node.labels(edges1)
  stablews_array <- array(NA, dim = c(length(labels), length(labels), n.sims))
  stable <- 0
  unstable <- 0

  while (stable < n.sims) {
    z <- sampler$select(0.5)
    z2 <- sampler$select2(prob)
    z3 <- sampler$select3(cdev[1], cdev[2])
    W <- sampler$community()
    if (!stable.community(W)){
      unstable <- unstable + 1
      next
    } else{
      stable <- stable + 1
      stablews_array[,,stable] <- W
    }
  }

  list(edges = edges1, stableweights = stablews_array)
}

# solve matrices given pressures

solver <- function(x, chosen_model, perturb){ # x is a biomodel matrix
  labels <- node.labels(chosen_model)
  index <- function(name) {
    k <- match(name, labels)
    if (any(is.na(k))) 
      warning("Unknown nodes:", paste(name[is.na(k)], collapse = " "))
    k}
  k.perturb <- index(c(perturb)) # nodes to perturb
  S.press <- double(length(labels))
  S.press[k.perturb] <- -rep(1, length(k.perturb))
  if(!is.null(solve(x, S.press))){
    out <- data.frame(node = labels, outcome = solve(x, S.press)) %>% 
      filter(node %in% c('LandwardMang', 'SeawardMang')) %>% 
      mutate(outcome = ifelse(outcome >= 0, 1, -1))
    return(out)
  }
}

# this function simulates matrices and gets stable ones and handles perturbations

system.sim_press <- function (n.sims, constrainedigraph, required.groups = c(0), from, to, class, arid, prob, cdev,
                              sampler = community.sampler_con2(constrainedigraph, required.groups, from, to, class, perturb, spatial, arid, cdev),  
                              perturb, spatial) {
  
  edges1 <- constrainedigraph$edges
  labels <- node.labels(edges1)
  index <- function(name) {
    k <- match(name, labels)
    if (any(is.na(k))) 
      warning("Unknown nodes:", paste(name[is.na(k)], collapse = " "))
    k
  }
  
  stableout <- list()
  stablews_array <- array(NA, dim = c(length(labels), length(labels), n.sims))
  stable <- 0
  unstable <- 0
  
  k.perturb <- index(names(perturb))
  S.press <- double(length(labels))
  S.press[k.perturb] <- -perturb
  
  while (stable < n.sims) {
    z <- sampler$select(0.5)
    z2 <- sampler$select2(prob)
    z3 <- sampler$select3(cdev[1], cdev[2])
    W <- sampler$community()
    if (!stable.community(W) & !is.null(solve(W, S.press))){
      unstable <- unstable + 1
      next
    } else{
      stable <- stable + 1
      stableout[[stable]] <- data.frame(nsim = stable, var = labels, outcome = solve(W, S.press))
      stablews_array[,,stable] <- W
    }
  }
  
  stableout <- do.call(rbind, stableout)
  stability <- data.frame(Num_unstable = unstable, Num_stable = stable, Potential_stability =  stable/(unstable+stable))
  list(edges = edges1, stability.df = stability, stableoutcome = stableout, stableweights = stablews_array)
}


# save DiagrammR plot

save_png <- function(plot, path){
  DiagrammeRsvg::export_svg(plot) %>%
    charToRaw() %>%
    rsvg::rsvg() %>%
    png::writePNG(path)
}

# Construct bounding sets used to order weights to meet the edge weight constraints.
# These bounds sets are used to order a set of random edge weights
# so that they meet the imposed constraints.

bound.sets <- function(constrained) {
  a <- constrained$a
  b <- constrained$b
  
  ## Find all weights bounded above by the root weight r
  bounded <- function(r) {
    v <- unique(a[b %in% r])
    repeat {
      va <- setdiff(a[b %in% v],v)
      if(length(va)==0) break;
      v <- c(va,v)
    }
    v
  }
  
  bounds <- list()
  n <- 0
  
  ## Unvisited weights
  us <- b
  ## Unbounded weights
  vs <- setdiff(b,a)
  while(length(vs) > 0) {
    for(v in vs) {
      ## Find all weights bounded by v
      bs <- bounded(v)
      if(v %in% bs) warning("Cyclic constraints found")
      bounds[[n <- n+1]] <- c(v,bs)
    }
    ## Unvisited weights
    us <- setdiff(us,vs)
    ## Weights bounded only by visited weights
    vs <- setdiff(us,a[b %in% us])
  }
  bounds
}

# Order a set of absolute edge weights to meet imposed edge weights
# constraints.

constraint.order <- function(w,bounds) {
  for(bs in bounds) {
    k <- which.max(w[bs])
    if(k!=1) w[bs[c(1,k)]] <- w[bs[c(k,1)]]
  }
  w
}

# community sampler that allows weights to be constrained for certain edges
# and allows relative strengths of edges to be constrained
# based on 'community.ordering.sampler' function from Wotherspoon
# does not handle perturbations
community.sampler_con1 <- function (constrainedigraph, required.groups = c(0), from, to, class, spatial, arid, cdev)# from, to, and class arguments for constraining edges as 'High', "Med', or 'Low', just a vector
{
  edges <- constrainedigraph$edges
  if (length(from) > 0){ # here add new column to edges with high, medium or low classification
    constrain <- data.frame(From = from, To = to, Class = class)
    if(spatial == 'N'){ # in non-spatial model, only if coastal development is being perturbed do we constrain SeaLevelRise -> LandwardMang edge
      if(length(which(names(perturb) == 'CoastalDev')) != 0){ 
        edges$Class <- dplyr::left_join(edges, constrain)$Class
      }else{
        edges$Class <- dplyr::left_join(edges, constrain)$Class
        edges <- mutate(edges, Class = ifelse(To == 'LandwardMang', 'NA', Class))
      }
    }else{ # in spatial model, always constrain SeaLevelRise -> LandwardMang edge according to pop. density, i.e. coastal squeeze
      edges$Class <- dplyr::left_join(edges, constrain)$Class
    }
  }
  n.nodes <- length(node.labels(edges))
  weight.labels <- edge.labels(edges)
  n.edges <- nrow(edges)
  W <- matrix(0, n.nodes, n.nodes)
  lower <- ifelse(edges$Type == "U" | edges$Type == "N", -1L, 
                  0L)
  upper <- ifelse(edges$Type == "U" | edges$Type == "P" , 1L, 
                  0L)
  ## set up constraints for high, medium, low
  lower[which(edges$Class == 'H' & edges$Type == 'P')] <- 0.66667
  upper[which(edges$Class == 'H' & edges$Type == 'N')] <- -0.66667
  lower[which(edges$Class == 'M' & edges$Type == 'P')] <- 0.33334
  upper[which(edges$Class == 'M' & edges$Type == 'P')] <- 0.66666
  lower[which(edges$Class == 'M' & edges$Type == 'N')] <- -0.66666
  upper[which(edges$Class == 'M' & edges$Type == 'N')] <- -0.33334
  upper[which(edges$Class == 'L' & edges$Type == 'P')] <- 0.33333
  lower[which(edges$Class == 'L' & edges$Type == 'N')] <- -0.33333
  k.edges <- as.vector(unclass(edges$To) + (unclass(edges$From) - 
                                              1) * n.nodes)
  uncertain <- which(!(edges$Group %in% required.groups))
  expand <- match(edges$Pair[uncertain], unique(edges$Pair[uncertain]))
  n.omit <- max(0, expand)
  bounds <- bound.sets(constrainedigraph)
  zs <- rep(1, length(uncertain))
  zss <- 1
  zsss <- c(1,1)
  community <- if (n.omit > 0) {
    function() {
      r <- runif(n.edges, lower, upper)
      r <- sign(r) * constraint.order(abs(r), bounds)
      r[uncertain] <- r[uncertain] * zs
      r[3] <- r[3] * zss # here making land propagule to land mang link uncertain
      r[c(9,10)] <- r[c(9,10)] * zsss # here making coastal development to landward & seaward mangrove uncertain
      W[k.edges] <- r
      W
    }
  }
  else {
    function() {
      r <- runif(n.edges, lower, upper)
      r <- sign(r) * constraint.order(abs(r), bounds)
      W[k.edges] <- r
      W
    }
  }
  select <- if (n.omit > 0) {
    function(p) {
      zs <<- rbinom(n.omit, 1, p)[expand]
    }
  }
  else {
    function(p = 0) {
      zs
    }
  }
  select2 <- if (arid == 'Y') {
    function(p2) {
      zss <<- rbinom(1, 1, p2)
    }
  }  
  else {
    function(p2 = 0) {
      zss
    }
  }
  select3 <- if (n.omit > 0) {
    function(p3, p4){
      zsss <<- c(rbinom(1, 1, runif(1,p3,p4)), if(p3 == 0 & p4 == 0){0}else{rbinom(1,1,runif(1,0,0.33))}) # coastal development to seaward mangrove is always uncertain, i.e. has a low probability of 0 to 0.33
    }
  }else{
    function(p3 = 0, p4 = 0){
      zsss
    }
  }
  list(community = community, select = select, select2 = select2, select3 = select3, weights = function(W) {
    W[k.edges]
  }, 
  weight.labels = weight.labels, uncertain.labels = weight.labels[uncertain])
}

# this one does handle perturbations
community.sampler_con2 <- function (constrainedigraph, required.groups = c(0), from, to, class, perturb, spatial, arid, cdev) # from, to, and class arguments for constraining edges as 'High', "Med', or 'Low', just a vector
{
  edges <- constrainedigraph$edges
  if (length(from) > 0){ # here add new column to edges with high, medium or low classification
    constrain <- data.frame(From = from, To = to, Class = class)
    if(spatial == 'N'){ # in non-spatial model, only if coastal development is being perturbed do we constrain SeaLevelRise -> LandwardMang edge
      if(length(which(names(perturb) == 'CoastalDev')) != 0){ 
        edges$Class <- dplyr::left_join(edges, constrain)$Class
      }else{
        edges$Class <- dplyr::left_join(edges, constrain)$Class
        edges <- mutate(edges, Class = ifelse(To == 'LandwardMang', 'NA', Class))
      }
    }else{ # in spatial model, always constrain SeaLevelRise -> LandwardMang edge according to pop. density, i.e. coastal development
      edges$Class <- dplyr::left_join(edges, constrain)$Class
    }
  }
  n.nodes <- length(node.labels(edges))
  weight.labels <- edge.labels(edges)
  n.edges <- nrow(edges)
  W <- matrix(0, n.nodes, n.nodes)
  lower <- ifelse(edges$Type == "U" | edges$Type == "N", -1L, 
                  0L)
  upper <- ifelse(edges$Type == "U" | edges$Type == "P" , 1L, 
                  0L)
  ## set up constraints for high, medium, low
  lower[which(edges$Class == 'H' & edges$Type == 'P')] <- 0.66667
  upper[which(edges$Class == 'H' & edges$Type == 'N')] <- -0.66667
  lower[which(edges$Class == 'M' & edges$Type == 'P')] <- 0.33334
  upper[which(edges$Class == 'M' & edges$Type == 'P')] <- 0.66666
  lower[which(edges$Class == 'M' & edges$Type == 'N')] <- -0.66666
  upper[which(edges$Class == 'M' & edges$Type == 'N')] <- -0.33334
  upper[which(edges$Class == 'L' & edges$Type == 'P')] <- 0.33333
  lower[which(edges$Class == 'L' & edges$Type == 'N')] <- -0.33333
  k.edges <- as.vector(unclass(edges$To) + (unclass(edges$From) - 
                                              1) * n.nodes)
  uncertain <- which(!(edges$Group %in% required.groups))
  expand <- match(edges$Pair[uncertain], unique(edges$Pair[uncertain]))
  n.omit <- max(0, expand)
  bounds <- bound.sets(constrainedigraph)
  zs <- rep(1, length(uncertain))
  zss <- 1
  zsss <- c(1,1)
  community <- if (n.omit > 0) {
    function() {
      r <- runif(n.edges, lower, upper)
      r <- sign(r) * constraint.order(abs(r), bounds)
      r[uncertain] <- r[uncertain] * zs
      r[3] <- r[3] * zss # here making land propagule to land mang link uncertain
      r[c(9,10)] <- r[c(9,10)] * zsss # here making coastal development to landward & seaward mangrove uncertain
      W[k.edges] <- r
      W
    }
  }
  else {
    function() {
      r <- runif(n.edges, lower, upper)
      r <- sign(r) * constraint.order(abs(r), bounds)
      W[k.edges] <- r
      W
    }
  }
  select <- if (n.omit > 0) {
    function(p) {
      zs <<- rbinom(n.omit, 1, p)[expand]
    }
  }
  else {
    function(p = 0) {
      zs
    }
  }
  select2 <- if (arid == 'Y') {
    function(p2) {
      zss <<- rbinom(1, 1, p2)
    }
  }  
  else {
    function(p2 = 0) {
      zss
    }
  }
  select3 <- if (n.omit > 0) {
    function(p3, p4){
    zsss <<- c(rbinom(1, 1, runif(1,p3,p4)), if(p3 == 0 & p4 == 0){0}else{rbinom(1,1,runif(1,0,0.33))}) # coastal development to seaward mangrove is always uncertain, i.e. has a low probability of 0 to 0.33
    }
    }else{
      function(p3 = 0, p4 = 0){
    zsss
      }
  }
  weights <- function(W) {
    W[k.edges]
  }
  list(community = community, select = select, select2 = select2, select3 = select3, weights = weights, 
       weight.labels = weight.labels, uncertain.labels = weight.labels[uncertain])
}


