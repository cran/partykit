
nodeprune.party <- function(x, ids, ...) {
  
  ### map names to nodeids
  if (!is.numeric(ids))
    ids <- match(ids, names(x))
  stopifnot(ids %in% nodeids(x))
  
  ### compute indices path to each node
  ### to be pruned off
  idxs <- lapply(ids, .get_path, obj = node_party(x))
  
  ### [[.party is NOT [[.list
  cls <- class(x)
  x <- unclass(x)
  ni <- which(names(x) == "node")
  
  for (i in 1:length(idxs)) {
    
    idx <- c(ni, idxs[[i]])
    ### check if we already pruned-off this node
    tmp <- try(x[[idx]], silent = TRUE)
    if (inherits(tmp, "try-error"))
      next()
    
    ### node ids of off-pruned daugther nodes
    idrm <- nodeids(x[[idx]])[-1]
    
    ### prune node by introducing a "new" terminal node
    x[[idx]] <- partynode(id = id_node(x[[idx]]),
                          info = info_node(x[[idx]]))
    
    ### constparty only: make sure the node ids in
    ### fitted are corrected
    if (length(idrm) > 0) {
      if(!is.null(x$fitted) && 
         "(fitted)" %in% names(x$fitted)) {
        j <- x$fitted[["(fitted)"]] %in% idrm
        x$fitted[["(fitted)"]][j] <- ids[i]
      }
    }
  }
  
  ### reindex to 1:max(nodeid)
  class(x) <- cls
  oldids <- nodeids(x)
  newids <- 1:length(nodeids(x))
  nodeids(x) <- newids ### this takes also care of $fitted!

  return(x)
}
