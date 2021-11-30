simSCProfiles_modified <- function(
  object,
  cell.ID.column,
  cell.type.column,
  n.cells,
  suffix.names = "_Simul",
  cell.types = NULL,
  file.backend = NULL,
  name.dataset.backend = NULL,
  compression.level = NULL,
  block.processing = FALSE,
  block.size = 1000,
  chunk.dims = NULL,
  verbose = TRUE
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("Object provided is not of DigitalDLSorter class")
  }
  if (is.null(digitalDLSorteR::zinb.params(object))) {
    stop("'zinb.params' slot is empty. To simulate single-cell profiles, ",
         "DigitalDLSorter object must contain estimated parameters ",
         "for data by estimateZinbwaveParams function")
  }
  if (is.null(digitalDLSorteR::single.cell.real(object))) {
    stop("'single.cell.real' slot is empty. To simulate single-cell ", 
         "profiles, DigitalDLSorter object must contain the original ", 
         "data. See ?loadSCProfiles")
  }
  if (!is.null(digitalDLSorteR::single.cell.simul(object = object))) {
    warning("'single.cell.simul' slot already has a SingleCellExperiment 
            object. Note that it will be overwritten\n", 
            call. = FALSE, immediate. = TRUE)
  }
  if (missing(cell.ID.column) || is.null(cell.ID.column)) {
    stop("'cell.ID.column' argument is needed. Please, see ",
         "?simSCProfiles")
  } else if (missing(cell.type.column) || is.null(cell.type.column)) {
    stop("'cell.type.column' argument is needed. Please, see ?simSCProfiles")
  } else if (missing(n.cells) || is.null(n.cells)) {
    stop("'n.cells' argument is needed. Please, see ?simSCProfiles")
  }
  # check if parameters related to hdf5 file are correct
  if (!is.null(file.backend)) {
    if (!requireNamespace("DelayedArray", quietly = TRUE) || 
        !requireNamespace("HDF5Array", quietly = TRUE)) {
      stop("digitalDLSorteR provides the possibility of using HDF5 files as back-end
         when data are too big to be located in RAM. It uses DelayedArray, 
         HDF5Array and rhdf5 to do it. Please Please install both packages to 
         use this functionality")
    } 
    hdf5Params <- digitalDLSorteR:::.checkHDF5parameters(
      file.backend = file.backend, 
      name.dataset.backend = name.dataset.backend, 
      compression.level = compression.level
    )
    name.dataset.backend <- hdf5Params[[1]]
    compression.level <- hdf5Params[[2]]
  }
  # extract data and model
  list.data <- digitalDLSorteR:::.extractDataFromSCE(
    SCEobject = digitalDLSorteR::single.cell.real(object),
    cell.ID.column = cell.ID.column,
    new.data = FALSE
  )
  zinb.object <- digitalDLSorteR::zinb.params(object)
  # check if cell.ID.column and cell.type.column are correct
  mapply(
    function(x, y) {
      digitalDLSorteR:::.checkColumn(
        metadata = list.data[[2]],
        ID.column = x,
        type.metadata = "cells.metadata",
        arg = y
      )
    },
    c(cell.ID.column, cell.type.column),
    c("'cell.ID.column'", "'cell.type.column'")
  )
  # for cases where estimation was performed in some cell types
  cell.types.used <- unique(list.data[[2]][rownames(zinb.object@model@X), 
                                    cell.type.column])
  # generate metadata for simulated cells
  colnames(list.data[[1]]) <- paste(list.data[[2]][, cell.type.column],
                                    list.data[[2]][, cell.ID.column],
                                    sep = "_") 
  list.data[[2]]$simCellName <- paste(list.data[[2]][, cell.type.column],
                                      list.data[[2]][, cell.ID.column],
                                      sep = "_")
  list.data[[2]]$Simulated <- FALSE
  # cell types in model
  cell.set.names <- NULL
  model.cell.types <- grep(
    pattern = cell.type.column,
    x = colnames(zinb.object@model@X),
    value = T
  )
  cell.type.names <- sub(
    pattern = cell.type.column,
    replacement = "",
    x = model.cell.types
  )
  if (!is.null(cell.types)) {
    if (!all(cell.types %in% cell.types.used)) {
      stop("Cell type(s) provided in 'cell.types' not found in ZINB-WaVE model.",
           "\n  Only cell types that have been used during estimation of ",
           "parameters can be simulated")
    }
    cell.sel <- cell.type.names %in% cell.types
    cell.type.names <- cell.type.names[cell.sel]
    model.cell.types <- model.cell.types[cell.sel]
  }
  names(cell.type.names) <- model.cell.types
  
  names(n.cells) <- model.cell.types
  for (s in model.cell.types) {
    cell.type.name <- cell.type.names[s]
    cell.index <- rownames(zinb.object@model@X)[which(zinb.object@model@X[, s] == 1)]
    nams <- sample(cell.index, size = n.cells[s], replace = T)
    if (is.null(cell.set.names)) {
      cell.set.names <- nams
      names(cell.set.names) <- paste(
        cell.type.name, suffix.names, seq(from = 1, to = n.cells[s]), sep = ""
      )
    } else {
      ns <- names(cell.set.names)
      cell.set.names <- c(cell.set.names, nams)
      names(cell.set.names) <- c(
        ns, paste(
          cell.type.name, suffix.names, 
          seq(from = length(ns) + 1, to = length(ns) + n.cells[s]), sep = ""
        )
      )
    }
  }
  intercept.celltype <- FALSE
  if (!is.null(cell.types)) {
    inter.cell.type <- setdiff(cell.types, cell.type.names)
    if (length(inter.cell.type) != 0) intercept.celltype <- TRUE
  } else {
    inter.cell.type <- setdiff(cell.types.used, cell.type.names)
    intercept.celltype <- TRUE
  }
  if (intercept.celltype) {
    # to get the intercept cell type the rowSum of all FinalCellType columns should be 0
    cell.index <- rownames(zinb.object@model@X)[rowSums(
      zinb.object@model@X[, grep(cell.type.column, 
                                 colnames(zinb.object@model@X), 
                                 value = T), drop = FALSE]
    ) == 0]
    nams <- sample(cell.index, size = n.cells[s], replace = T)
    ns <- names(cell.set.names)
    cell.set.names <- c(cell.set.names, nams)
    names(cell.set.names) <- c(
      ns, paste(inter.cell.type, suffix.names,
                seq(from = length(ns) + 1,
                    to = length(ns) + n.cells[s]),
                sep = "")
    )
  }
  # getting parameters from zinb-wave model
  mu <- zinbwave::getMu(zinb.object@model) # rows are cells
  pi <- zinbwave::getPi(zinb.object@model) # rows are cells
  theta <- zinbwave::getTheta(zinb.object@model) # for genes
  # setting dimensions of simulated matrix
  n <- length(cell.set.names) # as.numeric(nCells)
  J <- zinbwave::nFeatures(zinb.object@model)
  if (verbose) {
    # message about parameters
    message("=== Getting parameters from model:")
    message("    - mu: ", paste(dim(mu), collapse = ", "))
    message("    - pi: ", paste(dim(pi), collapse = ", "))
    message("    - Theta: ", length(theta), "\n")
    # message about selected cell types
    message(paste0("=== Selected cell type(s) from ZINB-WaVE model (", 
                   length(c(cell.type.names, inter.cell.type)), 
                   " cell type(s)):"))
    message(paste0("    - ", c(cell.type.names, inter.cell.type), 
                   collapse = "\n"), "\n")
    # messages about simulated matrix
    message("=== Simulated matrix dimensions:")
    message("    - n (cells): ", n)
    message("    - J (genes): ", J)
    message("    - i (# entries): ", n * J)
  }
  if (block.processing && is.null(file.backend)) {
    stop("'block.processing' is only compatible with the use of HDF5 files ", 
         "as back-end ('file.backend' argument)")
  } else if (block.processing && !is.null(file.backend)) {
    if (!requireNamespace("DelayedArray", quietly = TRUE) || 
        !requireNamespace("HDF5Array", quietly = TRUE)) {
      stop("digitalDLSorteR provides the possibility of using HDF5 files as back-end
         when data are too big to be located in RAM. It uses DelayedArray, 
         HDF5Array and rhdf5 to do it. Please install both packages to 
         use this functionality")
    } 
    if (n < block.size) {
      block.size <- n
      warning("The number of simulated cells is less than 'block.size'. ",
              "Only one block will be performed", 
              call. = FALSE, immediate. = TRUE)
    }
    if (verbose) 
      message("\n=== Simulating and writing new single-cell profiles by blocks")
    if (is.null(chunk.dims)){
      chunk.dims <- c(J, 1)
    } else {
      if (any(chunk.dims > c(J, n))) {
        warning("'chunk.dims' must be equal to or less than dimension of ", 
                "data. Setting default value", call. = FALSE, immediate. = TRUE)
        chunk.dims <- c(J, 1)
      }
    }
    if (!file.exists(file.backend)) rhdf5::h5createFile(file = file.backend) 
    rhdf5::h5createDataset(
      file = file.backend, dataset = name.dataset.backend, 
      dims = c(J, block.size), maxdims = c(J, n), 
      chunk = chunk.dims, storage.mode = "integer"
    )
    # pointers for hdf5 file and zinb-wave parameters
    r.i <- 0
    r.j <- 0
    # iteration over cells 
    for (iter in seq(ceiling(n / block.size))) {
      if ((block.size * iter) - n > 0) { # && dif < block.size
        dif <- (block.size * iter) - n
        block.size <- block.size - dif
      } else {
        dif <- block.size
      }
      sub.i <- seq(from = r.i + 1, to = r.i + block.size)
      sub.j <- seq(from = r.j + 1, to = r.j + block.size * J)
      r.i <- r.i + block.size
      r.j <- r.j + block.size
      datanb <- rnbinom(length(sub.j), mu = mu[cell.set.names[sub.i], ], 
                        size = rep(theta[1], length(sub.j))) # ceiling(sub.i/block.size)
      data.nb <- matrix(datanb, nrow = block.size)
      datado <- rbinom(length(sub.j), size = 1, prob = pi[cell.set.names[sub.i], ])
      data.dropout <- matrix(datado, nrow = block.size)
      sim.counts <- t(data.nb * (1 - data.dropout))
      if (iter == 1) {
        rhdf5::h5write(
          obj = sim.counts, file = file.backend, 
          name = name.dataset.backend, level = compression.level
        )
      } else {
        # check number of cells in the next loop
        rhdf5::h5set_extent(
          file = file.backend, dataset = name.dataset.backend, dims = c(J, n)
        )
        rhdf5::h5write(
          obj = sim.counts, 
          file = file.backend, name = name.dataset.backend, 
          index = list(seq(J), seq((dif * (iter - 1)) + 1, 
                                   (dif * (iter - 1)) + ncol(sim.counts))),
          level = compression.level
        )
      }
      if (verbose) message("    - Block ", iter, " written")
    }
    rhdf5::H5close()
    # HDF5Array object for SingleCellExperiment class
    sim.counts <- HDF5Array::HDF5Array(
      file = file.backend, 
      name = name.dataset.backend
    )
    dimnames(sim.counts) <- list(rownames(zinb.object@model@V), 
                                 names(cell.set.names))
  } else if (!block.processing) {
    i <- seq(n * J)
    mu <- mu[cell.set.names, ]
    pi <- pi[cell.set.names, ]
    datanb <- rnbinom(n * J, mu = mu[i], size = theta[ceiling(i/n)])
    data.nb <- matrix(datanb, nrow = n)
    
    datado <- rbinom(length(i), size = 1, prob = pi[i])
    data.dropout <- matrix(datado, nrow = n)
    
    sim.counts <- t(data.nb * (1 - data.dropout))
    sim.counts <- Matrix::Matrix(
      sim.counts, dimnames = list(
        rownames = rownames(zinb.object@model@V),
        colnames = names(cell.set.names))
    )  
  }
  sim.cells.metadata <- list.data[[2]][cell.set.names, ]
  sim.cells.metadata$simCellName <- NULL # it is not used
  sim.cells.metadata[, cell.ID.column] <- names(cell.set.names)
  rownames(sim.cells.metadata) <- names(cell.set.names)
  sim.cells.metadata$Simulated <- TRUE
  if (any(colnames(sim.cells.metadata) == "suffix"))
    warning("\n'suffix' column in cells metadata is going to be overwritten")
  sim.cells.metadata$suffix <- suffix.names
  
  sim.sce <- digitalDLSorteR:::.createSCEObject(
    counts = sim.counts,
    cells.metadata = sim.cells.metadata,
    genes.metadata = list.data[[3]][rownames(sim.counts), ],
    file.backend = file.backend,
    name.dataset.backend = name.dataset.backend,
    compression.level = compression.level,
    chunk.dims = chunk.dims,
    block.processing = block.processing,
    verbose = verbose
  )
  digitalDLSorteR::single.cell.simul(object) <- sim.sce
  if (verbose) message("\nDONE")
  return(object)
}