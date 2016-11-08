#'
#' ilr coordinates
#'
#' Calculate the ilr coordinates with respect basis given by function
#' ilr_basis
#'
#' @param X compositional dataset. Either a matrix, a data.frame or a vector
#' @return coordinates with respect basis ilr_basis
#' @examples
#' X = rnormalmultinomial(100, 100, rep(0,4), diag(4))
#' nm_fit(X, verbose = T)
#' @export
#' @export
ilr_coordinates = function(X){
  class_type = class(X)
  is_vector = is.vector(X)
  is_data_frame = is.data.frame(X)
  RAW = X
  if(is_vector){
    RAW = matrix(X, nrow=1)
  }
  if(is_data_frame){
    RAW = as.matrix(X)
  }
  ILR = .Call('normalmultinomial_ilr_coordinates', PACKAGE = 'normalmultinomial', RAW)
  colnames(ILR) = paste0('ilr_',1:ncol(ILR))

  if(is_vector){
    ILR = ILR[1,]
    names(ILR) = paste0('ilr_',1:length(ILR))
  }
  if(is_data_frame){
    ILR = as.data.frame(ILR)
  }
  class(ILR) = class_type
  ILR
}

#' @export
inv_ilr_coordinates = function(X, components.name = NULL){
  class_type = class(X)
  is_vector = is.vector(X)
  is_data_frame = is.data.frame(X)
  ILR = X
  if(is_vector){
    ILR = matrix(X, nrow=1)
  }
  if(is_data_frame){
    ILR = as.matrix(X)
  }
  RAW = .Call('normalmultinomial_inv_ilr_coordinates', PACKAGE = 'normalmultinomial', ILR)
  if(is.null(components.name)){
    components.name = paste0('c_', 1:(ncol(RAW)))
  }
  colnames(RAW) = components.name

  if(is_vector){
    RAW = RAW[1,]
    names(RAW) = components.name
  }
  if(is_data_frame){
    RAW = as.data.frame(RAW)
  }
  class(RAW) = class_type
  RAW
}

#' @export
clr_coordinates = function(X){
  components.name = colnames(X)
  class_type = class(X)
  is_vector = is.vector(X)
  is_data_frame = is.data.frame(X)
  RAW = X
  if(is_vector){
    RAW = matrix(X, nrow=1)
    components.name = names(X)
  }
  if(is_data_frame){
    RAW = as.matrix(X)
    components.name = names(X)
  }
  if(is.null(components.name)){
    components.name = paste0('c',1:ncol(RAW))
  }
  CLR = .Call('normalmultinomial_clr_coordinates', PACKAGE = 'normalmultinomial', RAW)
  colnames(CLR) = paste0('clr_',components.name)
  if(is_vector){
    CLR = CLR[1,]
    names(CLR) = paste0('clr_',components.name)
  }
  if(is_data_frame){
    CLR = as.data.frame(CLR)
  }
  class(CLR) = class_type
  CLR
}

#' @export
inv_clr_coordinates = function(X){
  components.name = colnames(X)
  class_type = class(X)
  is_vector = is.vector(X)
  is_data_frame = is.data.frame(X)
  CLR = X
  if(is_vector){
    CLR = matrix(X, nrow=1)
    components.name = names(X)
  }
  if(is_data_frame){
    CLR = as.matrix(X)
    components.name = names(X)
  }
  RAW = .Call('normalmultinomial_inv_clr_coordinates', PACKAGE = 'normalmultinomial', CLR)
  colnames(RAW) = gsub('clr_','',components.name)

  if(is_vector){
    RAW = RAW[1,]
    names(RAW) = gsub('clr_','',components.name)
  }
  if(is_data_frame){
    RAW = as.data.frame(RAW)
  }
  class(RAW) = class_type
  RAW
}
