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
  colnames(ILR) =
  if(is_vector){
    ILR = ILR[1,]
    names(ILR) = paste0('ilr_',1:length(ILR))
  }
  if(is_data_frame){
    ILR = as.data.frame(ILR)
    names(ILR) = paste0('ilr_',1:ncol(ILR))
  }
  class(ILR) = class_type
  ILR
}

#' @export
inv_ilr_coordinates = function(X, components.name = paste0('c_', 1:(ncol(X)+1))){
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
  colnames(RAW) = components.name
  if(is_vector){
    RAW = RAW[1,]
    names(RAW) = components.name
  }
  if(is_data_frame){
    RAW = as.data.frame(RAW)
    names(RAW) = components.name
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
    names(CLR) = paste0('clr_',components.name)
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
    names(RAW) = gsub('clr_','',components.name)
  }
  class(RAW) = class_type
  RAW
}
