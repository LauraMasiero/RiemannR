extractResponse <- function(response_list, data_set, events=1){ #data_set: non_adaptive_tr 
  if(events==1){
    labels <- as.vector(unlist(lapply(data_set, function(y) y$labels==1)))
    resp_list <- unlist(response_list)[labels==T]
    }
  else { if(events==0){
    labels <- as.vector(unlist(lapply(data_set, function(y) y$labels==0)))
    resp_list <- unlist(response_list)[labels==T] 
    }
  }
  resp_list <- resp_list 
  return(resp_list)
}