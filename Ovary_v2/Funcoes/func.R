func_pose = function(pos_list){
  pos_1 = cbind(as.data.frame(pos_list[[1]]), as.data.frame(pos_list[[2]]))
  pos_2 = cbind(as.data.frame(pos_list[[3]]), as.data.frame(pos_list[[4]]))
  pos_3 = cbind(as.data.frame(pos_list[[5]]), as.data.frame(pos_list[[6]]))
  pos_4 = cbind(as.data.frame(pos_list[[7]]), as.data.frame(pos_list[[8]]))
  pos_5 = cbind(as.data.frame(pos_list[[9]]), as.data.frame(pos_list[[10]]))
  pos_6 = cbind(as.data.frame(pos_list[[11]]), as.data.frame(pos_list[[12]]))
  
  pos_E = list(pos_1,pos_2,pos_3,pos_4,pos_5, pos_6)
  return(as.data.frame(pos_E))
}
