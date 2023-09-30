## ----------------------------------------------------------------------
##       ***** function for interleaving Adj and SP learning *****
##       Author: *************** December, 12, 2021 ***************
## ----------------------------------------------------------------------

`eemb` <- function(dataset,alpha,test = "zf", cluster = NULL,
                   max.sx = ncol(dataset), debug = FALSE){
  
  data.info = bnlearn:::check.data(dataset, allow.missing = TRUE)
  complete=data.info$complete.nodes
  
  fea_colname = colnames(dataset)
  
  Spouse = vector("list",ncol(dataset))
  names(Spouse) = fea_colname
  Adjs = vector("list",ncol(dataset))
  names(Adjs) = fea_colname
  MB = vector("list",ncol(dataset))
  names(MB) = fea_colname
  sepset <- lapply(colnames(dataset), function(.) {vec<-vector("list",ncol(dataset))
  names(vec)<-colnames(dataset)
  return(vec)})
  names(sepset)<-colnames(dataset)
  ci_num = 0
  
  # Calculate Adjacency set and spouse set
  for (target in fea_colname) {
    # if(target == "X3"){
    #   print("target is :")
    #   print(target)
    #   print(alpha)
    # }
    sp_num = vector("list",ncol(dataset))
    names(sp_num) = fea_colname
    
    
    SPs = vector("list",ncol(dataset))
    names(SPs) = fea_colname
    citest_num = 0
    cadj = character()
    adj = character()
    adj_num = c()
    dep_adj = numeric(length(fea_colname))
    dep_sp = numeric(length(fea_colname))
    NotAdj = c()
    
    # remove the nodes independent of the target node conditioning on empty set 
    # sort pvalue
    for (i in fea_colname) {
      if(!(identical(i,target))){
        s = character()
        
        citest_num = citest_num + 1
        pval1 <- bnlearn:::indep.test(target, i, sx = s, test = test, data = dataset,
                                      B = 0L, alpha = alpha, complete = complete)
        i_num <- (as.integer(gsub("X","",i)))
        dep_adj[i_num] = pval1[1]
        if(pval1[1] >= alpha){
          sepset[[target]][[i]] = s
        }
        else{
          cadj = append(cadj,i)
          adj_num = append(adj_num,i_num)
        }
      }
    }
    
    if(length(cadj)==0){
      next
    }
    
    depSort_adj = data.frame(cadj = cadj,dep=dep_adj[adj_num],num = adj_num)
    depSort_adj = depSort_adj[order(depSort_adj$dep),] 
    NotAdj = setdiff(fea_colname,union(cadj,target))
    
    
    ###############----------Phase1----------###############
    
    
    for (i in 1:length(cadj)) {
      adj = append(adj,depSort_adj$cadj[i])
      adj_length = length(adj)
      adj_tmp = adj
      last_break_flag = 0
      
      # remove false adj
      for (j in adj_length:1) {
        Y = adj[j]
        canadj = setdiff(adj_tmp,Y)
        cutSetSize = 1
        other_break_flag = 0
        
        
        while(length(canadj)>=cutSetSize){
          SS = combn(canadj,cutSetSize)
          
          for (si in 1:(ncol(SS))) {
            Z = SS[,si]
            if(depSort_adj$cadj[i] != Y){
              if(!(length(which(Z == depSort_adj$cadj[i])))){
                next
              }
            }
            citest_num = citest_num + 1
            pval = bnlearn:::indep.test(Y,target,sx=Z,test = test, data = dataset,
                                        B = 0L, alpha = alpha, complete = complete)
            
            # find spouses with regard to the Adj without Y
            
            if(pval[1] >= alpha){
              adj_tmp = canadj
              sepset[[target]][[Y]] = Z
              
              
              NotAdj = union(NotAdj,Y)
              for (o in 1:length(canadj)) {
                adj_var = canadj[o]
                if(length(which(sepset[[target]][[Y]] == adj_var))){next}
                
                S = union(sepset[[target]][[Y]],adj_var)
                citest_num = citest_num + 1
                pval1 = bnlearn:::indep.test(Y,target,sx=S,test = test, data = dataset,
                                             B = 0L, alpha = alpha, complete = complete)
                
                k_num <- (as.integer(gsub("X","",Y)))
                dep_sp[k_num] = pval1[1]
                
                
                if(pval1[1] <= alpha){
                  
                  SPs[[adj_var]] = append(SPs[[adj_var]],Y)
                  sp_num[[adj_var]] = append(sp_num[[adj_var]],k_num)
                  
                }
              }
              if(depSort_adj$cadj[i] == Y){last_break_flag=1}
              SPs[Y] = vector("list",1)
              names(SPs[Y]) = Y
              other_break_flag = 1
              break
            }
          }
          if(last_break_flag || other_break_flag){break}
          cutSetSize = cutSetSize + 1
        }
        if(last_break_flag == 1){break}
        
        # the new added node is belong to Adj
        # find spouses with regard to this node
        
        if(depSort_adj$cadj[i] == Y){
          for (k in NotAdj) {
            X = k
            S = union(sepset[[target]][[X]],Y)
            citest_num = citest_num + 1
            pval2 <- bnlearn:::indep.test(target, X, sx = S, test = test, data = dataset,
                                          B = 0L, alpha = alpha, complete = complete)
            
            k_num <- (as.integer(gsub("X","",X)))
            dep_sp[k_num] = pval2[1]
            
            if(pval2[1] <= alpha){
              SPs[[Y]] = append(SPs[[Y]],X)
              sp_num[[Y]] = append(sp_num[[Y]],k_num)
              
            }
          }
        }
      }
      adj = adj_tmp
    }
    
    
    # ################----------Phase2----------###############
    # step1:Remove false spouse positive
    
    for (i in 1:length(adj)) {
      Y = adj[i]
      
      if(length(SPs[[Y]])==0){
        Adjs[[target]] = adj_tmp
        next
      }
      
      sp = character()
      depSort_sp = data.frame(sp = SPs[[Y]],dep=dep_sp[sp_num[[Y]]],num = sp_num[[Y]])
      depSort_sp = depSort_sp[order(depSort_sp$dep),]
      
      
      #
      for (f in depSort_sp$sp) {
        sp = append(sp,f)
        sp_length = length(sp)
        sp_tmp = sp
        SP_break_flag = 0
        
        for (c in sp_length:1) {
          X = sp[c]
          cansp = setdiff(sp_tmp,X)
          cutSetSize = 1
      
          
          if(length(cansp) <= 0){
            condset = c()
            condset = union(condset,Y)
            citest_num = citest_num + 1
            pval3 <- bnlearn:::indep.test(X, target, sx=condset, test = test, data = dataset,
                                          B = 0L, alpha = alpha, complete = complete)
            
            if(pval3[1] >= alpha){
              sp_tmp = cansp
              
              if(f == X){last_SP_break_flag = 1}
              break
            }
          }else{
            SS = combn(cansp,cutSetSize)
            for(si in 1:(ncol(SS))){
              Z = SS[,si]
              condset = union(Z,Y)
              
              if(f != X){
                if(!(length(which(Z == f)))){
                  next
                }
              }
              
              citest_num = citest_num + 1
              pval3 <- bnlearn:::indep.test(X, target, sx=condset, test = test, data = dataset,
                                            B = 0L, alpha = alpha, complete = complete)
              
              
              if(pval3[1] >= alpha){
                sp_tmp = cansp
                
                if(f == X){SP_break_flag = 1}
                break
              }
            }
          }
          
          
          if(SP_break_flag == 1){break}
        }
        sp = sp_tmp
      }
      SPs[[Y]] = sp
      
      
      
      # 
      
      sp = character()
      for (f in SPs[[Y]]) {
        sp = append(sp,f)
        sp_length = length(sp)
        sp_tmp = sp
        SP_break_flag2 = 0
        
        for (c in sp_length:1) {
          X = sp[c]
          cansp = setdiff(sp_tmp,X)
          cutSetSize = 0
          other_SP_break_flag = 0
          
          while(length(cansp) >= cutSetSize){
            SS = combn(cansp,cutSetSize)
            for(si in 1:(ncol(SS))){
              Z = SS[,si]
              condset = Z
              
              if(f != X){
                if(!(length(which(Z == f)))){
                  next
                }
              }
              
              citest_num = citest_num + 1
              pval4 = bnlearn:::indep.test(X, Y, sx=condset, test = test, data = dataset,
                                           B = 0L, alpha = alpha, complete = complete)
              
              if(pval4[1] >= alpha){
                sp_tmp = setdiff(sp_tmp,X)
                
                if(f == X){SP_break_flag2=1}
                other_SP_break_flag=1
                break
              }
            }
            if((SP_break_flag2 == 1)||(other_SP_break_flag == 1)){break}
            cutSetSize = cutSetSize + 1
          }
          if(SP_break_flag2 == 1){break}
        }
        sp = sp_tmp
      }
      
      if(length(sp)==0){
        SPs[Y] = vector("list",1)
      }else{
        SPs[[Y]] = sp
      }
      
      # 
      
      sp = character()
      for (f in SPs[[Y]]) {
        sp = append(sp,f)
        sp_length = length(sp)
        sp_tmp = sp
        SP_break_flag3 = 0
        
        for (c in sp_length:1) {
          X = sp[c]
          cansp = setdiff(union(adj_tmp,sp_tmp),X)
          cutSetSize = 1
          other_SP_break_flag3 = 0
          
          while(length(cansp) >= cutSetSize){
            SS = combn(cansp,cutSetSize)
            for(si in 1:(ncol(SS))){
              Z = SS[,si]
              condset = unique(union(Z,Y))
              
              if(length(intersect(adj,condset)) == 1){next}
              
              
              if(f != X){
                if(!(length(which(Z == f)))){
                  next
                }
              }
              
              citest_num = citest_num + 1
              pval5 = bnlearn:::indep.test(X, target, sx=condset, test = test, data = dataset,
                                           B = 0L, alpha = alpha, complete = complete)
              
              if(pval5[1] >= alpha){
                sp_tmp = setdiff(sp_tmp,X)
                
                
                if(f == X){SP_break_flag3=1}
                other_SP_break_flag3=1
                break
              }
            }
            if((SP_break_flag3 == 1)||(other_SP_break_flag3 == 1)){break}
            cutSetSize = cutSetSize + 1
          }
          if(SP_break_flag3 == 1){break}
        }
        sp = sp_tmp
      }
      if(length(sp)==0){
        SPs[Y] = vector("list",1)
      }else{
        SPs[[Y]] = sp
      }
    }
    
    # phase3:Remove false positive from CAdj
    
    adj_tmp = adj
    for (Y in adj) {
      canadj = setdiff(adj_tmp,Y)
      cutsetsize = 1
      other_adj_break_flag=0
      while(length(canadj) >= cutSetSize){
        SS = combn(canadj,cutSetSize)
        for(si in 1:(ncol(SS))){
          Z = SS[,si]
          
          spouse_test = c()
          for (k in Z) {
            adj_var = k
            spouse_test = union(spouse_test,SPs[[adj_var]])
          }
          TestSet = union(Z,spouse_test)
          
          citest_num = citest_num + 1
          pval6 = bnlearn:::indep.test(Y, target, sx=TestSet, test = test, data = dataset,
                                       B = 0L, alpha = alpha, complete = complete)
          
          
          
          if(pval6[1] >= alpha){
            adj_tmp = canadj
            SPs[Y] = vector("list",1)
            names(SPs[Y]) = Y
            
            break
          }
          other_adj_break_flag=1
        }
        if(other_adj_break_flag == 1){break}
        cutSetSize = cutSetSize + 1
      }
    }
    Adjs[[target]] = adj_tmp
    
    
    for (ms in SPs){
      for (m in ms){
        if((!(is.null(m))) && (!(m %in% Spouse[[target]]))){
          Spouse[[target]] = append(Spouse[[target]],m)
        }
      }
    }
    ci_num = ci_num + citest_num
    # Adjs[[target]] = adj
  }
  
  # for (target in fea_colname) {
  #   for (X in Spouse[[target]]) {
  #     if(!(target %in% Spouse[[X]])){
  #       Spouse[[target]] = setdiff(Spouse[[target]],X)
  #     }
  #   }
  # }

  for (target in fea_colname) {
    for (X in Adjs[[target]]) {
      if(!(target %in% Adjs[[X]])){
        Adjs[[target]] = setdiff(Adjs[[target]],X)
      }
    }
  }
  
  
  

  return(list(adjs=Adjs,sps=Spouse))
}
