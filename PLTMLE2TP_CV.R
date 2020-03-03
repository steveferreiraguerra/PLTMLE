PLTMLE2TP_CV = function(training,validation){
  
  ind22 = 0
  ind12 = 0
  ind11 = 0
  
  e_22 = 0
  e_12 = 0
  e_11 = 0

  dat = training
  
  t = 2 # number of time-points
  d = 3 # number of regimes
  
  #### Indicators of treatment ####
  
  #t = 2
  
  I2_11 = ifelse(is.na(dat$A0) | is.na(dat$A4) ,NA,
                 ifelse(dat$A0 == 1 & dat$A4 == 1 , 1, 0))
  
  I2_01 = ifelse(is.na(dat$A0) | is.na(dat$A4) ,NA,
                 ifelse(dat$A0 == 0 & dat$A4 == 1 , 1, 0))
  
  I2_00 = ifelse(is.na(dat$A0) | is.na(dat$A4) ,NA,
                 ifelse(dat$A0 == 0 & dat$A4 == 0 , 1, 0))
  
  I2 = c(I2_11,I2_01,I2_00)
  
  #t = 1
  
  I1_11 = ifelse(is.na(dat$A0),NA,
                 ifelse(dat$A0 == 1, 1, 0))
  
  I1_01 = ifelse(is.na(dat$A0) ,NA,
                 ifelse(dat$A0 == 0, 1, 0))
  
  I1_00 = ifelse(is.na(dat$A0),NA,
                 ifelse(dat$A0 == 0, 1, 0))
  
  I1 = c(I1_11,I1_01,I1_00)
  
  ###### Set regimes data ########
  
  dat_11 = dat
  dat_11[,grep("A", names(dat_11))] = 1
  
  dat_01 = dat
  dat_01[,grep("A", names(dat_01))] = c(rep(0,nrow(dat_01)),rep(1,nrow(dat_01)))
  
  dat_00 = dat
  dat_00[,grep("A", names(dat_00))] = 0
  
  dat_rep = rbind(dat_11,dat_01,dat_00)
  
  # Estimation of g models
  
  #A4
  g1 = glm(A4 ~ L0 + A0 + L4 , family = binomial,data = dat)
  g1_a = ifelse(dat_rep$A4 == 1, predict(g1,type = "response", newdata = dat_rep), 
                1-predict(g1,type = "response", newdata = dat_rep)) 
  
  #A0
  g0 = glm(A0 ~ L0, family = binomial,data = dat)
  g0_a = ifelse(dat_rep$A0 == 1, predict(g0,type = "response", newdata = dat_rep), 
                1-predict(g0,type = "response", newdata = dat_rep))
  
  # Cumulative g
  
  g10_a = g1_a*g0_a
  
  # Truncation
  
  bound_up = quantile(c(g10_a,g0_a), c(.95), na.rm = T)
  bound_low = quantile(c(g10_a,g0_a), c(.05), na.rm = T)
  
  g10_a[which(g10_a > bound_up)] = bound_up
  g0_a[which(g0_a > bound_up)] = bound_up
  
  g10_a[which(g10_a < bound_low)] = bound_low
  g0_a[which(g0_a < bound_low)] = bound_low
  
  # Step 1 - Initial estimation of Q22 -> This is for t = 2, j = 2
  
  Q22 = glm(Y12 ~ L0 + A0 + L4 + A4, family = "binomial",data = dat)
  dat_rep$Q22 = predict(Q22,type = "response",newdata = dat_rep)
  
  # Step 2 - Update to Q22*
  
  # repeated copy of Y2 = Q32*
  
  Y = dat_rep$Y12
  
  # create weighted covariate
  
  S1 = 1
  S2 = 2
  S3 = c(rep(1,nrow(dat)),rep(1,nrow(dat)),rep(0,nrow(dat)))
  off = logit(dat_rep$Q22)
  
  update_data_22 = as.data.frame(cbind(Y,S1,S2,S3,off))
  colnames(update_data_22) = c("Y","S1","S2","S3","off")
  
  update_data_22$I2 = I2
  update_data_22$H = (1/g10_a)*update_data_22$I2
  
  # Estimation of epsilon
  
  eps_22 = glm(Y ~ -1 + S1 + S2 + S3 + offset(off) , family=quasibinomial,data = update_data_22, 
               weights=scale(H,center = F), subset = (I2 == 1))
  
  # Generate Q22*
  
  Q22_star = predict(eps_22, type = "response", newdata = update_data_22)
  Q22_star[which(is.na(Q22_star))] = 1
  
  # See if solved score --> if not, attempt solving numerically --> if not, NA
  
  score = (((Y - Q22_star)*cbind(S1,S2,S3))/g10_a)[I2==1,]
  if(sum(abs(colSums(score,na.rm = T))) > 0.001){

  # Create score function to be solved  
        
    ScoreFunction = function(e) {
      Q22_star = plogis(off + cbind(S1,S2,S3)%*% e)
      Q22_star[which(is.na(Q22_star))] = 1
      score = (((Y - Q22_star[,1])*cbind(S1,S2,S3))/g10_a)[I2==1,]
      return(sum(abs(colSums(score,na.rm = T)))) #each column has to add to zero
    }
    
  # Use same solver (nlminb) as in the ltmle package  
    
    FindEps = function(minimizer) {
      num.tries = 20
      init.eps = c(0,0,0) # initiate epsilon to 0
      for (i in 1:num.tries) {
        m = nlminb(start=init.eps, objective=ScoreFunction, control=list(abs.tol=0.0001, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
        e = m$par
        obj.val = m$objective
        if (obj.val < 0.0001) {
          return(list(e=e, solved=TRUE, m=m))
        }
        init.eps = rnorm(3) # random initial epsilon, if epsilon = 0 did not work  
      }
      return(list(e=c(0,0,0), solved=FALSE)) # returns Q22* :=  Q22 (not updated)
    }
    
    mm = FindEps(nlminb)
    
    if (mm$solved){
      Q22_star = plogis(off + cbind(S1,S2,S3)%*% mm$e)[,1]
      Q22_star[which(is.na(Q22_star))] = 1
    }   else {
    paramList = list("varIC_final" = c(NA,NA,NA))
    return(paramList)
    stop("Did not converge")
    }
    ind22 = 1
    e_22 = mm$e 
  }
  
  ### Step 3 - Estimation of Q12 -> For t = 2, j = 1

  dat$Q22_star_11 = Q22_star[1:nrow(dat)]
  dat$Q22_star_01 = Q22_star[(nrow(dat) + 1):((d-1)*nrow(dat))]
  dat$Q22_star_00 = Q22_star[((d-1)*nrow(dat) + 1):(d*nrow(dat))]
  
  QQ12_11 = glm(Q22_star_11 ~ L0 + A0 , family = "quasibinomial",data = dat)
  QQ12_01 = glm(Q22_star_01 ~ L0 + A0 , family = "quasibinomial",data = dat)
  QQ12_00 = glm(Q22_star_00 ~ L0 + A0 , family = "quasibinomial",data = dat)
  
  Q12_11 = predict(QQ12_11,type = "response",newdata = dat_11)
  Q12_01 = predict(QQ12_01,type = "response",newdata = dat_01)
  Q12_00 = predict(QQ12_00,type = "response",newdata = dat_00)
  
  dat_rep$Q12 = c(Q12_11,Q12_01,Q12_00)
  
  # Step 4 - Update to Q12*
  
  # Construction of weighted covariate
  
  Y = Q22_star
  S1 = 1
  S2 = 2
  S3 = c(rep(1,nrow(dat)),rep(1,nrow(dat)),rep(0,nrow(dat)))
  off = logit(dat_rep$Q12)
  
  update_data_12 = as.data.frame(cbind(Y,S1,S2,S3,off))
  colnames(update_data_12) = c("Y","S1","S2","S3","off")
  
  update_data_12$I1 = I1
  update_data_12$H = (1/g0_a)*update_data_12$I1
  
  # MLE estimation of epsilon
  
  eps_12 = glm(Y ~ -1 + S1 + S2 + S3 + offset(off) , family=quasibinomial,data = update_data_12, 
               weights=scale(H,center = F), subset = (I1 == 1))
  
  # Generate Q22*
  
  Q12_star = predict(eps_12, type = "response", newdata = update_data_12)
  
  # See if solved score --> if not, attempt solving numerically --> if not, NA
  
  score = (((Y - Q12_star)*cbind(S1,S2,S3))/g0_a)[I1==1,]
  if(sum(abs(colSums(score,na.rm = T))) > 0.01){
  
  # Create score function to be solved  
    
    ScoreFunction = function(e) {
      Q12_star = plogis(off + cbind(S1,S2,S3)%*% e)
      Q12_star[which(is.na(Q12_star))] = 1
      score = (((Y - Q12_star[,1])*cbind(S1,S2,S3))/g0_a)[I1==1,]
      return(sum(abs(colSums(score,na.rm = T)))) #each column has to add to zero
    }
    
  # Use same solver as in the ltmle package  
    
    FindEps = function(minimizer) {
      num.tries = 20
      init.eps = c(0,0,0) # initiate epsilon to 0
      for (i in 1:num.tries) {
        m = nlminb(start=init.eps, objective=ScoreFunction, control=list(abs.tol=0.0001, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
        e = m$par
        obj.val = m$objective
        if (obj.val < 0.0001) {
          return(list(e=e, solved=TRUE, m=m))
        }
        init.eps = rnorm(3) # random initial epsilon, if epsilon = 0 did not work  
      }
      return(list(e=c(0,0,0), solved=FALSE)) # returns Q12* :=  Q12 (not updated)
    }
    
    mm = FindEps(nlminb)
    
    if (mm$solved){
      Q12_star = plogis(off + cbind(S1,S2,S3)%*% mm$e)[,1]
      Q12_star[which(is.na(Q12_star))] = 1
    }   else {
    paramList = list("varIC_final" = c(NA,NA,NA))
    return(paramList)
    stop("Did not converge")}
    
    ind12 = 1
    e_12 = mm$e 
  }
  
  ###  For t = 1, j = 1 ###
  
  # Step 1 - Initial estimation of Q11
  
  Q11 = glm(Y4 ~ L0 + A0, family = "binomial",data = dat)
  dat_rep$Q11 = predict(Q11,type = "response",newdata = dat_rep)
  
  # Step 2 - Update to Q11*
  
  Y = dat_rep$Y4
  
  S1 = 1
  S2 = 1
  S3 = c(rep(1,nrow(dat)),rep(0,nrow(dat)*2))
  off = logit(dat_rep$Q11)
  
  update_data_11 = as.data.frame(cbind(Y,S1,S2,S3,off))
  colnames(update_data_11) = c("Y","S1","S2","S3","off")
  
  update_data_11$I1 = I1
  update_data_11$H = (1/g0_a)*update_data_11$I1
  
  eps_11 = glm(Y ~ -1 + S1 + S2 + S3 + offset(off) , family=quasibinomial,data = update_data_11, 
               weights=scale(H,center = F), subset = (I1 == 1))

  Q11_star = predict(eps_11, type = "response", newdata = update_data_11)
  
  # See if solved score --> if not, attempt solving numerically --> if not, NA
  
  score = (((Y - Q11_star)*cbind(S1,S2,S3))/g0_a)[I1==1,]
  if(sum(abs(colSums(score,na.rm = T))) > 0.01){
  
  # Create score function to be solved  
    
    ScoreFunction = function(e) {
      Q11_star = plogis(off + cbind(S1,S2,S3)%*% e)
      Q11_star[which(is.na(Q11_star))] = 1
      score = (((Y - Q11_star[,1])*cbind(S1,S2,S3))/g0_a)[I1==1,]
      return(sum(abs(colSums(score,na.rm = T)))) #each column has to add to zero
    }
    
  # Use same solver as in the ltmle package  
    
    FindEps = function(minimizer) {
      num.tries = 20
      init.eps = c(0,0,0) # initiate epsilon to 0
      for (i in 1:num.tries) {
        m = nlminb(start=init.eps, objective=ScoreFunction, control=list(abs.tol=0.0001, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
        e = m$par
        obj.val = m$objective
        if (obj.val < 0.0001) {
          return(list(e=e, solved=TRUE, m=m))
        }
        init.eps = rnorm(3) # random initial epsilon, if epsilon = 0 did not work  
      }
      return(list(e=c(0,0,0), solved=FALSE)) # returns Q11* :=  Q11 (not updated)
    }
    
    mm = FindEps(nlminb)
    
    if (mm$solved){
      Q11_star = plogis(off + cbind(S1,S2,S3)%*% mm$e)[,1]
      Q11_star[which(is.na(Q11_star))] = 1
    }   else {
    paramList = list("varIC_final" = c(NA,NA,NA))
    return(paramList)
    stop("Did not converge")}

    
    ind11 = 1
    e_11 = mm$e 
  }
  
  Q1_star = c(Q11_star, Q12_star)
  
  
  ### Final Step - MSM
  
  S1 = 1 
  S2 = c(rep(1,nrow(dat)*d), rep(2,nrow(dat)*d))
  S3 = c(rep(1,nrow(dat)),rep(0,nrow(dat)*2),rep(1,nrow(dat)*(2)),rep(0,nrow(dat)))
  
  msm_data = as.data.frame(cbind(Q1_star,S1,S2,S3))
  colnames(msm_data) = c("Q1_star","S1","S2","S3")
  
  msm = glm(Q1_star ~ -1 + S1 + S2 + S3, family = "quasibinomial", data = msm_data)
  
  ####################################### VALIDATION ##############################################
  
  dat = validation
  
  #### Indicators of treatment ####
  
  #t = 2
  
  I2_11 = ifelse(is.na(dat$A0) | is.na(dat$A4) ,NA,
                 ifelse(dat$A0 == 1 & dat$A4 == 1 , 1, 0))
  
  I2_01 = ifelse(is.na(dat$A0) | is.na(dat$A4) ,NA,
                 ifelse(dat$A0 == 0 & dat$A4 == 1 , 1, 0))
  
  I2_00 = ifelse(is.na(dat$A0) | is.na(dat$A4) ,NA,
                 ifelse(dat$A0 == 0 & dat$A4 == 0 , 1, 0))
  
  I2 = c(I2_11,I2_01,I2_00)
  
  #t = 1
  
  I1_11 = ifelse(is.na(dat$A0),NA,
                 ifelse(dat$A0 == 1, 1, 0))
  
  I1_01 = ifelse(is.na(dat$A0) ,NA,
                 ifelse(dat$A0 == 0, 1, 0))
  
  I1_00 = ifelse(is.na(dat$A0),NA,
                 ifelse(dat$A0 == 0, 1, 0))
  
  I1 = c(I1_11,I1_01,I1_00)
  
  ###### Set regimes data ########
  
  dat_11 = dat
  dat_11[,grep("A", names(dat_11))] = 1
  
  dat_01 = dat
  dat_01[,grep("A", names(dat_01))] = c(rep(0,nrow(dat_01)),rep(1,nrow(dat_01)))
  
  dat_00 = dat
  dat_00[,grep("A", names(dat_00))] = 0
  
  dat_rep = rbind(dat_11,dat_01,dat_00)
  
  # Estimation of g models
  
  #A4
  g1_a = ifelse(dat_rep$A4 == 1, predict(g1,type = "response", newdata = dat_rep), 
                1-predict(g1,type = "response", newdata = dat_rep)) 
  
  #A0
  g0_a = ifelse(dat_rep$A0 == 1, predict(g0,type = "response", newdata = dat_rep), 
                1-predict(g0,type = "response", newdata = dat_rep))
  
  # Cumulative g
  
  g10_a = g1_a*g0_a
  
  # Truncation
  
  bound_up = quantile(c(g10_a,g0_a), c(.95), na.rm = T)
  bound_low = quantile(c(g10_a,g0_a), c(.05), na.rm = T)
  
  g10_a[which(g10_a > bound_up)] = bound_up
  g0_a[which(g0_a > bound_up)] = bound_up
  
  g10_a[which(g10_a < bound_low)] = bound_low
  g0_a[which(g0_a < bound_low)] = bound_low
  
  # t = 2, j = 2

  dat_rep$Q22 = predict(Q22,type = "response",newdata = dat_rep)

  Y = dat_rep$Y12
  
  S1 = 1
  S2 = 2
  S3 = c(rep(1,nrow(dat)),rep(1,nrow(dat)),rep(0,nrow(dat)))
  off = logit(dat_rep$Q22)
  
  update_data_22 = as.data.frame(cbind(Y,S1,S2,S3,off))
  colnames(update_data_22) = c("Y","S1","S2","S3","off")
  
  update_data_22$I2 = I2
  update_data_22$H = (1/g10_a)*update_data_22$I2
  
  # Generate Q22* - Using either eps from logistic model or solved numerically

  if (ind22 == 0){
    Q22_star = predict(eps_22, type = "response", newdata = update_data_22)
    Q22_star[which(is.na(Q22_star))] = 1
  } else {      
    Q22_star = plogis(off + cbind(S1,S2,S3)%*% e_22)[,1]
    Q22_star[which(is.na(Q22_star))] = 1
  }  
  
  # Storing for 1st component of variance of EIC 
  
  temp_score22_11 = (((Y - Q22_star)*cbind(S1,S2,S3))/g10_a)[1:nrow(dat),]*I2_11 
  #colSums(temp_score22_11,na.rm = T)
  temp_score22_01 = (((Y - Q22_star)*cbind(S1,S2,S3))/g10_a)[(nrow(dat) + 1):((d-1)*nrow(dat)),]*I2_01
  #colSums(temp_score22_01,na.rm = T)
  temp_score22_00 = (((Y - Q22_star)*cbind(S1,S2,S3))/g10_a)[((d-1)*nrow(dat) + 1):(d*nrow(dat)),]*I2_00
  
  # t = 2, j = 1
  
  Q12_11 = predict(QQ12_11,type = "response",newdata = dat_11)
  Q12_01 = predict(QQ12_01,type = "response",newdata = dat_01)
  Q12_00 = predict(QQ12_00,type = "response",newdata = dat_00)
  
  dat_rep$Q12 = c(Q12_11,Q12_01,Q12_00)

  Y = Q22_star
  S1 = 1
  S2 = 2
  S3 = c(rep(1,nrow(dat)),rep(1,nrow(dat)),rep(0,nrow(dat)))
  off = logit(dat_rep$Q12)
  
  update_data_12 = as.data.frame(cbind(Y,S1,S2,S3,off))
  colnames(update_data_12) = c("Y","S1","S2","S3","off")
  
  update_data_12$I1 = I1
  update_data_12$H = (1/g0_a)*update_data_12$I1
  
  # Using either eps from logistic model or solved numerically
  
  if (ind12 == 0){
    Q12_star = predict(eps_12, type = "response", newdata = update_data_12)
    Q12_star[which(is.na(Q12_star))] = 1
  } else {      
    Q12_star = plogis(off + cbind(S1,S2,S3)%*% e_12)[,1]
    Q12_star[which(is.na(Q12_star))] = 1
  } 
  
  # Storing for 1st component of variance of EIC
  
  temp_score12_11 = (((Y - Q12_star)*cbind(S1,S2,S3))/g0_a)[1:nrow(dat),]*I1_11
  #colSums(temp_score12_11,na.rm = T)
  temp_score12_01 = (((Y - Q12_star)*cbind(S1,S2,S3))/g0_a)[(nrow(dat) + 1):((d-1)*nrow(dat)),]*I1_01
  #colSums(temp_score12_01,na.rm = T)
  temp_score12_00 = (((Y - Q12_star)*cbind(S1,S2,S3))/g0_a)[((d-1)*nrow(dat) + 1):(d*nrow(dat)),]*I1_00
  
  # t = 1, j = 1
  
  dat_rep$Q11 = predict(Q11,type = "response",newdata = dat_rep)
  
  Y = dat_rep$Y4

  S1 = 1
  S2 = 1
  S3 = c(rep(1,nrow(dat)),rep(0,nrow(dat)*2))
  off = logit(dat_rep$Q11)
  
  update_data_11 = as.data.frame(cbind(Y,S1,S2,S3,off))
  colnames(update_data_11) = c("Y","S1","S2","S3","off")
  
  update_data_11$I1 = I1
  update_data_11$H = (1/g0_a)*update_data_11$I1

  # Using either eps from logistic model or solved numerically
  
  if (ind11 == 0){
    Q11_star = predict(eps_11, type = "response", newdata = update_data_11)
    Q11_star[which(is.na(Q11_star))] = 1
  } else {      
    Q11_star = plogis(off + cbind(S1,S2,S3)%*% e_11)[,1]
    Q11_star[which(is.na(Q11_star))] = 1
  } 
  
  # Storing for 1st component of variance of EIC
  
  temp_score11_11 = (((Y - Q11_star)*cbind(S1,S2,S3))/g0_a)[1:nrow(dat),]*I1_11 #pour le regime (1,1), temps = 1, k = 1
  #colSums(temp_score11_11,na.rm = T)
  temp_score11_01 = (((Y - Q11_star)*cbind(S1,S2,S3))/g0_a)[(nrow(dat) + 1):((d-1)*nrow(dat)),]*I1_01
  #colSums(temp_score11_01,na.rm = T)
  temp_score11_00 = (((Y - Q11_star)*cbind(S1,S2,S3))/g0_a)[((d-1)*nrow(dat) + 1):(d*nrow(dat)),]*I1_00
  
  ########## Variance ##############
  
  #Sum over all a_bar
  
  # For a specific t and j : Sum regimes (1,1), (0,1),  (0,0) for all summary measures (columns)
  
  temp_score11 = cbind(rowSums(cbind(temp_score11_11[,1],temp_score11_01[,1],temp_score11_00[,1]), na.rm=TRUE),
                       rowSums(cbind(temp_score11_11[,2],temp_score11_01[,2],temp_score11_00[,2]), na.rm=TRUE),
                       rowSums(cbind(temp_score11_11[,3],temp_score11_01[,3],temp_score11_00[,3]), na.rm=TRUE))
  
  temp_score12 = cbind(rowSums(cbind(temp_score12_11[,1],temp_score12_01[,1],temp_score12_00[,1]), na.rm=TRUE),
                       rowSums(cbind(temp_score12_11[,2],temp_score12_01[,2],temp_score12_00[,2]), na.rm=TRUE),
                       rowSums(cbind(temp_score12_11[,3],temp_score12_01[,3],temp_score12_00[,3]), na.rm=TRUE))
  
  temp_score22 = cbind(rowSums(cbind(temp_score22_11[,1],temp_score22_01[,1],temp_score22_00[,1]), na.rm=TRUE),
                       rowSums(cbind(temp_score22_11[,2],temp_score22_01[,2],temp_score22_00[,2]), na.rm=TRUE),
                       rowSums(cbind(temp_score22_11[,3],temp_score22_01[,3],temp_score22_00[,3]), na.rm=TRUE))
  
  #Sum for j & t
  
  temp_score = temp_score11 + temp_score12 + temp_score22
  
  ####### 2nd component of EIC ######
  
  Q1_star = c(Q11_star, Q12_star)
  
  ### MSM
  
  S1 = 1 
  S2 = c(rep(1,nrow(dat)*d), rep(2,nrow(dat)*d))
  S3 = c(rep(1,nrow(dat)),rep(0,nrow(dat)*2),rep(1,nrow(dat)*(2)),rep(0,nrow(dat)))
  
  msm_data = as.data.frame(cbind(Q1_star,S1,S2,S3))
  colnames(msm_data) = c("Q1_star","S1","S2","S3")

  Y_hat = predict(msm, type = "response", newdata = msm_data) 
  pred_diff = (Q1_star - Y_hat)*cbind(S1,S2,S3)
  
  # Split Q1_star for t = 1, 2
  
  pred_diff11_11 = pred_diff[1:nrow(dat),]
  pred_diff11_01 = pred_diff[(nrow(dat) + 1):((d-1)*nrow(dat)),]
  pred_diff11_00 = pred_diff[((d-1)*nrow(dat) + 1):(d*nrow(dat)),]
  pred_diff12_11 = pred_diff[(d*nrow(dat) + 1):((d+1)*nrow(dat)),]
  pred_diff12_01 = pred_diff[((d+1)*nrow(dat) + 1):((d+2)*nrow(dat)),]
  pred_diff12_00 = pred_diff[((d+2)*nrow(dat) + 1):((d+3)*nrow(dat)),]
  
  # Sum over all a_bar
  
  # For a specific t and j : Sum regimes (1,1), (0,1),  (0,0) for all summary measures (columns)
  
  pred_diff11 = cbind(rowSums(cbind(pred_diff11_11[,1],pred_diff11_01[,1],pred_diff11_00[,1]), na.rm=TRUE),
                      rowSums(cbind(pred_diff11_11[,2],pred_diff11_01[,2],pred_diff11_00[,2]), na.rm=TRUE),
                      rowSums(cbind(pred_diff11_11[,3],pred_diff11_01[,3],pred_diff11_00[,3]), na.rm=TRUE))
  
  pred_diff12 = cbind(rowSums(cbind(pred_diff12_11[,1],pred_diff12_01[,1],pred_diff12_00[,1]), na.rm=TRUE),
                      rowSums(cbind(pred_diff12_11[,2],pred_diff12_01[,2],pred_diff12_00[,2]), na.rm=TRUE),
                      rowSums(cbind(pred_diff12_11[,3],pred_diff12_01[,3],pred_diff12_00[,3]), na.rm=TRUE))
  
  #colSums(pred_diff11,na.rm = T)
  #colSums(pred_diff12,na.rm = T)
  
  # Sum for t
  
  pred_diff_tot = pred_diff11 + pred_diff12
  
  # Add 2 components together
  
  IC = pred_diff_tot + temp_score
  
  varIC = var(IC)
  
  #### BREAD ###
  
  bread = t(cbind(S1,S2,S3))%*%(Y_hat*(1-Y_hat)*cbind(S1,S2,S3))/nrow(dat)
  
  breadinv = solve(bread)
  
  varIC_final = diag(breadinv%*%varIC%*%t(breadinv)/nrow(dat))
  
  paramList = list("varIC_final" = varIC_final)
  
  return(paramList)
  
}
