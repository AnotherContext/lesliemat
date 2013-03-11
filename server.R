require(shiny)
require(popbio)

shinyServer(function(input, output) {
  Leslie_calc <- reactive(function(){   
    
    model_type <- switch(input$model_type,
                         model1 = "model1",
                         model2 = "model2"
    );
    
    h  <- input$h;
    sj <- input$sj;
    sa <- input$sa; 
    ff <- input$ff;  
    #ff<-75
    #sa<-0.15
    #sj<-0.1
    #h<-1
    p1 <- sj;
    p2 <- (1+h)/2;
    p3 <- (1-h)/(h+1);
    F1 <- 0;
    F2 <- sa*ff*(1-h)/2;
    F3 <- sa*ff*2*h/(h+1);
    F4 <- sa*ff;
    
    A0 <- matrix(0, nrow=4, ncol=4)
    #A0[1, 1] <- F1
    A0[1, 2] <- F2 
    A0[1, 3] <- F3 
    A0[1, 4] <- F4 
    A0[2, 1] <- p1 
    A0[3, 2] <- p2 
    A0[4, 3] <- p3 
    
    #documentation of eigen() say they are sorted in decreasing order.
    #so the first item of ..$values is the largest eigenvalue and thus the first column of ..$vectors is the corresponding eigenvector. 
    
    #dominant eigenvalue
    lambda <- eigen(A0)$values[1]
    #right eigenvector
    right_vec <- eigen(A0)$vectors[, 1] 
    #left eigenvector
    left_vec <- eigen(t(A0))$vectors[, 1] 
    
    stable_age <- matrix(0, ncol=4)
    #stable age distribution    
    for(i in 1:4){
      stable_age[i]<- right_vec[i]/sum(right_vec)
    }
    
    #sensitivity matrix
    sm <- array(dim=c(4, 4)) 
    
    for (i in 1:4){
      for (j in 1:4){
        sm[i, j] <- A0[i, j] * left_vec[i] * right_vec[j] / (lambda * t(left_vec) %*% right_vec)
      }
    }
    
    sm <- matrix(as.numeric(sm), nrow=4, ncol=4)
    rownames(sm) <- c("Class 1", "Class 2", "Class 3", "Class 4")
    colnames(sm) <- rownames(sm)
    rownames(A0) <- c("A[1, ]", "A[2, ]", "A[3, ]", "A[4, ]")
    colnames(A0) <- c("A[ , 1]", "A[ , 2]", "A[ , 3]", "A[ , 4]")
    
    #ngens <- 500
    #n <- matrix(0, nrow=4,ncol=ngens+1);
    
    #c determines the strength of the density dependence for model2
    c_val <- input$c_value;
    
    # initial pop'n
    n0 <- c(1, 0, 0, 0)  
    n <- array(NA, dim=c(4, input$ngens))
    n[,1] <- n0
    
    #total population, project future growth
    Nt <- array(NA, dim=ncol(n)) 
    
    #proportion in each age class
    pn <- array(dim=c(4, ncol(n))) 
    
    if (model_type == "model1"){
      
      for (t in 2:input$ngens){
        n[,t] <- A0 %*% n[, t-1]
        Nt[t] <- sum(n[, t])
        pn[,t] <- n[, t]/Nt[t]
      }
    } else {
      for (t in 2:input$ngens){
        A0[1, 2] <- F2*exp(-c_val*(sum(n[2:4,t-1],1)))
        n[,t] <- A0 %*% n[, t-1]
        Nt[t] <- sum(n[, t])
        pn[,t] <- n[, t]/Nt[t]  
      }
    }
    
    list(A0=A0, sm=sm, x=n, Nt=Nt, px=pn, Lambda=lambda, stable_age=stable_age)
  })
  
  output$Lesliemat <- reactiveTable(function(){
    Leslie_calc()$A0
  })
  
  output$stableage<-reactivePrint(function(){
    Leslie_calc()$stable_age    
  })
  
  output$elastTab <- reactiveTable(function(){
    Leslie_calc()$sm
  })
  
  output$elastfunc<-reactivePrint(function(){
    eigen.analysis(Leslie_calc()$A0 , zero=T)   
  })
  
  output$plotmodel <- reactivePlot(function(){
    scale_type <- switch(input$scale_type,
                         logar = "logar",
                         linear = "linear"
    )
    
    out <- Leslie_calc()
    with(out, {
      title <- substitute(paste("Dominant eigenvalue = ", Lambda,
                                sep=""), list(Lambda=as.numeric(Lambda)))
      
      if (scale_type=="logar"){
        plot(x=1:ncol(x), y=Nt, type="l", xlab="Time", ylab="Log number of individuals", lty=2, log="y", 
             main=title)
      } else {
        plot(x=1:ncol(x), y=Nt, type="l", xlab="Time", ylab="Number of individuals", lty=2, 
             main=title)
      }
      lines(x=1:ncol(x), y=x[1, ], col="darkseagreen")
      lines(x=1:ncol(x), y=x[2, ], col="deepskyblue4")
      lines(x=1:ncol(x), y=x[3, ], col="palevioletred2")
      lines(x=1:ncol(x), y=x[4, ], col="red")
      ypos <- max(abs(Nt), na.rm=T)
      legend(x=400, y=ypos, c("N(t)","Class 1", "Class 2", "Class 3", "Class 4"),
             col=c("black","darkseagreen", "deepskyblue4", "palevioletred2", "red"),
             text.col="black", lty=c(1, 1, 1, 1))
    })  
  })
  
  
  output$plotstableage <- reactivePlot(function(){
    out <- Leslie_calc()
    with(out, {
      title <- substitute(paste("Dominant eigenvalue = ", Lambda,
                                sep=""), list(Lambda=as.numeric(Lambda)))
      plot(1:4,as.numeric(stable_age), xlab="Age", ylab="Stable Age Distribution", main=title, xaxt="n", col="white")
      axis(side=1, at=0:4)
      polygon(c(0.6, 0.6, 1.4, 1.4), c(0, as.numeric(stable_age[1]), as.numeric(stable_age[1]),0), col="darkseagreen")
      polygon(c(1.6, 1.6, 2.4, 2.4), c(0, as.numeric(stable_age[2]), as.numeric(stable_age[2]),0), col="deepskyblue4")
      polygon(c(2.6, 2.6, 3.4, 3.4), c(0, as.numeric(stable_age[3]), as.numeric(stable_age[3]),0), col="palevioletred2")
      polygon(c(3.6, 3.6, 4.4, 4.4), c(0, as.numeric(stable_age[4]), as.numeric(stable_age[4]),0), col="red")
      abline(h=0)
    })
  })
})