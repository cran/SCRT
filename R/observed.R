observed <-
function(design,statistic,data=read.table(file.choose(new=FALSE))){

  statAB<-function(A,B)
  {
    return(mean(A,na.rm=TRUE)-mean(B,na.rm=TRUE))
  }
  
  statBA<-function(A,B)
  {
    return(mean(B,na.rm=TRUE)-mean(A,na.rm=TRUE))
  }
  
  statabsAB<-function(A,B)
  {
    return(abs(mean(A,na.rm=TRUE)-mean(B,na.rm=TRUE)))
  }
  
  statPAPB<-function(A1,B1,A2,B2)
  {
    return(mean(c(mean(A1,na.rm=TRUE),mean(A2,na.rm=TRUE)),na.rm=TRUE)-mean(c(mean(B1,na.rm=TRUE),mean(B2,na.rm=TRUE)),na.rm=TRUE))
  }
  
  statPBPA<-function(A1,B1,A2,B2)
  {
    return(mean(c(mean(B1,na.rm=TRUE),mean(B2,na.rm=TRUE)),na.rm=TRUE)-mean(c(mean(A1,na.rm=TRUE),mean(A2,na.rm=TRUE)),na.rm=TRUE))
  }
  
  statabsPAPB<-function(A1,B1,A2,B2)
  {
    return(abs(mean(c(mean(A1,na.rm=TRUE),mean(A2,na.rm=TRUE)),na.rm=TRUE)-mean(c(mean(B1,na.rm=TRUE),mean(B2,na.rm=TRUE)),na.rm=TRUE)))
  }
  
  if(design=="CRD"|design=="RBD"|design=="ATD"|design=="AB"|design=="Custom"){
    observed.a<-data[,2][data[,1]=="A"]
    observed.b<-data[,2][data[,1]=="B"]
    
    if(statistic=="A-B"){
      observed.statistic<-statAB(observed.a,observed.b)
    }
    else if(statistic=="B-A"){
      observed.statistic<-statBA(observed.a,observed.b)
    }
    else if(statistic=="|A-B|"){
      observed.statistic<-statabsAB(observed.a,observed.b)
    }
    else{
      A<-observed.a
      B<-observed.b
      observed.statistic<-eval(parse(text=statistic))
    }
  }
  
  if(design=="ABA"){
    observed.a1<-data[,2][data[,1]=="A1"]
    observed.b1<-data[,2][data[,1]=="B1"]	
    observed.a2<-data[,2][data[,1]=="A2"]	
    observed.a<-c(observed.a1,observed.a2)	
    
    if(statistic=="A-B"){
      observed.statistic<-statAB(observed.a,observed.b1)
    }
    else if(statistic=="B-A"){
      observed.statistic<-statBA(observed.a,observed.b1)
    }
    else if(statistic=="|A-B|"){
      observed.statistic<-statabsAB(observed.a,observed.b1)
    }
    else if(statistic=="PA-PB"){
      observed.statistic<-statPAPB(observed.a1,observed.b1,observed.a2,NA)
    }
    else if(statistic=="PB-PA"){
      observed.statistic<-statPBPA(observed.a1,observed.b1,observed.a2,NA)
    }
    else if(statistic=="|PA-PB|"){
      observed.statistic<-statabsPAPB(observed.a1,observed.b1,observed.a2,NA)
    }
    else{
      A1<-observed.a1
      B1<-observed.b1	
      A2<-observed.a2
      A<-observed.a
      B<-observed.b1
      observed.statistic<-eval(parse(text=statistic))
    }
  }
  
  if(design=="ABAB"){
    observed.a1<-data[,2][data[,1]=="A1"]
    observed.b1<-data[,2][data[,1]=="B1"]	
    observed.a2<-data[,2][data[,1]=="A2"]	
    observed.b2<-data[,2][data[,1]=="B2"]	
    observed.a<-c(observed.a1,observed.a2)	
    observed.b<-c(observed.b1,observed.b2)	
    
    if(statistic=="A-B"){
      observed.statistic<-statAB(observed.a,observed.b)
    }
    else if(statistic=="B-A"){
      observed.statistic<-statBA(observed.a,observed.b)
    }
    else if(statistic=="|A-B|"){
      observed.statistic<-statabsAB(observed.a,observed.b)
    }
    else if(statistic=="PA-PB"){
      observed.statistic<-statPAPB(observed.a1,observed.b1,observed.a2,observed.b2)	
    }
    else if(statistic=="PB-PA"){
      observed.statistic<-statPBPA(observed.a1,observed.b1,observed.a2,observed.b2)
    }
    else if(statistic=="|PA-PB|"){
      observed.statistic<-statabsPAPB(observed.a1,observed.b1,observed.a2,observed.b2)
    }
    else{
      A1<-observed.a1
      B1<-observed.b1	
      A2<-observed.a2
      B2<-observed.b2
      A<-observed.a
      B<-observed.b
      observed.statistic<-eval(parse(text=statistic))
    }
  }
  
  if(design=="MBD"){
    observed.a<-list() 
    for(it in 1:(ncol(data)/2)){
      observed.a[[it]]<-data[,it*2][data[,(it*2)-1]=="A"]
    }
    observed.b<-list()
    for(it in 1:(ncol(data)/2)){
      observed.b[[it]]<-data[,it*2][data[,(it*2)-1]=="B"]
    }
    differences<-numeric(ncol(data)/2)
    
    if(statistic=="A-B"){
      for(it in 1:(ncol(data)/2)){
        differences[it]<-statAB(observed.a[[it]],observed.b[[it]])
      }
    }
    else if(statistic=="B-A"){
      for(it in 1:(ncol(data)/2)){
        differences[it]<-statBA(observed.a[[it]],observed.b[[it]])
      }
    }
    else if(statistic=="|A-B|"){
      for(it in 1:(ncol(data)/2)){
        differences[it]<-statabsAB(observed.a[[it]],observed.b[[it]])
      }
    }
    else{
      for(it in 1:(ncol(data)/2)){
        A<-observed.a[[it]]
        B<-observed.b[[it]]
        differences[it]<-eval(parse(text=statistic))
      }
    }
    observed.statistic<-mean(differences,na.rm=TRUE)
  }
  
  return(observed.statistic)
}
