distribution.random <-
function(
  design,
  statistic,
  save="no",
  number,
  limit,
  data=read.table(file.choose(new=FALSE)),
  starts=file.choose(new=FALSE),
  assignments=file.choose(new=FALSE)
){
  
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
  
  if(design=="CRD"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    file.a<-tempfile(pattern="ascores",tmpdir=tempdir())
    file.b<-tempfile(pattern="bscores",tmpdir=tempdir())
    
    for(it in 1:number)
    {
      index<-if(it==1) data[,1] else sample(data[,1],replace=FALSE)
      scores.a<-data[,2][index=="A"]
      scores.a<-as.matrix(scores.a)
      scores.a<-t(scores.a)
      write.table(scores.a,file=file.a,append=TRUE,col.names=FALSE,row.names=FALSE)
      scores.b<-data[,2][index=="B"]
      scores.b<-as.matrix(scores.b)
      scores.b<-t(scores.b)
      write.table(scores.b,file=file.b,append=TRUE,col.names=FALSE,row.names=FALSE)
    }
    ascores<-read.table(file.a)
    bscores<-read.table(file.b)
    unlink(file.a,recursive=FALSE)
    unlink(file.b,recursive=FALSE)
    
    ascores<-as.matrix(ascores)
    bscores<-as.matrix(bscores)
    distribution<-numeric(number)
    
    if(statistic=="A-B"){
      for(it in 1:number){
        distribution[it]<-statAB(ascores[it,],bscores[it,])
      }
    }
    else if(statistic=="B-A"){
      for(it in 1:number){
        distribution[it]<-statBA(ascores[it,],bscores[it,])
      }
    }
    else if(statistic=="|A-B|"){
      for(it in 1:number){
        distribution[it]<-statabsAB(ascores[it,],bscores[it,])
      }
    }
    else{
      for(it in 1:number){
        A<-ascores[it,]
        B<-bscores[it,]
        distribution[it]<-eval(parse(text=statistic))
      }
    }
    
    distribution<-sort(distribution)
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE)
    }
    
    return(distribution)
  }
  
  if(design=="RBD"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    file.a<-tempfile(pattern="ascores",tmpdir=tempdir())
    file.b<-tempfile(pattern="bscores",tmpdir=tempdir())
    MT<-nrow(data)
    ab<-c("A","B")
    
    for(it in 1:number)
    {
      if(it==1) 
        index<-data[,1] 
      else{
        index<-numeric()
        repeat{
          index<-c(index,sample(ab,2,replace=FALSE))
          if(length(index)==MT)break
        }
      }
      scores.a<-data[,2][index=="A"]
      scores.a<-as.matrix(scores.a)
      scores.a<-t(scores.a)
      write.table(scores.a,file=file.a,append=TRUE,col.names=FALSE,row.names=FALSE)
      scores.b<-data[,2][index=="B"]
      scores.b<-as.matrix(scores.b)
      scores.b<-t(scores.b)
      write.table(scores.b,file=file.b,append=TRUE,col.names=FALSE,row.names=FALSE)
    }
    ascores<-read.table(file.a)
    bscores<-read.table(file.b)
    unlink(file.a,recursive=FALSE)
    unlink(file.b,recursive=FALSE)
    
    ascores<-as.matrix(ascores)
    bscores<-as.matrix(bscores)
    distribution<-numeric(number)
    
    if(statistic=="A-B"){
      for(it in 1:number){
        distribution[it]<-statAB(ascores[it,],bscores[it,])
      }
    }
    else if(statistic=="B-A"){
      for(it in 1:number){
        distribution[it]<-statBA(ascores[it,],bscores[it,])
      }
    }
    else if(statistic=="|A-B|"){
      for(it in 1:number){
        distribution[it]<-statabsAB(ascores[it,],bscores[it,])
      }
    }
    else{
      for(it in 1:number){
        A<-ascores[it,]
        B<-bscores[it,]
        distribution[it]<-eval(parse(text=statistic))
      }
    }
    
    distribution<-sort(distribution)
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE)
    }
    
    return(distribution)
  }
  
  if(design=="ATD"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    file.a<-tempfile(pattern="ascores",tmpdir=tempdir())
    file.b<-tempfile(pattern="bscores",tmpdir=tempdir())
    MT<-nrow(data)
    N<-c(rep(0,MT/2),rep(1,MT/2))
    
    for(i in 1:number)
    {
      if(i==1)
        index<-data[,1] 
      else{
        repeat{
          index<-rbind(sample(N,MT,replace=FALSE))
          check<-numeric()
          for(itr in 1:(MT-limit)){
            check2<-0
            for(it in itr:(itr+limit)){
              check2<-check2+index[,it]
            }
            check<-cbind(check,check2)
          }
          
          if(sum(check==(limit+1)|check==0)==0){
            for(it in 1:(length(index))){
              if(index[,it]==0){
                index[,it]<-"A"
              }
              else{
                index[,it]<-"B"
              }
            }
            break
          }
        }
      }
      
      scores.a<-data[,2][index=="A"]
      scores.a<-as.matrix(scores.a)
      scores.a<-t(scores.a)
      write.table(scores.a,file=file.a,append=TRUE,col.names=FALSE,row.names=FALSE)
      scores.b<-data[,2][index=="B"]
      scores.b<-as.matrix(scores.b)
      scores.b<-t(scores.b)
      write.table(scores.b,file=file.b,append=TRUE,col.names=FALSE,row.names=FALSE)
    }
    ascores<-read.table(file.a)
    bscores<-read.table(file.b)
    unlink(file.a,recursive=FALSE)
    unlink(file.b,recursive=FALSE)
    
    ascores<-as.matrix(ascores)
    bscores<-as.matrix(bscores)
    distribution<-numeric(number)
    
    if(statistic=="A-B"){
      for(it in 1:number){
        distribution[it]<-statAB(ascores[it,],bscores[it,])
      }
    }
    else if(statistic=="B-A"){
      for(it in 1:number){
        distribution[it]<-statBA(ascores[it,],bscores[it,])
      }
    }
    else if(statistic=="|A-B|"){
      for(it in 1:number){
        distribution[it]<-statabsAB(ascores[it,],bscores[it,])
      }
    }
    else{
      for(it in 1:number){
        A<-ascores[it,]
        B<-bscores[it,]
        distribution[it]<-eval(parse(text=statistic))
      }
    }
    
    distribution<-sort(distribution)
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE)
    }
    
    return(distribution)
  }
  
  if(design=="AB"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    
    observed<-data[,2]
    MT<-nrow(data)
    quantity<-choose(MT-2*limit+1,1)
    selection<-sample(1:quantity,number-1,replace=TRUE)
    index.a<-limit:(MT-limit)
    
    scores.a<-list(data[,2][data[,1]=="A"])
    for(it in 1:(number-1)){
      scores.a[[it+1]]<-c(observed[1:index.a[selection[it]]])
    }
    scores.b<-list(data[,2][data[,1]=="B"])
    for(it in 1:(number-1)){
      scores.b[[it+1]]<-c(observed[(1+index.a[selection[it]]):length(observed)])
    }
    distribution<-numeric(number)
    
    if(statistic=="A-B"){
      for(it in 1:number){
        distribution[it]<-statAB(scores.a[[it]],scores.b[[it]])
      }
    }
    else if(statistic=="B-A"){
      for(it in 1:number){
        distribution[it]<-statBA(scores.a[[it]],scores.b[[it]])
      }
    }
    else if(statistic=="|A-B|"){
      for(it in 1:number){
        distribution[it]<-statabsAB(scores.a[[it]],scores.b[[it]])
      }
    }
    else{
      for(it in 1:number){
        A<-scores.a[[it]]
        B<-scores.b[[it]]
        distribution[it]<-eval(parse(text=statistic))
      }
    }
    
    distribution<-sort(distribution)
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
    }
    return(distribution)
  }
  
  if(design=="ABA"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    
    observed<-data[,2]
    MT<-nrow(data)
    quantity<-choose(MT-3*limit+2,2)
    selection<-sample(1:quantity,number-1,replace=TRUE)
    index1<-1:(MT-3*limit+1)
    index2<-rev(index1)
    index.a<-numeric()
    for(it in 1:length(index1)){
      index.a<-c(index.a,(rep((limit-1+index1[it]),index2[it])))
    }
    index.b<-numeric()
    for(itr in index1){
      for(it in itr:(MT-3*limit+1)){
        index.b<-c(index.b,2*limit-1+it)
      }
    }
    
    scores.a1<-list(data[,2][data[,1]=="A1"])
    for(it in 1:(number-1)){
      scores.a1[[it+1]]<-c(observed[1:(index.a[selection[it]])])
    }
    scores.b1<-list(data[,2][data[,1]=="B1"])
    for(it in 1:(number-1)){
      scores.b1[[it+1]]<-c(observed[(1+index.a[selection[it]]):(index.b[selection[it]])])
    }
    scores.a2<-list(data[,2][data[,1]=="A2"])
    for(it in 1:(number-1)){
      scores.a2[[it+1]]<-c(observed[(1+index.b[selection[it]]):(MT)])
    }
    scores.a<-list()
    for(it in 1:number){
      scores.a[[it]]<-c(scores.a1[[it]],scores.a2[[it]])
    }
    distribution<-numeric(number)
    
    if(statistic=="A-B"){
      for(it in 1:number){
        distribution[it]<-statAB(scores.a[[it]],scores.b1[[it]])
      }
    }
    else if(statistic=="B-A"){
      for(it in 1:number){
        distribution[it]<-statBA(scores.a[[it]],scores.b1[[it]])
      }
    }
    else if(statistic=="|A-B|"){
      for(it in 1:number){
        distribution[it]<-statabsAB(scores.a[[it]],scores.b1[[it]])
      }
    }
    else if(statistic=="PA-PB"){
      for(it in 1:number){
        distribution[it]<-statPAPB(scores.a1[[it]],scores.b1[[it]],scores.a2[[it]],NA)
      }
    }
    else if(statistic=="PB-PA"){
      for(it in 1:number){
        distribution[it]<-statPBPA(scores.a1[[it]],scores.b1[[it]],scores.a2[[it]],NA)
      }
    }
    else if(statistic=="|PA-PB|"){
      for(it in 1:number){
        distribution[it]<-statabsPAPB(scores.a1[[it]],scores.b1[[it]],scores.a2[[it]],NA)
      }
    }
    else{
      for(it in 1:number){
        A1<-scores.a1[[it]]
        B1<-scores.b1[[it]]
        A2<-scores.a2[[it]]
        A<-scores.a[[it]]
        B<-scores.b1[[it]]
        distribution[it]<-eval(parse(text=statistic))
      }
    }
    
    distribution<-sort(distribution)
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
    }
    return(distribution)
  }
  
  if(design=="ABAB"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    
    observed<-data[,2]
    MT<-nrow(data)
    quantity<-choose(MT-4*limit+3,3)
    selection<-sample(1:quantity,number-1,replace=TRUE)
    index1<-1:(MT-4*limit+1)
    index2<-rev(cumsum(index1))
    
    index.a1<-numeric()
    for(it in 1:length(index1)){
      index.a1<-c(index.a1,(rep((limit+index1[it]-1),index2[it])))
    }
    scores.a1<-list(data[,2][data[,1]=="A1"])
    for(it in 1:(number-1)){
      scores.a1[[it+1]]<-c(observed[1:(index.a1[selection[it]])])
    }
    index.b1<-numeric()
    for(itr in index1){
      for(it in (itr-1):(MT-4*limit)){
        index.b1<-c(index.b1,rep((2*limit+it),(MT-4*limit+1-it)))
      }
    }
    scores.b1<-list(data[,2][data[,1]=="B1"])
    for(it in 1:(number-1)){
      scores.b1[[it+1]]<-c(observed[(1+index.a1[selection[it]]):index.b1[selection[it]]])
    }
    indexa2<-numeric()
    for(it in 1:length(index1)){
      indexa2<-c(indexa2,(index1[it]:length(index1)))
    }
    index.a2<-numeric()
    for(it in 1:length(indexa2)){
      index.a2<-c(index.a2,(4*limit-limit-1+(indexa2[it]:length(index1))))
    }
    scores.a2<-list(data[,2][data[,1]=="A2"])
    for(it in 1:(number-1)){
      scores.a2[[it+1]]<-c(observed[(1+index.b1[selection[it]]):index.a2[selection[it]]])
    }
    scores.b2<-list(data[,2][data[,1]=="B2"])
    for(it in 1:(number-1)){
      scores.b2[[it+1]]<-c(observed[(1+index.a2[selection[it]]):MT])
    }
    scores.a<-list()
    for(it in 1:number){
      scores.a[[it]]<-c(scores.a1[[it]],scores.a2[[it]])
    }
    scores.b<-list()
    for(it in 1:number){
      scores.b[[it]]<-c(scores.b1[[it]],scores.b2[[it]])
    }
    distribution<-numeric(number)
    
    if(statistic=="A-B"){
      for(it in 1:number){
        distribution[it]<-statAB(scores.a[[it]],scores.b[[it]])
      }
    }
    else if(statistic=="B-A"){
      for(it in 1:number){
        distribution[it]<-statBA(scores.a[[it]],scores.b[[it]])
      }
    }
    else if(statistic=="|A-B|"){
      for(it in 1:number){
        distribution[it]<-statabsAB(scores.a[[it]],scores.b[[it]])
      }
    }
    else if(statistic=="PA-PB"){
      for(it in 1:number){
        distribution[it]<-statPAPB(scores.a1[[it]],scores.b1[[it]],scores.a2[[it]],scores.b2[[it]])
      }
    }
    else if(statistic=="PB-PA"){
      for(it in 1:number){
        distribution[it]<-statPBPA(scores.a1[[it]],scores.b1[[it]],scores.a2[[it]],scores.b2[[it]])
      }
    }
    else if(statistic=="|PA-PB|"){
      for(it in 1:number){
        distribution[it]<-statabsPAPB(scores.a1[[it]],scores.b1[[it]],scores.a2[[it]],scores.b2[[it]])
      }
    }
    else{
      for(it in 1:number){
        A1<-scores.a1[[it]]
        B1<-scores.b1[[it]]
        A2<-scores.a2[[it]]
        B2<-scores.b2[[it]]
        A<-scores.a[[it]]
        B<-scores.b[[it]]
        distribution[it]<-eval(parse(text=statistic))
      }
    }
    
    distribution<-sort(distribution)
    if(save=="yes"|save=="check"){
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
    }
    return(distribution)
  }
  
  if(design=="MBD"){
    N<-ncol(data)/2
    MT<-nrow(data)
    readLines(con=starts,n=N)->startpoints
    limits<-strsplit(startpoints,"\\s")
    limits<-lapply(limits,function(x){x[x!=""]})
    
    numbers<-numeric(N)
    for(it in 1:N){
      numbers[it]<-length(limits[[it]])
    }
    fileCOMBSTARTPTS<-tempfile(pattern="startpoints",tmpdir=tempdir())
    repeat{
      startpt<-numeric(N)
      for(it in 1:N){
        if(numbers[it]!=1){
          startpt[it]<-sample(limits[[it]],1)
        }
        else{
          startpt[it]<-limits[[it]]
        }
      }
      selectdesign<-sample(startpt,replace=FALSE)
      selectdesign1<-rbind(selectdesign)
      write.table(selectdesign1,file=fileCOMBSTARTPTS,append=TRUE,col.names=FALSE,row.names=FALSE)
      combstartpts<-read.table(fileCOMBSTARTPTS)
      if(nrow(combstartpts)==(number-1))break
    }
    unlink(fileCOMBSTARTPTS,recursive=FALSE)
    
    ascores<-list()
    for(it in 1:N){
      ascores[[it]]<-data[,it*2][data[,it*2-1]=="A"]
    }
    scores.a<-list(ascores)
    for(iter in 1:(number-1)){
      ascores<-list()
      for(it in 1:N){
        ascores[[it]]<-data[1:(combstartpts[iter,it]-1),it*2]
      }
      scores.a[[iter+1]]<-ascores
    }
    bscores<-list()
    for(it in 1:N){
      bscores[[it]]<-data[,it*2][data[,it*2-1]=="B"]
    }
    scores.b<-list(bscores)
    for(iter in 1:(number-1)){
      bscores<-list()
      for(it in 1:N){
        bscores[[it]]<-data[combstartpts[iter,it]:MT,it*2]
      }
      scores.b[[iter+1]]<-bscores
    }
    differs<-list()
    
    if(statistic=="A-B"){
      for(iter in 1:number){
        differ<-numeric(N)
        for(it in 1:N){
          differ[it]<-statAB(scores.a[[iter]][[it]],scores.b[[iter]][[it]])
        }
        differs[[iter]]<-differ
      }
    }
    else if(statistic=="B-A"){
      for(iter in 1:number){
        differ<-numeric(N)
        for(it in 1:N){
          differ[it]<-statBA(scores.a[[iter]][[it]],scores.b[[iter]][[it]])
        }
        differs[[iter]]<-differ
      }
    }
    else if(statistic=="|A-B|"){
      for(iter in 1:number){
        differ<-numeric(N)
        for(it in 1:N){
          differ[it]<-statabsAB(scores.a[[iter]][[it]],scores.b[[iter]][[it]])
        }
        differs[[iter]]<-differ
      }
    }
    else{
      for(iter in 1:number){
        differ<-numeric(N)
        for(it in 1:N){
          A<-scores.a[[iter]][[it]]
          B<-scores.b[[iter]][[it]]
          differ[it]<-eval(parse(text=statistic))
        }
        differs[[iter]]<-differ
      }
    }
    
    distribution<-numeric(number)
    for(it in 1:number){
      distribution[it]<-mean(differs[[it]],na.rm=TRUE)
    }
    
    distribution<-sort(distribution)
    if(save=="yes"){
      fileSAVE<-file.choose(new=FALSE)
    }
    if(save=="yes"|save=="check"){
      write.table(distribution,file=fileSAVE,col.names=FALSE,row.names=FALSE,append=FALSE)
    }
    
    return(distribution)
  }

  if(design=="Custom"){
    if(save=="yes"){
      file<-file.choose(new=FALSE)
    }
    
    observed<-data[,2]
    assignments<-read.table(assignments)
    selection<-sample.int(nrow(assignments),number-1,replace=TRUE)
    
    scores.a<-list(data[,2][data[,1]=="A"])
    for(it in 1:(number-1))
      scores.a[[it+1]]<-c(observed[assignments[selection[it],]=="A"])
    
    scores.b<-list(data[,2][data[,1]=="B"])
    for(it in 1:(number-1))
      scores.b[[it+1]]<-c(observed[assignments[selection[it],]=="B"])
    
    distribution<-numeric(number)
    if(statistic=="A-B"){
      for(it in 1:number)
        distribution[it]<-statAB(scores.a[[it]],scores.b[[it]])
    }
    else if(statistic=="B-A"){
      for(it in 1:number)
        distribution[it]<-statBA(scores.a[[it]],scores.b[[it]])
    }
    else if(statistic=="|A-B|"){
      for(it in 1:number)
        distribution[it]<-statabsAB(scores.a[[it]],scores.b[[it]])
    }
    else{
      for(it in 1:number){
        A<-scores.a[[it]]
        B<-scores.b[[it]]
        distribution[it]<-eval(parse(text=statistic))
      }
    }
    
    distribution<-sort(distribution)
    if(save=="yes")
      write.table(distribution,file=file,col.names=FALSE,row.names=FALSE,append=FALSE)
    
    return(distribution)
  }
  
}
