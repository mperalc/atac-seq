# sort df by two strings 
# the order of the strings is important
# 1st by "first", 2nd by "second"


sort_by_2_strings=function(df,first,second){
  order=c() # to store column order
  for (f in first){
    for(s in second){
      x=intersect(grep(paste("\\b",f,"\\b",sep=""), names(df)),
                  grep(paste("\\b",s,"\\b",sep=""), names(df))) # \\b are boundary anchors that allow to match whole words contained within a pattern
      # searching column numbers that match the two patterns
      order=c(order,x)
    }
  }
  df=df[,order]
  return(df)
}
