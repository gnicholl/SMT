

biaggregate = function(namedmatrix) {
  if(!is.matrix(namedmatrix)) namedmatrix=as.matrix(namedmatrix)
  blah2 = aggregate(namedmatrix,list(rownames(namedmatrix)),sum)
  rownames(blah2) = blah2[,1]
  blah2 = blah2[2:(ncol(namedmatrix)+1)]
  colnames(blah2) = colnames(namedmatrix)
  blah3 = aggregate(t(blah2),list(colnames(blah2)),sum)
  rownames(blah3) = blah3[,1]
  blah3 = blah3[2:(nrow(blah2)+1)]
  return(t(blah3))
}
