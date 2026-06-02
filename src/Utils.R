SGDPLS <- function(X_train, Y_train, n_components=2, epochs=5, eta=0.000001) {
    write.csv(X_train, file="WEC_X_train.csv", row.names=FALSE)
    write.csv(Y_train, file="WEC_Y_train.csv", row.names=FALSE)
    system(paste("python ../src/SGDPLS.py --m", n_components, "--epochs", epochs, "--eta", eta))
    Beta <- as.matrix(read.csv("SGDPLS_Beta.csv", header=FALSE))
    time <- as.numeric(readLines("SGDPLS_time.txt"))
    return(list(Beta=Beta, time=time))
}

CIPLS <- function(X_train, Y_train, n_components=2) {
  write.csv(X_train, file="WEC_X_train.csv", row.names=FALSE)
  write.csv(Y_train, file="WEC_Y_train.csv", row.names=FALSE)
  system(paste("python ../src/CIPLS.py --m", n_components))
  Beta <- as.matrix(read.csv("CIPLS_Beta.csv", header=FALSE))
  time <- as.numeric(readLines("CIPLS_time.txt"))
  return(list(Beta=Beta, time=time))
}
