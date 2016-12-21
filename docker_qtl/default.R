#* @get /hello
hello <- function(name="World") {
  return (list("Hello", name))
}

function(req){
  print(paste0(date(), " - ",
               req$REMOTE_ADDR, " - ",
               req$REQUEST_METHOD, " ",
               req$PATH_INFO))
  forward()
}
