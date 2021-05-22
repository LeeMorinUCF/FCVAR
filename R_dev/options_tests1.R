

# Function to test passing dots argument to
# update default options with arguments passed through dots.


options_practice <- function(...) {


  # Goal is to allow the user to pass a list of arguments in ...
  # then these will be provided.
  # Otherwise, the defaults will be used.

  dots <- list(...)
  print("dots = ")
  print(dots)

  print("names(dots) = ")
  print(names(dots))

  for (i in names(dots)) {
    print(i)
    print(unname(unlist(dots[i])))
  }

  if ('dogs' %in% names(dots)) {
    dogs_in <- unname(unlist(dots['dogs']))
  } else {
    dogs_in <- 23
  }

  if ('cats' %in% names(dots)) {
    cats_in <- unlist(unname(dots['cats']))
  } else {
    cats_in <- 9
  }

  if ('guinea_pigs' %in% names(dots)) {
    guinea_pigs_in <- unname(unlist(dots['guinea_pigs']))
  } else {
    guinea_pigs_in <- 'yummy'
  }

  if ('rhinoceroses' %in% names(dots)) {
    rhinoceroses_in <- unname(unlist(dots['rhinoceroses']))
  } else {
    rhinoceroses_in <- 'huge'
  }


  out <- list(dogs = dogs_in,
              cats = cats_in,
              guinea_pigs = guinea_pigs_in,
              rhinoceroses = rhinoceroses_in)

  return(out)
}


test_out <- options_practice(cats = 7, dogs = '5', guinea_pigs = 'fat')
test_out

attributes(test_out)




test_out <- options_practice(dogs = 999, guinea_pigs = 'funny', rhinoceroses = 42)
test_out


test_out <- options_practice()
test_out



