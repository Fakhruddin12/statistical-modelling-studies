code <- readLines("analysis.R")
start <- grep("^parse_stop_hour <-", code)
end <- which(seq_along(code) > start & code == "}")[1]
eval(parse(text = code[start:end]))
stopifnot(identical(parse_stop_hour(c("930", "0930", "09:30", "30", "00:30", "2359", "24:00", "12:70", "bad")), c(9L, 9L, 9L, 0L, 0L, 23L, NA_integer_, NA_integer_, NA_integer_)))
cat("Nine time-parser examples passed. No model was fitted.\n")
