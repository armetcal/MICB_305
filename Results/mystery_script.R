generate = function(){
  
list_of_options = list(
  c("   ***   ",
    "  *   *  ",
    " *     * ",
    "*  ^ ^  *",
    "*  ___  *",
    "*  \\_/  *", 
    " *     * ",
    "  *****  "
  ),
  c(
    "  /\\_/\\  ",
    " ( o o ) ",
    "  > ^ <  "
  ),
  c(
    "   ><(((('>  ",
    "             ",
    "   ><(((('>  "
  ),
  c(
    " ,___,  ",
    " [O.o]  ",
    " /)__)  ",
    " --\"--  "
  ),
  c(
    "   .--.  ",
    "  |o_o | ",
    "  |:_/ | ",
    " //   \\ \\ ",
    "(|     | )",
    "/'\\_   _/`\\",
    "\\___)=(___/"
  ),
  c(
    "  ʕ•ᴥ•ʔ ",
    "  (    )  ",
    "  (____)  "
  ),
  c(
    "  / \\__  ",
    " (    @\\___  ",
    " /         O  ",
    "/   (_____ /  ",
    "/_____/   U   "
  )
)

to_print = sample(list_of_options,1)

cat(paste0(to_print[[1]],collapse='\n'))

# We're just printing to screen - no return value needed
# Therefore no return statement!
}
