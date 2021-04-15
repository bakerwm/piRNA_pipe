#!/usr/bin/env Rscripts

## make qc stat
## save to a html file

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3){
  print("Usage: Rscript report.R <project.dir> <out.dir> <template>")
  print("")
  print("Option:")
  print("  project.dir     The directory of piRNA_pipe output")
  print("     out.html     The filename of html report")
  print("     template     The Rmarkdown template")
  stop("arguments failed")
}

input    <- args[1]
out_html <- args[2]
template <- args[3]

rmarkdown::render(input       = template,
                  output_file = out_html,
                  params      = list(input_dir = input))
## EOF
