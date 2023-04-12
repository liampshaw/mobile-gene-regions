'Combine two plots together.

Usage:
  combine_two_plots.R <plot_one> <plot_two> [--output_pdf=<pdf>] [--width=<w>] [--height=<h>] [--rel_width=<rw>]

Options:
  -h --help     Show this screen.
  --version     Show version.
  --output_pdf=<pdf> Name of output pdf. [default: plot_output.pdf]
  --width=<w> Width of output pdf. [default: 8]
  --height=<h> Height of output pdf. [default: 6]
  --rel_width=<rw> Relative width of plot one compared to plot two. [default: 1]

' -> doc

library(docopt)
library(ggplot2)
library(cowplot)

args <- docopt(doc, version = 'Combine_two_plots v1.0')

# Get file types (extensions)
file_type_one = gsub(".*\\.", "", args$plot_one)
file_type_two = gsub(".*\\.", "", args$plot_two)


if (file_type_one=="pdf"){
  p.left <- cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(args$plot_one, density = 600))
}
if (file_type_one=="png"){
    p.left <- cowplot::ggdraw() + cowplot::draw_image(args$plot_one)
}

if (file_type_two=="pdf"){
  p.right <- cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(args$plot_two, density = 600))
}
if (file_type_two=="png"){
    p.right <- cowplot::ggdraw() + cowplot::draw_image(args$plot_two)
}

# Combine
p.final = cowplot::plot_grid(p.left, p.right, 
  nrow=1, align='h',
  rel_widths=as.numeric(args$rel_width))
ggsave(p.final, file=args$output_pdf, width=as.numeric(args$width), height=as.numeric(args$height))

