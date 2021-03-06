% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{metadataMetrics_violin}
\alias{metadataMetrics_violin}
\title{Plot the violin plot of the distribution per group (y) of a desired metric in the metadata}
\usage{
metadataMetrics_violin(object, x, y, output_image, width = 17, height = 15)
}
\arguments{
\item{object}{Seurat object}

\item{x}{Groups to divide the violin plot. Violing use this color.}

\item{y}{metric in metadata to plot distribution}

\item{output_image}{file name of the image}

\item{width}{image width}

\item{height}{image height}
}
\value{
PNG file. violin plots.
}
\description{
Plot the violin plot of the distribution per group (y) of a desired metric in the metadata
}
