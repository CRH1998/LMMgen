% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_kinship_data.R
\name{pair_parents_with_children}
\alias{pair_parents_with_children}
\title{pair_parents_with_children}
\usage{
pair_parents_with_children(parents, unassigned_individuals)
}
\arguments{
\item{parents}{Data frame. Data frame with four columns containing dadid, sex, momid, sex}

\item{unassigned_individuals}{Data frame. Data frame containing id and sex}
}
\value{
list containing the current family dataframe (sub_family_df) and the remaining unassigned individuals (unassigned_individuals)
}
\description{
pairs parents with unassigned individuals
}
