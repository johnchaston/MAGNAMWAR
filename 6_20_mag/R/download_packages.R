# #' Download Required Packages #' #' Automatically downloads all the required packages for full analysis #' @export
# download_packages <- function() { try(install.packages('lme4', dependencies=T), F) try(install.packages('multcomp',
# dependencies=T), F) try(install.packages('seqinr', dependencies=T), F) try(install.packages('qpcR', dependencies=T),
# F) } library('devtools') use_package('lme4') use_package('multcomp') use_package('seqinr') use_package('plyr')
# use_package('qqman') use_package('ape') use_package('ggplot2')
