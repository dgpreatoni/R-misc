###############################################################################
# UAGRA R Scripts                                          PCA_utilities.R
###############################################################################
# PCA utilities - Commodity functions to analyse results of PCA
###############################################################################
#
# Version 1.0
#
# Description: functions to deal with prcomp/princomp outputs
#             could be adapted to deal with ade4 outputs...
#
# Usage: in your script, just source() this file: 
#        See function bodies for details.
#
# Requires:  
#
# Returns: 
#
###############################################################################
# created prea@uninsubria.it 20410210
# updated 
# revision history:
#   <author name> <timestamp> - <concise description of the changes made, one per line>
#
###############################################################################
# 
#    Copyright (C) 2014 Damiano G. Preatoni
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    A copy of the GNU General Public License is available at 
#    http://www.gnu.org/licenses/gpl.txt, or can be requested to the Free
#    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#    Boston, MA 02110-1301 USA.
#
###############################################################################

cat('Trying to load ggbiplot library...\n')
if(!require(ggbiplot)) {
  if(!require(devtools)) {install.packages("devtools", dep=TRUE)}
  install_github("ggbiplot", "vqv")
}


################################################################################
# makes a scree plot, looks better than ggbiplot::ggscreeplot
# http://marcoplebani.com/pca/
screePlot <- function(pcaobj, xlabs=NULL, ylabs=NULL) {
  # extract principal components variances
  vars <- data.frame(var=pcaobj$sdev**2, pc=1:length(pcaobj$sdev))
  # make a plot
  scree <- ggplot(vars, aes(pc, var)) + geom_bar(stat="identity", fill="gray") + geom_line() + xlab("PC Axis no.") + ylab('Explained Variance')
  print(scree)
  return(scree)
}

