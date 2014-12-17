##  THIS NEEDS TO BE INITIALIZED WITH THE SCRIPT TEMPLATE USED FOR EVERY R SCRIPT
##
##  Copyright (c) 2013, 2014, Boston University. All rights reserved.
##  
##  Redistribution and use in source and binary forms, with or without
##  modification, are permitted provided that the following conditions are met: 
##  
##  1. Redistributions of source code must retain the above copyright notice, this
##     list of conditions and the following disclaimer. 
##  2. Redistributions in binary form must reproduce the above copyright notice,
##     this list of conditions and the following disclaimer in the documentation
##     and/or other materials provided with the distribution. 
##  
##  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
##  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
##  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
##  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
##  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
##  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
##  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
##  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
##  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
##  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##  
##  The views and conclusions contained in the software and documentation are those
##  of the authors and should not be interpreted as representing official policies, 
##  either expressed or implied, of Boston University.
##  
##  Authors:
##    name1 [1], name2 [1], name3 [2],
##    Stefano Monti [1,2]
##  
##  [1] Bioinformatics Program, Boston University
##  [2] Center for Computational Biomedicine, Boston University  

if ( !exists("CBMDEV") ) {
  CBMDEV <- Sys.getenv('CBMDEV')
  if (CBMDEV=="") stop( "Use 'setenv CBMDEV ..' to set CBMrepository's base directory" )
}
if ( !exists("CBMMLAB") ) {
  CBMMLAB <- Sys.getenv('CBMMLAB')
  if (CBMMLAB=="") stop( "Use 'setenv CBMMLAB ..' to set CBMrepository's base directory" )
}
source( paste(CBMMLAB, "R/misc.R", sep="/") )
#source( paste(RHOME, "", sep="/") )
#require(Biobase)

