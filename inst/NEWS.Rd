\name{NEWS}
\title{News for Package 'sglOptim'}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in sglOptim version 1.2.0 (2015-09-16)}{
  \itemize{
  \item Improved flexibility and reliability :
  \itemize{
    \item added handling of q = 1 cases (each feature has one parameter)
    \item correct handling of non penalized parameters when computing lambda.max (fit non penalized parameters before computing critical lambdas)
    \item increased speed when optimizing with a high number of groups/classes (> 100) 
    \item some improvements on error reporting and argument validation
   }
   
   \item Fixes :
     \itemize{
        \item fixed bug related to grouping 
        \item fixed seval minor bugs
        \item fixed startup messages 
     }
  }
}