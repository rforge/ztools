\name{flush_output}
\alias{flush_output}

\title{Outputs and flushes the given text}

\description{
 Outputs and flushes the given text so that it stays in one line.
}
\usage{
flush_output(dynamic.text,fixed.text)
}

\arguments{
     \item{dynamic.text}{ Any text or number that changes dynamically. Typically an index in a loop.}
     \item{fixed.text}{ dynamic.text and fixed.text are pasted together in one line. Typically contains the end of a loop index.}
}

\details{
The function handles the often annoying situation of providing an output in a loop where the information is constantly being shifted out of view with each step. The inputs \code{dynamic.text} and \code{fixed.text} are pasted together, where the first ist constantly changing and the second stays...well, fixed. The output is flushed using \code{\link{flush.console}} with each step, so that the display of output in the console is current. The pasting action relies on \code{\link{paste0}} and uses no separator, so you have to provide the approriate text.

}

\value{
Nothing. 
}

\examples{
n <- 10000
for(i in 1:n){
flush_output(i,paste0(" of ",n))
}
}