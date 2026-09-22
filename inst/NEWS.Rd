\name{news}
\title{Release information for ASMap}
\section{Changes in version 1.1-0}{
\subsection{README}{
  \itemize{
    \item This release restores \pkg{ASMap} to CRAN following the
    archiving of version 1.0-8. There are no user visible changes to the
    functionality of the package.
  }
}
\subsection{BUG FIXES}{
  \itemize{
    \item The package vignette used the \code{tabu} LaTeX package, which
    is no longer maintained and is incompatible with current versions of
    the \code{array} package. Re-building the vignette therefore failed
    with an \dQuote{Undefined control sequence} error, which led to the
    archiving of version 1.0-8 on CRAN. The vignette now uses
    \code{tabularx} from the standard LaTeX \code{tools} bundle, which
    provides the same functionality.
    \item \code{pushCross} and \code{pullCross} used the \code{type}
    argument before it had been resolved by \code{match.arg}. In
    \code{pushCross} this tested a condition of length four, which R
    4.2.0 and later treat as an error, so the function failed whenever
    the default argument was used. In \code{pullCross} it caused
    recursive indexing. Both now resolve \code{type} before first use.
    \item \code{pushCross} checked the number of genotypes in the pushed
    marker set before establishing that any such markers existed, which
    produced an uninformative error. It now reports that there are no
    markers of that type to push back.
    \item \code{mstmap.data.frame} failed rather than warning when marker
    or genotype names contained spaces, and assigned corrected genotype
    names to the marker names. Spaces are now replaced with a
    \code{"-"} separator in the correct component, with a warning.
    \item \code{alignCross} raised an uninformative error when the
    \code{maps} argument was not a list, and when \code{object} had a
    single class. Both are now reported correctly.
    \item \code{mstmap.data.frame} did not reject numeric input for
    \code{RILn} population types, because the test used
    \code{is.numeric} on a \code{data.frame}.
    \item \code{mstmap.cross} and \code{mstmap.data.frame} could pick up
    an unrelated object named \code{omit.list} from the calling
    environment.
    \item Corrected the axis labels produced by \code{pValue} and a
    cross reference in a \code{statMark} error message.
  }
}
\subsection{INTERNAL CHANGES}{
  \itemize{
    \item Two memory leaks in the MSTmap C++ code have been repaired. A
    working array in the marker relocation step was allocated on every
    call and never released, and the object holding the genotype data
    and the pairwise distance matrix was never released at the end of a
    map construction. For a linkage group of \eqn{n} markers the latter
    retained an \eqn{n \times n} matrix of doubles per call. The base
    class destructor is now virtual, as the object is owned through a
    base class pointer.
    \item The distance function member is now initialised, and an
    unrecognised distance function is reported rather than leading to
    the use and release of an uninitialised pointer.
    \item Corrected undefined behaviour when testing whether a genotype
    string is numeric, which passed a possibly negative \code{char} to
    \code{isdigit}.
    \item Added a guard, equivalent to the one already used for doubled
    haploid populations, for the case where two markers in a RIL
    population share no scored individuals. This previously divided by
    zero and propagated \code{NaN} into the distance matrix.
    \item These changes do not alter the linkage maps produced. Marker
    orders and genetic distances are unchanged on all bundled datasets.
  }
}
}
\section{Changes in version 1.0-x}{
\subsection{README}{
  \itemize{
    \item Version 1.0-x of the \pkg{ASMap} package celebrates the
    Journal of Statistical Software publication, Taylor and Butler
    (2017). The paper succinctly summarizes the important contents of
    the package vignette and provides a useful quick reference for the
    functionality of the package. If you use \pkg{ASMap} in your
    publication please reference it using the bibliographic reference obtainable using
    \code{citation("ASMap")} and given below.
    \item Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic
    Linkage Map Construction and Diagnosis. Journal of Statistical Software,
    \bold{79}(6), 1--29.
  }
}
\subsection{NEW FEATURES}{
    \itemize{
      \item There is now a p-value function to emulate Figure 1.1 in
      the vignette. The function allows users to specify their own set
      of threshold distances and population sizes, making it easier to
      determine the correct choice of \code{p.value} to use in
      \code{mstmap.cross} and \code{mstmap.data.frame}.
      \item The function \code{combineMap} has been updated to include
      merging of linkage maps by \code{"marker"} or
      \code{"genotype"}. Merging by \code{"genotype"} is identical to
      the previous version of the function. Merging by
      \code{"marker"} provides the ability to combine linkage maps with
     different sets of genotypes through common markers.
    }
  }
\subsection{BUG FIXES}{
  \itemize{
    \item \code{fixClones} muddled linkage groups when an inputted map
    had more than one linkage group. This has now been amended.
    \item The \code{fixClones} function crashed when \code{consensus =
      FALSE} was set. This has now been amended.
    \item The \code{cross.scheme} of the object was lost when pushing
    back markers using \code{pushCross}. This has now been amended.
    \item The function \code{quickEst} crashed when there were linkage
    groups with one marker. This has now been amended.
    \item The function \code{alignCross} crashed with a non-descript
    error when there were no matching markers in supplied reference
    maps. This has now been amended to provide a more descript error.
   }
 }
 }
\section{Changes in version 0.4-x}{
  \subsection{README}{
  \itemize{
  \item This is the first documented NEWS release of ASMap. For this
  reason it contains information about features of the package as well
  as recent bug fixes.
 }
}
\subsection{NEW FEATURES}{
    \itemize{
      \item The package now contains a function to assist in identifying
      and orienting linkage groups (see \code{?alignCross}).
      \item The heat map in the package now uses a diverging colour
      palette chosen from the default palettes available in
      \pkg{RColorBrewer}.
      \item The package now includes a PDF vignette. This can be found
      on any operating system using the command

      \code{> vignette("ASMap")}

      The vignette includes extensive descriptions of the package
      functions as well as fully worked example that showcases the
      package functions throuh pre-diagnosis, construction and post
      diagnosis of a linkage map. It also describes how the package can
      be used efficiently for techniques such as fine mapping or combining
      maps.
      \item The package now includes three fully documented marker data
      sets for users to explore. See the vignette for more details.
}
}
\subsection{BUG FIXES}{
  \itemize{
    \item Fixed a nasty bug in \code{combineMap()} that switched
      genotype names for certain linkage maps
     \item Fixed bug that returned inflated distances for \code{"ARIL"}
     populations when using \code{mstmap.data.frame()}.
     \item If the constructed map contained linkage groups with one
      marker \code{profileMark()} graphically displaced the interval
      statistics. This has now been amended.
     \item \code{pushCross()} can now push large volumes of co-located
      markers back into the constructed map very efficiently.
    \item \code{profileMark()} was not outputting \code{crti.val} column
      when \code{crit.val = "bonf"} argument was chosen. This has now
      been amended.
      \item \code{quickEst()} has been improved and now provides an even
      more accurate and efficient distance calculations for linkage maps.
    }
}
}


