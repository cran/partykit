
\name{NEWS}
\title{NEWS file for the \pkg{partykit} package}

\section{Changes in Version 1.2-5 (2019-07-17)}{
  \subsection{Bugfixes}{
  \itemize{
    \item{Trying to split in a variable where all observations were missing
          nevertheless produced a split, as reported by Kevin Ummel.}
  }
}
}

\section{Changes in Version 1.2-4 (2019-05-17)}{
  \subsection{Bugfixes}{
  \itemize{
    \item{update reference output, fix RNGversion}
  }
}
}


\section{Changes in Version 1.2-3 (2019-01-28)}{
 \subsection{New Features}{
    \itemize{
      \item{\code{varimp} runs in parallel mode, optionally.}
      \item{\code{weights} in \code{cforest} can now be used
            to specify a matrix of weights (number of observations
            times number of trees) to be used for tree induction (this was
            always possible in \code{party::cforest}. This was advertised in
            the documentation but actually not implemented so far.}
    }
  }
  \subsection{Bugfixes}{
  \itemize{
    \item{\code{predict} did not pay attention to \code{xlev}; this caused
          problems when empty factor levels were removed prior to tree fitting. See
         \url{https://stackoverflow.com/questions/51043870/lmtree-suspect-behavior-with-factors}.}
    \item{\code{nodeprune} may have got fitted terminal node numbers wrong,
          spotted by Jason Parker.}
  }
}
}


\section{Changes in Version 1.2-2 (2018-06-05)}{
  \subsection{Bugfixes}{
  \itemize{
    \item{In \code{mob()} using the \code{cluster} argument with a \code{factor}
      variable sometimes lead to \code{NA}s in the covariance matrix estimate
      if empty categories occured in subgroups. The problem had been introduced
      in version 1.2-0 and has been fixed now.}
    \item{Methods for the \code{sctest} generic from the \pkg{strucchange}
      package are now dynamically registered if \pkg{strucchange} is attached.
      Alternatively, the methods can be called directly using their full names
      \code{sctest.constparty} and \code{sctest.modelparty}.}
    \item{The \code{prune.modelparty} function is now fully exported but
      it is also registered with the \code{prune} generic from \pkg{rpart}.}
  }
  }
}

\section{Changes in Version 1.2-1 (2018-04-20)}{
  \subsection{New Features}{
    \itemize{
      \item{New \code{scale} argument for \code{predict} in \code{cforest}.
            For simple regression forests, predicting
            the conditional mean by nearest neighbor weights with
            \code{scale = TRUE} is now
            equivalent to the aggregation of means. The unscaled version
            proposed in <doi:10.1002/sim.1593> can be obtained with
            \code{scale = FALSE}.}
    }
  }
  \subsection{Bugfixes}{
  \itemize{
    \item{Bug fix for case weights in \code{mob()} in previous version (1.2-0)
      introduced a bug in the handling of proportionality weights. Both cases
      are handled correctly now.}
    \item{\code{glmtree} can now handle \code{caseweights = TRUE} correctly
      for \code{vcov} other than the default \code{"opg"}. Internally,
      the \code{glm} objects are adjusted by correcting the dispersion
      estimate and the degrees of freedom.}
    \item{\code{lookahead} did not work in the presence of missing values.}
    \item{Calling \code{partykit::ctree} did not work when partykit was not
      attached.}
    \item{\code{node_inner} now allows to set a different \code{gpar(fontsize = ...)}
      in the inner nodes compared to the overall tree.}
    \item{\code{splittest} asked for Monte-Carlo p-values, even
          when the test statistic was used as criterion.}
  }
  }
}

\section{Changes in Version 1.2-0 (2017-12-18)}{
  \subsection{New Features}{
    \itemize{
      \item{We welcome Heidi Seibold as co-author!}
      \item{Internal re-organisation for \code{ctree} by means of new extensible
            tree infrastructure (available in  \code{extree_data} and 
            \code{extree_fit}). Certain parts of the new infrastructure
            are still experimental. \code{ctree} is fully backward
            compatible.}
      \item{Use \pkg{libcoin} for computing linear test statistics 
            and p-values for \code{ctree}.}
      \item{Use \pkg{inum} for binning (the new \code{nmax} argument).}
      \item{Quadratic test statistics for splitpoint selection are now available for \code{ctree} via 
            \code{ctree_control(splitstat = "quadratic")}.}
      \item{Maximally selected test statistics for variable selection are now available for \code{ctree} via 
            \code{ctree_control(splittest = TRUE)}.}
      \item{Missing values can be treated as a separate category, also
            for splits in numeric variables in \code{ctree} via
            \code{ctree_control(MIA = TRUE)}.}
      \item{Permutation variable importance, including conditional variable
            importance, was added to \pkg{partykit}.}
      \item{New \code{offset} argument in \code{ctree}.}
      \item{New \code{get_paths} for computing paths to nodes.}
      \item{\code{node_barplot} gained a \code{text} argument
            that can be used to draw text labels for the
            percentages displayed.}
      \item{The \code{margins} used in \code{plot.party} can now
            also be set by the user.}
    }
  }
  \subsection{Bugfixes}{
  \itemize{
    \item{Bug fix in \code{mob()} if \code{weights} are used and \code{caseweights = TRUE} (the default).
          The statistics for the parameter instability tests were computed incorrectly and
	  consequently the selection of splitting variables and also the stopping criterion
	  were affected/incorrect.} 
    \item{Avoid log(p) values of \code{-Inf} inside \code{mob()} by replacing weighted averaging
          with naive averaging in the response surface regression output in case
	  the p values are below machine precision.}
    \item{The \code{as.party} method for \code{rpart} objects
          without any splits only returned a naked \code{partynode}
	  rather than a full \code{party}. This has been corrected
	  now.}
    \item{\code{nodeapply} did not produce the same results for permutations of \code{ids}.
          Spotted by Heidi Seibold.}
    \item{Out-of-bag predictions in \code{predict.cforest} were incorrect.}
    \item{\code{perm} in \code{predict} was only considered when
          \code{newdata} was given. Spotted by Heidi Seibold.}
    \item{Don't try to search for binary splits in unordered factors with more than 31
          levels. This potentially caused an integer overrun in previous versions.
          \code{party::ctree()} uses an approximation for binary split
          searches in unordered factors; thus, using \pkg{party} might be an
          alternative.}
  }
  }
}

\section{Changes in Version 1.1-1 (2016-09-20)}{
  \itemize{
    \item{Proper support of quasi-families in \code{glmtree}
          and hence \code{palmtree}.}
    \item{NA handling by following the majority was 
          potentially incorrect in \code{ctree}.}
    \item{Minor speed improvements.}
    \item{Breaking ties before variable selection was
          suboptimal for very small log-p-values.}
  }
}

\section{Changes in Version 1.1-0 (2016-07-14)}{
  \itemize{
    \item{Added a new function \code{palmtree} that fits partially
      additive (generalized) linear model trees. These employ
      model-based recursive partitioning (\code{mob}) based on (generalized)
      linear models with some local (i.e., leaf-specific) and
      some global (i.e., constant throughout the tree) regression
      coefficients.}
    \item{Splits in ordinal variables are now represented correctly 
      in the (still internal) \code{list.rule} method.}
    \item{Kaplan-Meier curves in \code{"constparty"} trees were
      plotted incorrectly due to use a wrong scaling of the
      x-axis. Spotted by Peter Calhoun <calhoun.peter@gmail.com>.}
    \item{Use \code{quote(stats::model.frame)} instead of
      \code{as.name("model.frame")}.}
    \item{The \code{as.party} methods for \code{rpart} and
      \code{Weka_tree} now have a \code{data = TRUE} argument
      so that by default the data is preserved in the \code{party}
      object (instead of an empty model frame).}
    \item{The \code{predict} method for \code{cforest} objects
      did not work for one-row data frames, fixed now.}
    \item{Added \code{rot} and \code{just} arguments to \code{node_barplot}
      for more fine control of x-axis labeling (e.g., with 45 degree
      rotation).}
  }
}

\section{Changes in Version 1.0-5 (2016-02-05)}{
  \itemize{
    \item{The \pkg{partykit} package has now been published in
      \emph{Journal of Machine Learning Research}, 16, 3905-3909.
      \url{http://jmlr.org/papers/v16/hothorn15a.html}}
    \item{Added support for setting background in panel functions.}
    \item{The \code{as.list()} method for \code{partynode} objects
      erroneously created an object \code{thisnode} in the calling
      environment which is avoided now.}
  }
}

\section{Changes in Version 1.0-4 (2015-09-29)}{
  \itemize{
    \item{Bug fix in \code{plot()} method for \code{constparty} objects.
          In the previous partykit version clipping was accidentally 
	  also applied to the axes labels.}
    \item{For \code{constparty} objects \code{plot(..., type = "simple")}
          did not work correctly whereas \code{plot(as.simpleparty(...))}
	  yielded the desired visualization. Now internally \code{as.simpleparty()}
	  is called also in the former case.}
    \item{The \code{as.simpleparty()} method now preserves p-values from
          \code{constparty} objects (if any).}
    \item{Added a \code{getCall()} method for \code{"party"} objects.}
    \item{In the \code{predict()} method for \code{"lmtree"} and \code{"glmtree"}
          objects the \code{offset} (if any) was sometimes ignored. It is
	  now always used in the prediction.}
  }
}

\section{Changes in Version 1.0-3 (2015-08-14)}{
  \itemize{
    \item{Import \code{logrank_trafo} from \pkg{coin}.}
  }
}

\section{Changes in Version 1.0-2 (2015-07-28)}{
  \itemize{
    \item{\code{nodeprune(..., ids = 1)} did not prune the tree to the root
      node. Fixed now.}
    \item{\code{predict.cforest} used \code{na.omit} instead of \code{na.pass}.}
    \item{\code{predict.party} now features new \code{perm} argument for
          permuting splits in specific variables (useful for computing
          permutation variable importances).}
   \item{\code{NAMESPACE} updates.}
  }
}

\section{Changes in Version 1.0-1 (2015-04-07)}{
  \itemize{
    \item{The support for (generalized) linear model trees with just
      a constant regressor has been improved. Now \code{lmtree(y ~ x1 + x2)}
      is short for \code{lmtree(y ~ 1 | x1 + x2)}, analogously for \code{glmtree()}.
      Plotting now also works properly in this case.}

    \item{The \code{as.party()} method for \code{"rpart"} objects did not
      work if one of the partitioning variables was a \code{"character"}
      variable rather than a \code{"factor"}. A suitable work-around has
      been added.}
      
    \item{The \code{node_barplot()} panel function can now also be used
      for multivariate responses, e.g., when all responses are numeric and
      on the same scale.}

    \item{The package now also includes a new data set \code{HuntingSpiders}
      which is essentially a copy of the \code{spider} data from the package
      \pkg{mvpart} that is currently archived on CRAN. The documentation has
      been improved somewhat and is likely to change further to explain how
      the data has been transformed in De'ath (2002).}

    \item{The survival tree example for the GBSG2 data was broken due to
      the response being (incorrectly) also part of the explanatory variables.
      Fixed by using the latest \pkg{Formula} package (at least version 1.2-1).}
  }
}

\section{Changes in Version 1.0-0 (2015-02-20)}{
  \itemize{
    \item{Version 1.0-0 published. This version is described in the MLOSS paper 
          accepted for publication by the Journal of Machine Learning Research today.}
    \item{The unbiased version of \code{cforest} (with \code{replace = FALSE}) is 
          now the default (as in \pkg{party}).}
    \item{Register all S3 methods in \code{NAMESPACE}.}
  }
}

\section{Changes in Version 0.8-4 (2015-01-06)}{
  \itemize{
    \item{Extended \code{mob()} interface by a \code{cluster} argument. This can be
    a vector (numeric, integer, factor) with cluster IDs that are
    then passed on to the 'fit' function (if supported) and used
    for clustering the covariance matrix in the parameter stability
    tests. \code{lmtree()} and \code{glmtree()} hence both gained a \code{cluster}
    argument which is used only for cluster covariances but not
    for the model estimation (i.e., corresponding to a working
    independence model).}

    \item{Optionally, the parameters' variance-covariance matrix in \code{mob()}
    can now be estimated by the sandwich matrix instead of the default
    outer-product-of-gradients (OPG) matrix or the information matrix.}

    \item{Reimplementation of \code{cforest()} available with extended 
    prediction facilities. Both the internal representation and the user interface
    are still under development are likely to change in future versions.}
    
    \item{Added multicore support to \code{mob()}, \code{ctree()}, and \code{cforest()}.
    If control argument \code{cores} is specified (e.g., \code{cores = 4}) then the
    search for the best variable or split point (often involving numerous model fits in
    \code{mob()} or resampling in \code{ctree()}) is carried out using \code{parallel::mclapply()}
    rathern than sequential \code{for()} or \code{sapply()}. Additionally, other
    \code{applyfun}s can be provided, e.g., using networks of workstations etc.}

    \item{Bug fix in \code{mob()} that occurred when regressor variables and
    partitioning variables overlapped and were not sorted in the
    underlying model frame.}
  }
}

\section{Changes in Version 0.8-3 (2014-12-15)}{
  \itemize{
    \item{\pkg{mvpart} was archived 2014-12-15.}
  }
}

\section{Changes in Version 0.8-2 (2014-09-12)}{
  \itemize{
    \item{Fixed an uninitialized memory issue reported by valgrind.}
  }
}

\section{Changes in Version 0.8-1 (2014-09-08)}{
  \itemize{

  \item{partykit now depends on R version >= 3.1.0 in order to import the
        \code{depth()} generic from the \pkg{grid} package.}

  \item{The print methods for \code{party}/\code{partynode} objects with only a root node
    was modified. Now, the terminal panel function is also applied
    if there is only a root node (while previously it was not).}

  \item{ \code{ctree()} now catches \code{sum(weights) <= 1} situations before they 
    lead to an error.}

  \item{ Code from suggested packages is included by using \code{::} syntax as
    required by recent R versions.}
    
  \item{ Argument \code{ytrafo} of \code{ctree()} can now be a function which will be
    updated in every node.}
    
  \item{ A small demo briefly illustrating some memory and speed properties
    has been added. It can be run interactively via
    \code{demo("memory-speed", package = "partykit").}}
    
  \item{ Section 3 of the "constparty" vignette now shows how properties of
    a new tree algorithm can be assessed by using \pkg{partykit} building
    blocks.}
  }
}


\section{Changes in Version 0.8-0 (2014-03-27)}{
  \itemize{

  \item{Major improved version of \pkg{partykit}. The previously existing functions
    in the package were tested and enhanced, new functions and
    extensive vignettes added.}

  \item{Extended and improved introductory documentation. The basic classes
    and class constructors \code{partynode}/\code{partysplit}/\code{party} are introduced in 
    much more detail now in \code{vignette("partykit", package = "partykit")}.}
    
  \item{The class \code{constparty} (inheriting from \code{party}) for representing \code{party}
    objects with constant fits in the nodes (along with coercion methods
    for \code{rpart}, \code{J48}, etc.) is now described in more detail in the new
    \code{vignette("constparty", package = "partykit")}.}

  \item{The package now includes a reimplementation of the model-based
    recursive partitioning algorithm (MOB) using \pkg{partykit} infrastructure.
    The generic algorithm, the corresponding convenience interfaces
    \code{lmtree()} and \code{glmtree()} as well as various illustrations and possible
    extensions are described in detail in the new
    \code{vignette("mob", package = "partykit")}.}

  \item{Improved implementation of conditional inference trees (CTree), see
    the new \code{vignette("ctree", package = "partykit")} for details.}

  \item{New \code{nodeprune()} generic for pruning nodes in all \code{party} trees and
    \code{partynode} objects.}

  \item{Deal with empty levels in \code{ctree()} for \code{teststat = "quad"}
    (bug reported by Wei-Yin Loh <loh_at_stat.wisc.edu>).}

  \item{In \code{predict()} method for \code{constparty} objects, \code{type = "prob"} now returns
    ECDF for numeric responses and \code{type = "response"} the (weighted) mean.}

  \item{New panel function \code{node_ecdf()} for plotting empirical cumulative
    distribution functions in the terminal leaves of \code{constparty} trees.}
  }
}

\section{Changes in Version 0.1-6 (2013-09-03)}{
  \itemize{

  \item{Bug fix in \code{as.party()} method for J48 trees with ordered factors.}
 }
}

\section{Changes in Version 0.1-5 (2013-03-22)}{
  \itemize{

  \item{Fix C code problems reported by clang under OS X.}
  }
}

\section{Changes in Version 0.1-4 (2012-06-05)}{
  \itemize{

  \item{Added \code{node_surv()} for plotting survival ctrees. Accompanying
    infrastructure for survival trees was enhanced.}
    
  \item{\code{ctree()} now checks for (and does not allow) \code{x >= max(x)} splits.}
  }
}

\section{Changes in Version 0.1-3 (2012-01-11)}{
  \itemize{

  \item{Added \pkg{ipred} to the list of suggested packages due to usage of
        GlaucomaM and GBSG2 data in tests/examples.}
  }
}


\section{Changes in Version 0.1-2 (2011-12-18)}{
  \itemize{

  \item{The \code{node_terminal()} panel-generating function is now customizable
    by a FUN argument that is passed to \code{formatinfo()}.}

  \item{The \code{plot()} method for \code{simpleparty} object now sets up a formatting
    function passed to \code{formatinfo()}, both in \code{print()} and \code{plot()}.}

  \item{Fixed bug in \code{pmmlTreeModel()} for processing label IDS in splits when
    not all levels are present.}
  
  \item{Cleaned up unused variables in C code and partial argument matching
    in R code.}
  }
}  

\section{Changes in Version 0.1-1 (2011-09-29)}{
  \itemize{

  \item{First CRAN release.}
  
  \item{See \code{vignette("partykit", package = "partykit")} for a (somewhat rough)
    introduction to the package and its classes/methods.}
  }
}
