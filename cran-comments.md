## Submission

This is a bug fix release. It follows 1.1-0 more closely than I would normally
submit an update, for the reason set out below. If you would prefer that I hold
it until a more usual interval has passed, please say so and I will do that.

Version 1.1-0 fixed the vignette failure that led to the archiving of 1.0-8.
While adding a regression test suite to the package immediately afterwards, I
found several latent defects that prevent documented functionality from being
used at all:

* `pushCross()` evaluated its `type` argument before `match.arg()` had resolved
  it, so a condition of length four was tested. Since R 4.2.0 that is an error
  rather than a warning, and the function therefore fails whenever the default
  argument is used.

* `pullCross()` indexed with the whole default `type` vector, which performs
  recursive indexing rather than selecting a component, and fails.

* Two `warning()` calls in `mstmap.data.frame()` contained unescaped quotes.
  R parses these as an arithmetic expression on two character strings, so they
  fail at run time rather than issuing the documented warning. The same code
  also assigned corrected genotype names to the marker names. A marker or
  genotype name containing a space therefore caused an error.

Because these defects are in code paths that already fail, no successful
computation is affected by correcting them.

Two memory leaks in the C++ sources are also fixed. A working array in the
marker ordering routine was allocated on every call and never released, and the
object holding the genotype data and the n by n pairwise distance matrix was
never released at the end of a map construction. I would prefer these were not
present when the package is next checked with valgrind.

The linkage maps produced by the package are unchanged. The regression tests
added in this version pin the marker orders and genetic distances for the
bundled datasets, and they are unchanged. For the largest bundled dataset I
also compared the output of the previous and the current build of the C++ code
directly, and the marker orders and genetic distances are identical.

## Test environments
* Local Windows 11 install, R 4.4.3
* Github Actions:
    - macOS: r-release
    - windows: r-devel, r-release, r-oldrel
    - ubuntu 24.04: r-release, r-oldrel

## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTEs:

```
* checking CRAN incoming feasibility ... NOTE
Days since last update: 1
```

This is the short interval explained above.

```
* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
             user system elapsed
mstmap.cross 0.67   1.45    7.45
```

Only the elapsed time exceeds the threshold. It is dominated by file I/O, as
the underlying MSTmap algorithm reads and writes temporary files. Combined user
and system CPU time remains well below 5 seconds.

The tests added in this version take approximately 18 seconds. The tests that
compare stored marker orders and genetic distances are skipped on CRAN, since
those are floating point results that may differ in the last bits between
platforms and compilers; they are run locally and in continuous integration.
