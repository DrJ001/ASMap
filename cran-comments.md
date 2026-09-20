## Submission

This submission restores ASMap to CRAN following its archival on 2026-06-29.

The archival resulted from a WARNING in the "re-building of vignette outputs"
check on the r-devel Fedora flavors:

```
! Undefined control sequence.
\tabu@cleanup ...bu@naturalX =\tabu@naturalX@save
l.109 \end{tabu}
```

The cause was the vignette's use of the `tabu` LaTeX package, which is no
longer maintained and is incompatible with current versions of the `array`
package. `tabu` was used in a single location (a title page layout) and has
been replaced with `tabularx`, which is part of the standard `tools` bundle
and provides the same `X` column specifier. The vignette now re-builds
cleanly.

The vignette PDF is also now compacted with `--compact-vignettes="gs+qpdf"`.

## Test environments
* Local Windows 10 install, R 4.4.3, TeX Live 2023
* Github Actions:
    - macOS: r-release
    - windows: r-devel, r-release, r-oldrel
    - ubuntu 24.04: r-release, r-oldrel

## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTEs:

```
* checking CRAN incoming feasibility ... NOTE
New submission
Package was archived on CRAN
```

This is expected, as the package is being resubmitted following archival.

```
* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
             user system elapsed
mstmap.cross 0.69   1.58    5.86
```

The elapsed time for this example is dominated by file I/O, as the underlying
MSTmap algorithm reads and writes temporary files. Combined user and system
CPU time is under 2.5 seconds.
