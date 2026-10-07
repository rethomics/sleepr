## Release summary

sleepr 0.4.0. `sleep_annotation()` now needs a sleep rule (`"classic"`, the previous
behaviour, or the new `"k"` rule), so existing calls must add `rule = "classic"` or
declare it once with `options(sleepr.sleep_rule = "classic")`. See NEWS.md.

## Maintainer

The current maintainer (cre) is Quentin Geissmann. This submission needs his approval,
or his email confirming a change of maintainer.

## Test environments

* local: Manjaro Linux, R (release)

## R CMD check results

0 errors | 0 warnings | 1 note

* "Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed":
  local environment only.

## Reverse dependencies

scopr and ggetho (rethomics) do not call `sleep_annotation()` themselves; users pass it
to `scopr::load_ethoscope(FUN = ...)`, which forwards extra arguments such as `rule`.
