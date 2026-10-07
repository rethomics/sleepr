# sleepr 0.4.0

**Breaking: `sleep_annotation()` needs a sleep rule.** To reproduce earlier results,
add `rule = "classic"`, or declare it once with `options(sleepr.sleep_rule = "classic")`
(or `SLEEPR_SLEEP_RULE=classic` in `.Renviron`). With scopr, extra arguments reach
`FUN`: `load_ethoscope(metadata, FUN = sleep_annotation, rule = "classic")`. Without a
rule, `sleep_annotation()` stops with an explanation.

* `rule = "k"` (k = 3 by default; `"k2"` is stricter): sleep scored from walking and
  sustained movement events, ignoring tracking noise. It matches video ground truth at
  night, and flies it scores asleep respond to air puffs like sleeping flies. The
  algorithm and defaults are identical to ethoscopy 3.0.0 (`sleep_rules`).
* `untracked = "immobile"` (default, as before) or `"break"` in `sleep_annotation()`.
* `motion_qc()`: per animal and light phase, how often a still animal crosses the
  classic threshold, and how much of the recording is untracked; with the helpers
  `find_still_bins()`, `find_spikes()`, `light_phase()`, `pixel_size()` and
  `default_still_shift()`.
* `classify_events()` and `k_rule_bins()` are exported.
