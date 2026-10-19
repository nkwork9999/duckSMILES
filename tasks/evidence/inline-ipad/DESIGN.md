# Inline iPad workbench

User-requested layout: keep compound information and status below search/category controls, in the lower half of the display, without scrolling or a compound-detail modal.

- Normal iPad widths (>=700 pt): compact brand/actions, search, category selector, paginated table in upper half; selected compound panel fixed below.
- Categories: All / Valid / Needs review. Search combines with category. Page controls keep large libraries accessible without a scrolling table.
- Detail: selected input name/formula, result status, labelled connectivity preview, all 12 descriptors, original SMILES, global calculation/reuse status.
- Selection follows filtered results; a no-match query clears the detail rather than leaving a stale molecule visible.
- iPhone and accessibility text-size layouts retain their existing adaptive detail flow in this iteration; the fixed two-half layout is scoped to the user's iPad screen.
- Connectivity preview remains explicitly distinct from a chemical 2D depiction. No chemistry or cache semantics change.

Verification contract: selecting a row leaves search/category controls hittable and shows inline properties; selecting Needs review shows the invalid row's original-input warning in the lower panel; searching Caffeine excludes Ethanol and updates the selected details. Verify portrait and landscape captures and retain the existing iPhone flow.
