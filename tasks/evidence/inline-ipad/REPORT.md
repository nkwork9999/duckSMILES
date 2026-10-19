# Fixed iPad compound panel — 2026-10-17

## Requested behavior

Keep search/categories and the library above a selected-compound display in the lower half, with properties and status visible without scrolling or a compound-detail modal.

## Implemented

- Normal iPad widths >=700 points: fixed upper/lower work surface.
- Compact import/search toolbar, All / Valid / Needs review filters, selected-row highlight.
- Paginated library; page capacity adapts to the available height. No scrolling container in this iPad surface.
- Lower panel: selected name/formula, per-input status, explicitly labelled connectivity preview, all 12 descriptors, original SMILES, local batch/cache status and elapsed time.
- Search and category filters update the selection; zero matches clear the detail. Graphs load only for the selected valid molecule and reuse the existing cache.
- Eight new interface strings added to all seven language files. Key/placeholder validation passes (60 keys per language).
- iPhone and accessibility text sizes retain their previous adaptive detail flow; fixed half-screen layout is scoped to the user's iPad screen.

## Verification

- Existing installed iPad app first failed the new behavior check: selecting Ethanol opened a modal and hid search. See `baseline-check.json`.
- Simulator compilation succeeded. iPad UI test passed, including inline selection, search remaining hittable, all descriptor labels within reach, invalid category warning and Caffeine selection after search.
- Landscape assertions passed: search and last descriptor remain hittable, and the lower panel is within the app window. Landscape screenshots still appear cropped/offset in the review viewer after waiting for window rotation and an additional second. Bounds/hit-target assertions pass, but clean landscape visual QA remains unresolved.
- Full iPad test suite: one iPad test passed, one iPhone-only flow skipped on iPad, zero failures. iPhone-only flow was not rerun in this iteration; chemistry/core source was not changed.
- Portrait screenshot: `ipad-portrait-final.png`; selected Caffeine screenshot: `ipad-selection.png`; landscape: `ipad-landscape.png`.

## Simulator recovery

The first test launch stalled before Chemob execution. A sampled Simulator installation proxy was in the dyld shared-cache fixup path. The task's iPad Simulator and CoreSimulator service were restarted without erasing device data; the unused task iPhone Simulator was shut down. Installation and tests subsequently succeeded. This environment delay is not a chemistry benchmark. The initial interrupted test attempt is preserved in the build/runtime notes.

## Remaining limits

The diagram is a connectivity preview, not a chemistry-aware 2D depiction. Translations remain drafts. No physical-device verification or new performance benchmark was performed. Very large accessibility sizes retain a scrolling layout rather than forcing all properties into half a screen.
