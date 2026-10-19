# Chemob interface and localization — 2026-10-17

## Changes

- Preserved paper/ink/lime identity; lime now marks the primary Add SMILES action.
- Replaced individual oversized cards and identical hexagon decorations with separated, compact compound rows.
- iPhone rows show compound name, formula, MW, cLogP and TPSA. Invalid input has a textual warning and icon.
- At widths >=700 points the list aligns MW, cLogP, TPSA and HBD/HBA into columns.
- Search and import remain above the scrolling normal-size library. Search has a clear action and an empty-results state.
- Accessibility text sizes use a fully scrolling vertical header/actions and expanded compound rows. A maximum-size inspection exposed clipping in the first version; the final implementation removes that fixed header at accessibility sizes.
- Native button semantics, stable UI-test identifiers, press feedback and original input data are retained.
- English source UI and errors, with Japanese, Indonesian, Hindi, German, French and Italian translations. Device preferred language selects the interface. English is the development/fallback language. Shared `.strings` files require no network or translation SDK.
- Machine-readable CSV identifiers, SMILES and user compound names retain their original text.

## Verification

- iPhone 17e / iOS 26.5 Simulator: UI test passes, including sample loading, six sample row hit targets, search excluding unrelated rows, clear search, descriptor detail, invalid input detail, paste input and cached-result reuse.
- Swift Core: 10 tests passed, 0 failures.
- Translation validator: all seven files have 52 identical keys, nonempty translations, unique keys and matching interpolation placeholders.
- Visible screenshots inspected: English phone; Swiss German phone; Hindi phone; English iPad Air 11-inch portrait and landscape; maximum accessibility text on phone.
- Simulator build succeeded with Swift 6. The AppIntents metadata extraction notice is present because the app has no AppIntents dependency; no new Swift compiler warning was found.

## Limits

- Translation wording is a first draft; native-speaker review is still required before release. Other Indian languages and Romansh are not included.
- This is Simulator verification; physical-device rendering/performance and 375-point width have not been verified in this iteration. Light appearance remains intentional; full dark styling is not implemented.
- The detail view still has an explicitly labelled connectivity preview. True chemical 2D coordinates/thumbnails, sorting and batch selection remain future work.
- Screenshot timing badges are individual warm-cache observations, not a new performance comparison.

## Evidence

- `iphone-en.png`: current compact English library.
- `iphone-de-CH.png`, `iphone-hi.png`: translated library controls.
- `ipad-en.png`, `ipad-landscape.png`: aligned property columns; Simulator raw landscape capture stores rotated pixel orientation.
- `iphone-large-text.png`: expanded, scrolling maximum accessibility layout.
- `localization-check.json`: translation-key and placeholder validation.
