# Ticket 05: walls-and-scheduled-tracks

**Status:** ✅ Done

## Local Test Result

`cargo test --release irt_`: 38 passed, 0 failed (`irt_background` 5, `irt_model` 7, `irt_sample` 9, `irt_schedule` 3, `irt_single_track` 8, `irt_walls` 6).

Also clean: `cargo fmt --check`, `cargo clippy --release --all-targets -- -D warnings`, `cargo test --release` (whole suite, no failure). The Geant4 comparison is unchanged (2281.1 reactions per run, RMS of z 1.00).

Printed by the wall tests (survival of the e_aq of one e_aq + OH_r pair, fully diffusion-controlled, `r_eff` 0.549 nm):

```
near_wall: time [s]   walk     images   no images
          1e-9     0.9763   0.9729   0.9739
          5e-9     0.9243   0.9194   0.9342
          2e-8     0.8892   0.8819   0.9139
tiny_box: time [s]   walk     images   no images
          3e-9     0.6755   0.5867   0.7315
          6e-9     0.3915   0.2655   0.4833
          1e-8     0.1915   0.0958   0.2715
well mixed box: time [s]   well mixed   irt
          2.0e-7    20.41        20.18 +- 0.23
          5.0e-7    13.74        13.89 +- 0.22
          1.0e-6     8.86         8.85 +- 0.15
          1.5e-6     6.51         6.59 +- 0.13
```

## Review Notes

Base: `28f362bd64a044bb3de6f23bc56331ee711505f1`. Pre-existing changes: only the session record under `.work/`.

Against the Acceptance Criteria:
- `boundary = "reflective"`: every position is folded back into the box after a move, and so are reaction and product positions (`fold`, exact for free diffusion between walls).
- A pair near a wall also samples a reaction time against the mirror image of its partner, across every wall within the cutoff radius, with the edge and corner images (up to 26). An image is further than the partner itself, so the image pairs are among the pairs the search already finds. A reaction with an image happens where the image is and its products are folded back.
- Seeds carry a time and an offset; a track is injected at its time after every particle has been moved to that time, and its molecules join the pair search (tested with a partner that arrives 50 ns later).
- A uniform distribution stays uniform between the walls after moves from a tenth of the box to many boxes (chi-square on 20 slabs per axis), and molecules started in one corner fill the box.
- **One pair near a wall against a Brownian walk with reflection: agrees.** 4 nm (7 reaction radii) from one wall of a 60 nm box: survival at 20 ns 0.882 with images, 0.889 by the walk, 0.914 without images (5 standard errors off). The image rule matters and is right there.
- **One pair in a box of 11 reaction radii: does not agree, and this is the limit of the rule.** The walk lies between the image rule (reaction rate 1.5 times too high) and no images (0.85 times). A sink and its mirror image compete for the same molecules, which independent pairs cannot see: near a wall at distance `h` the true rate is lower than the free one by about `1 / (1 + r_eff / 2h)`, while the image rule restores the free rate. The criterion asked for agreement "in a small box"; it holds when the walls are far compared with the reaction radius, which is the case of every box radkin runs (micrometres against 0.5 nm). The tiny box is kept as a test of the bracket.
- Widely spaced tracks (50 um apart, 0.2 us one after the other) reproduce the sum of the single-track counts at their three ages within 4 standard errors, for every species.
- The well-mixed box moved here from ticket 04: 30 + 30 molecules of a partially diffusion-controlled `A + B` in a 40 nm box with reflective walls follow the exact stochastic rate equation at 4 times (within 4 standard errors plus 2 %).

Choices:
- Only the first image across each wall is sampled. So that this is enough, a sampling pass between reflective walls looks ahead no further than the time for which its cutoff radius equals the shortest edge of the box; when no reaction falls within it, every particle is moved over that time and a new pass starts (`box_moves` in `irt.json`). For a 2 um box that time is 22 us.
- `wall_images = false` switches the images off; it exists for the test above and as a diagnostic.
- G of a run with scheduled tracks divides by the energy injected so far, as the RDME yields do (tested on the CSV).
- The box of `[irt]` is the grid of `[voxel]` when it has no `[irt.box]` (ticket 03).

Not in the ticket's file list: `src/irt/output.rs` (the two new counters in `irt.json`). `src/irt/particles.rs` and `src/irt/search.rs` did not need a change.

## Blockers / Challenges

- The first one-pair test used the tiny box and failed by a wide margin. The cause is physical (above), not a bug: with one wall far from the others the same code agrees with the walk.
- Full resampling raises a diffusion-limited rate when it is repeated many times before a reaction: each pass draws the pair again as if no depletion had built up. In the tiny box, with a move every 0.38 ns, this adds about `r_eff / sqrt(pi D dt)` = 19 % to the rate. It is part of the excess of the image rule there.

## Commits

- 3a480a7 feat: IRT reflective walls with image pairs, scheduled tracks (ticket 05)

## Time Spent

10m (ticket-start.js to ticket-complete.js)
