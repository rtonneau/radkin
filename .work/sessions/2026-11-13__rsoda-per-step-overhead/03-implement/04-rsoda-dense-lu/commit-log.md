# Ticket 04: rsoda-dense-lu

**Status:** ✅ Done

## Local Test Result

With the slice patch applied in rsoda: `cargo test --release` → 129 passed (golden included). radkin through a local `[patch]` path (not committed): A/B against the ticket-3 bench, 6x10. Dense -1.2% to -4.8% (about -3%). Sparse (does not use the code) +1% to +8%, i.e. the noise floor.
Microbench (22x22, 30% density): dgefa 2.9 → 1.9 us, dgesl 360 → 300 ns.
After the revert: the rsoda tree is clean on `perf/per-step-overhead` (08a2abf). radkin's checks are run in ticket 5 / at finish with the docs-only change.

## Review Notes

Base: `3558c05ff2c7269f76e5a9f0c6e40cc2751675ed`, rsoda `08a2abf`.
- Tried: column-slice AXPY loops in dense dgefa (multipliers, elimination) and dgesl (forward L, back U). Same per-element operations and order, so bit-identical (Rust does not contract to FMA).
- Outcome: below the 5% rule at solve level, so reverted, per the acceptance criteria. docs/perf records the numbers and the reason (dgesl is latency-bound on chained divisions; reciprocals would change the bits).
- Not tried: skipping zero multipliers in dgefa (factor is only 4-8%, and it changes the signs of zeros), banded LU (unused by radkin).

Dropped patch (rsoda `src/linalg.rs`):

```diff
diff --git a/src/linalg.rs b/src/linalg.rs
index c248a61..380a48b 100644
--- a/src/linalg.rs
+++ b/src/linalg.rs
@@ -74,19 +74,24 @@ pub fn dgefa(a: &mut [f64], lda: usize, n: usize, ipvt: &mut [usize]) -> usize {
 
         // Compute multipliers.
         let t = -1.0 / a[ix(k, k, lda)];
-        for i in (k + 1)..=n {
-            a[ix(i, k, lda)] *= t;
+        for v in &mut a[ix(k + 1, k, lda)..=ix(n, k, lda)] {
+            *v *= t;
         }
 
-        // Row elimination with column indexing.
+        // Row elimination with column indexing. The loops over rows
+        // `k + 1..=n` run on column slices (the same operations in the same
+        // order as the indexed Fortran loops, but free of bounds checks).
         for j in (k + 1)..=n {
             let t = a[ix(l, j, lda)];
             if l != k {
                 a[ix(l, j, lda)] = a[ix(k, j, lda)];
                 a[ix(k, j, lda)] = t;
             }
-            for i in (k + 1)..=n {
-                a[ix(i, j, lda)] += t * a[ix(i, k, lda)];
+            let (left, right) = a.split_at_mut(ix(1, j, lda));
+            let col_k = &left[ix(k + 1, k, lda)..=ix(n, k, lda)];
+            let col_j = &mut right[k..n];
+            for (x, m) in col_j.iter_mut().zip(col_k) {
+                *x += t * m;
             }
         }
     }
@@ -121,16 +126,18 @@ pub fn dgesl(a: &[f64], lda: usize, n: usize, ipvt: &[usize], b: &mut [f64], tra
                 b[l - 1] = b[k - 1];
                 b[k - 1] = t;
             }
-            for i in (k + 1)..=n {
-                b[i - 1] += t * a[ix(i, k, lda)];
+            let col = &a[ix(k + 1, k, lda)..=ix(n, k, lda)];
+            for (x, m) in b[k..n].iter_mut().zip(col) {
+                *x += t * m;
             }
         }
         // Solve U * x = y.
         for k in (1..=n).rev() {
             b[k - 1] /= a[ix(k, k, lda)];
             let t = -b[k - 1];
-            for i in 1..k {
-                b[i - 1] += t * a[ix(i, k, lda)];
+            let col = &a[ix(1, k, lda)..ix(k, k, lda)];
+            for (x, u) in b[..k - 1].iter_mut().zip(col) {
+                *x += t * u;
             }
         }
     } else {
```

## Blockers / Challenges

An early dgesl microbench let b underflow toward subnormals (repeated in-place solves). I fixed it by resetting b on every call.

## Commits

- 9479528 docs(perf): dense LU slice loops measured and dropped, below the 5% rule (ticket 04)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
