# Implementation Plan

**Session:** rdme-one-way-reverse-pairs
**Date:** 2026-11-15T14:51:49.348Z
**Estimated effort:** half a day

## Strategy

One ticket: infer the one-way reverse pairs at matrix build, feed them into the existing reversible-pair machinery, log them, and amend the docs.

## Tickets Overview

- **Ticket 1:** exact one-way reverse pairs (dynamic species only, unambiguous, one-way only) are inferred and get the shared factor. Unit-tested and documented.

## Sequencing Rationale

It is a single reviewable change on top of the previous session's ticket 01.

## Risks & Mitigation

- **Risk:** a false pairing in another network → **Mitigation:** the ambiguity rule, a log line for each inferred pair, and tests on CONFIG_PAPER_01 showing no new pairs where none are intended.
- **Risk:** a change in the reaction order or the CSV columns → **Mitigation:** pairs are metadata only; existing tests and golden checks must pass unchanged.

## Assumptions

- The previous session's `ReversiblePair` / `mesh_correct` / multi-grid logic accepts any (forward, backward) index pair.
