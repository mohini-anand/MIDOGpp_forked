---
description: Independently review a project premise — unit of analysis, a null, a mechanism claim, dataset semantics, a metric choice — against the data, derivation and literature rather than the repo's own decision record
argument-hint: [D5 | a premise or topic | blank for the standing register]
---

Launch the `premise-reviewer` subagent (Agent tool, `subagent_type: "premise-reviewer"`) on the
premise named below. Do not perform the review yourself in this session — its protocol lives in
`.claude/agents/premise-reviewer.md` and the whole point is that it runs with a fresh context and
treats this repo's documentation as the claim under test rather than as evidence.

Premise, decision, or topic: $ARGUMENTS

Rules for dispatching:

- If `$ARGUMENTS` is empty, do not dispatch the whole register — nine premises, each with a
  derivation and a literature search, is an unbounded run. Ask it for the **first three** in the
  register's own order (the matcher's mechanism, the exchangeable unit, the precision null), which
  are the cheapest to settle and the most widely depended on, and to report which it did not reach.
  Do not invent a premise of your own to hand it, and do not reorder the register.
- If the argument names a `DECISIONS.md` entry (`D5`, `D1`), pass that identifier through verbatim;
  the agent extracts the claim itself, including any amendments.
- Do not pass the project's reasoning along with the premise. Hand over the claim, not the case for
  it — a summary of why the repo believes it would defeat the agent's independence.
- Run it in the background. It derives, simulates and searches the literature, so it takes a while;
  the result arrives as a task notification.
- When it finishes, relay: the path of the register it wrote, the verdict per premise, any
  **unsound** verdict with its blast radius in one line, and which D-entries it believes now need
  an amendment. It will not have edited `DECISIONS.md` — that is the user's call.

This is the counterpart to `/audit-notebook`. That one asks whether a notebook's numbers are right
given the project's standards; this one asks whether the standards are.
