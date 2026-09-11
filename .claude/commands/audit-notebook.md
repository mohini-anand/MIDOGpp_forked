---
description: Independently audit an analysis notebook — recompute every number, check the math and config, rule on each conclusion
argument-hint: <path/to/notebook.ipynb> [extra instructions, e.g. "tier C" or "focus on the ROI-level null"]
---

Launch the `notebook-auditor` subagent (Agent tool, `subagent_type: "notebook-auditor"`) to audit
the notebook named below. Do not perform the audit yourself in this session — the agent's protocol
lives in `.claude/agents/notebook-auditor.md` and it must run with its own fresh context.

Target and any extra instructions: $ARGUMENTS

Rules for dispatching:

- If `$ARGUMENTS` is empty, do not guess a target. List the repo's notebooks newest-mtime first,
  and separately the untracked ones from `git status --porcelain`, then ask which to audit.
- If the argument is a bare notebook name rather than a path, resolve it by **exact match
  first**: `find . -name '<name>.ipynb' -not -path './.git/*'`. Only if that returns nothing, retry
  with `-name '<name>*.ipynb'`. This ordering matters because notebooks here are written by
  copying the previous one and named by suffixing it, so a short name is routinely a strict prefix
  of several longer ones and a glob-first lookup would silently resolve to a variant. Show the
  matches and ask which only when a lookup still returns more than one file — the same base name
  can also exist in more than one directory. Pass the full repo-relative path.
- Pass the resolved path to the agent verbatim, along with anything else in `$ARGUMENTS` as
  additional scope (for example a request to reach Tier C re-execution, or to concentrate on one
  section). Do not add scope the invoker did not ask for.
- Run it in the background. An audit reads several CSVs and re-derives every table, so it takes a
  while; the result arrives as a task notification.
- When it finishes, relay to the user: the path of the audit log it wrote, the per-conclusion
  verdicts as a compact list, each Tier 1 finding in one line, and the recompute tier it reached.
  Do not paste the whole log into the reply.
