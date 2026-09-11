# CLAUDE.md — pointer

**The project instructions live in [`markdown_files/CLAUDE.md`](markdown_files/CLAUDE.md).**
Read that file: it is the repository map, the physics switches, and the standing rules.

---

## Why this stub exists

`CLAUDE.md` was moved into `markdown_files/` on 2026-09-11 so that all tracked documentation sits in
one directory. But Claude Code auto-loads `CLAUDE.md` **from the repository root and parent
directories only** — a copy inside `markdown_files/` is loaded only when work happens inside that
directory. Moving it without leaving this pointer would silently stop the standing rules from
reaching any session started at the root, which is where essentially all work happens.

This stub keeps the mechanism working while the content lives where the rest of the documentation
does. It is deliberately two lines of substance: there is **one** source of truth, and it is
`markdown_files/CLAUDE.md`.

⚠ If you delete this file, the rules stop loading automatically. The alternative, if a stub is
unwanted, is to move the real file back to the root and leave a pointer in `markdown_files/`
instead.
