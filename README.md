# skills-science-upload

Claude skills for preparing and uploading data, and other materials, to public repositories and sites.

## Skills

| Skill | What it does | Version |
|-------|--------------|---------|
| [`upload-ncbi-sra`](upload-ncbi-sra/SKILL.md) | Fill NCBI BioSample and SRA submission templates from curated metadata, validate cross-links, and check and demultiplex FASTQ files | 1.00 |

## Install

Skills run with full agent permissions, so read a skill before using it.

### With npx (Node.js)

Installs into the current project (`./.claude/skills/`); tested with Claude Code. Omit `--skill` and `-a` to choose interactively.

```bash
npx skills add kenscripts/skills-science-upload --list
npx skills add kenscripts/skills-science-upload --skill upload-ncbi-sra -a claude-code -y
```

### Manually

Clone, then copy into `~/.claude/skills/` for all projects or `<project>/.claude/skills/` for one project. Use `ln -s` instead of `cp -r` to follow `git pull`.

```bash
git clone https://github.com/kenscripts/skills-science-upload.git
mkdir -p ~/.claude/skills
cp -r skills-science-upload/upload-ncbi-sra ~/.claude/skills/
```
