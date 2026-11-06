# skills-science-upload

Claude skills for preparing and uploading data, and other materials, to public repositories and sites.

## Skills

| Skill | What it does | Version |
|-------|--------------|---------|
| [`upload-ncbi-sra`](upload-ncbi-sra/SKILL.md) | Fill NCBI BioSample and SRA submission templates from curated metadata, validate cross-links, and check and demultiplex FASTQ files | 1.00 |

New upload skills for other repositories or sites are added as sibling folders, each with its own `SKILL.md`.

## Layout

Each skill is a folder named after the skill:

```
<skill-name>/
├── SKILL.md      instructions and frontmatter
├── references/   background notes loaded when needed
└── scripts/      helper scripts the skill runs
```

## Install

With [`npx`](https://docs.npmjs.com/cli/commands/npx) (Node.js):

```bash
npx skills add kenscripts/skills-science-upload                     # choose skills and agent interactively
npx skills add kenscripts/skills-science-upload --list              # show the skills without installing
npx skills add kenscripts/skills-science-upload --skill upload-ncbi-sra -a claude-code -y
```

This installs into the current project (`./.claude/skills/`). Tested with Claude Code; the tool also supports other agents.

Or manually, for all your projects:

```bash
git clone https://github.com/kenscripts/skills-science-upload.git
mkdir -p ~/.claude/skills
cp -r skills-science-upload/upload-ncbi-sra ~/.claude/skills/    # or ln -s to follow git pull
```

For one project only, use `<project>/.claude/skills/` instead.

## Requirements

Installing a skill does not install its dependencies.

- `upload-ncbi-sra`: Python 3.10+ with `openpyxl` and `pandas`; `cutadapt` for demultiplexing. The same list is in the skill's `compatibility` field.

```bash
pip install openpyxl pandas cutadapt
```

Skills run with full agent permissions, so read a skill before using it.

## Author

Kenneth Acosta
