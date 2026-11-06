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

Copy or symlink a skill folder into `~/.claude/skills/`, for example `ln -s "$PWD/upload-ncbi-sra" ~/.claude/skills/`.

## Author

Kenneth Acosta
