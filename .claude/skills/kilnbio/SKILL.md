---
name: kilnbio
description: Use the KilnBio research pipeline to plan assay processing, run configured multiomics studies, audit assay eligibility, or compare calls across analysis methods.
---

# Analyze a study with KilnBio

Use the version in the researcher's checkout. In a shell-capable session, start from the repository root and read `VERSION`, `CLI.txt` and the guide for the requested route. With the KilnBio MCP connection, call `list_guides` and `read_guide` to get those instructions from the installed version. In a chat without tools or a checkout, give the researcher the corresponding local commands and identify the inputs they need. The [research repository](https://github.com/varunorama/kilnbio-research) provides the guides and [access terms](https://github.com/varunorama/kilnbio-research/blob/main/ACCESS.md).

- For raw sequencing reads, use the assay's sample sheet and reference requirements in `docs/USER_GUIDE.md`. For existing measurements, start at the supported measurement or count-table entry point. Identify the biological unit, assay scale and contrast before choosing the command.
- For an execution-enabled MCP project, call `list_pipelines` to see the currently installed studies. When the researcher requests execution, use `run_pipeline(pipeline_id, run_id)` with a new descriptive run ID. Reuse that ID if delivery is uncertain. Poll `pipeline_status(run_id)` and report completion from the receipt with the local artifact locations. Read project-specific data only in an authorized environment. For a local CLI study, use `kilnbio study STUDY.json RUN_DIR --dry-run` to preview it. See `docs/STUDY_RUNNER.md` for execution and resume behavior.
- To count who remains eligible when adding an assay, use `assay_reach` with one authorized pseudonymous row per donor, or `kilnbio assay-audit reach` on the local roster. Distinguish donor availability from same-specimen pairing. For the latter, use `assay-audit pair` with specimen identifiers.
- To compare methods, use `compare_methods` for small authorized summary tables or `kilnbio method-agreement` for local files. Match the contrast direction, feature family, units, assembly and assay channel. Report planned and tested denominators, sign conflicts and missing calls. Concordance between methods applied to the same cohort is descriptive; independent validation requires additional data.

Keep restricted inputs in the researcher's approved environment. Use CLI paths for a local run instead of pasting participant data into chat. Base reported findings on completed outputs and record the version and analysis settings.
