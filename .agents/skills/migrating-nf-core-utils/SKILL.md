---
name: migrating-nf-core-utils
description: Migrates nf-core pipelines from vendored utility subworkflows and local helpers to the nf-core-utils Nextflow plugin. Use when replacing utils_nextflow_pipeline, utils_nfcore_pipeline or utils_nfschema_plugin, updating nf-core-utils or nf-schema pins, or diagnosing plugin function errors after an update.
license: MIT
---

# Migrate to nf-core-utils

Replace vendored helpers with calls supported by the pinned plugin versions. Preserve pipeline-specific validation, reporting and early-exit behaviour. Keep nf-test infrastructure changes out of the migration. Tool citations are a separate follow-up, covered in [citations.md](references/citations.md).

For diagnosis or pin-only work, fix only the confirmed mismatch and run relevant checks; do not migrate other helpers unless requested. For a partial migration, replace only the remaining helpers. Do not use this workflow for nf-test helper-only work.

## Inspect the pipeline and target versions

Read the pipeline's contribution rules before editing. Inventory the existing imports, local helpers, vendored subworkflows and plugin pins:

```bash
rg -n 'utils_(nextflow_pipeline|nfcore_pipeline|nfschema_plugin)' modules.json main.nf workflows subworkflows
rg -n 'plugin/|^def ' main.nf workflows subworkflows/local
rg -n 'plugins\s*\{|nf-core-utils@|nf-schema@' nextflow.config conf .github
```

Search any additional configs passed with `-c` or supplied by Seqera Platform. Compare local helpers with the exports of the target plugin version; keep helpers that contain pipeline-specific logic.

Read signatures from the exact releases being pinned, not from another pipeline's migration. For an installed plugin, substitute the target version for `VERSION`:

```bash
nextflow plugin install nf-core-utils@VERSION
P="$HOME/.nextflow/plugins/nf-core-utils-VERSION"
cat "$P/classes/META-INF/MANIFEST.MF"
javap -cp "$P/classes" nfcore.plugin.NfUtilsExtension
```

Inspect `nextflow.validation.ValidationExtension` in the corresponding `nf-schema-VERSION/classes` directory too. If the plugin cache layout differs or `javap` is unavailable, read the tagged plugin source. Check the plugin's Nextflow requirement against `manifest.nextflowVersion`, CI's minimum `NXF_VER` and the README badge; update them together if needed.

## Replace helpers and remove unused copies

Use [function-map.md](references/function-map.md) for version-labelled examples. Add explicit pins for both plugins in `nextflow.config`. Replace each helper's imports and calls before removing it.

- Preserve custom help text, summary fields, validation options and completion hooks.
- Preserve `cast_cli_params: false` if the old validation call used it.
- Preserve the git suffix in version output with `getWorkflowVersion(workflow.manifest.version, workflow.commitId)` when required.
- Pass a string path to `methodsDescriptionText`, for example `mqc_methods_yaml.toString()`.
- Keep version tuples intact when collecting them. Use `.toList()` when an empty stream must still emit a list. For legacy version-file channels, read each `Path` with `.text` before passing YAML content if the target API requires it; this also supports remote work directories.

Remove vendored subworkflows only when no callers remain. Prefer `nf-core subworkflows remove <name>` to keep `modules.json` consistent. Remove unused local helper definitions and imports, but retain customised citation text until it has an equivalent replacement.

## Pin launch-time configs

External `plugins {}` blocks can change the effective plugin list. Repeat the pipeline's exact `nf-core-utils` and `nf-schema` pins in each launch-time block, alongside any additional plugins. Check `.github/workflows/awstest.yml`, `awsfulltest.yml`, Platform configs and supplied `-c` files, where present.

Reproduce the external launch configuration and confirm the loaded versions in `.nextflow.log`; matching text in the repository alone does not verify a Platform run.

## Diagnose function errors

`Missing process or function X(...) -- Did you mean 'X' instead?` can indicate an overload mismatch, not just a misspelt name. Obtain the failing run's `.nextflow.log` or Platform/CI logs and check:

- Which plugin version actually loaded, including `Downloading plugin` or `Starting plugin` messages.
- Whether that version exports the function with the supplied argument count and types.
- Whether a local definition or external config changes the intended call or pin.

Correct the confirmed mismatch; do not bump plugins or rename functions based only on the error text. If the logs are unavailable, report the missing evidence.

## Verify the migration

For a helper migration, run the checks below. For a pin-only fix, run the failing launch and checks covering the affected calls.

- Run `nextflow lint` on changed Nextflow files when supported by the target Nextflow version, and `nf-core pipelines lint`.
- Run the pipeline's nf-tests with its minimum supported Nextflow version and an available container profile.
- Run `nextflow run . --help`, `nextflow run . --help_full` and `nextflow run . --version`; check their output and successful early exit without launching analysis tasks.
- Exercise a small test-profile run with the external config used by CI or Platform. Confirm plugin versions and inspect the generated software-version YAML and MultiQC methods description.
- Review snapshot changes before regenerating them through nf-test. The plugin may quote versions previously represented as numbers; distinguish quoting changes from missing or altered versions. Do not hand-edit snapshots.

Before calling a failure pre-existing, run the same command on the pipeline's unchanged base revision. Report checks that could not run. Follow the pipeline's changelog and PR conventions; summarise plugin versions, external config pins and any helpers kept local. Opening or merging a PR requires the user's authority.

## Attribution

Adapted from [Maxime U Garcia's nf-core-utils-migration](https://github.com/maxulysse/nf-core-utils-migration/tree/844f10cfc45f63154f027f28e0f20bc5825bf656), under the [MIT licence](LICENSE). The original guidance draws on migrations in [rnavar#233](https://github.com/nf-core/rnavar/pull/233), [rnavar#252](https://github.com/nf-core/rnavar/pull/252), [seqinspector#241](https://github.com/nf-core/seqinspector/pull/241), [rnaseq#1977](https://github.com/nf-core/rnaseq/pull/1977) and [createpanelrefs#107](https://github.com/nf-core/createpanelrefs/pull/107).
