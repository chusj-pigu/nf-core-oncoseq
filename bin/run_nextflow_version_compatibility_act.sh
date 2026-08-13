#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "$script_dir/.." && pwd)"
workflow="$repo_root/.github/workflows/nextflow-version-compatibility.yml"

event_name="${ACT_EVENT:-workflow_dispatch}"
job_name="${ACT_JOB:-test}"
runner_image="${ACT_RUNNER_IMAGE:-catthehacker/ubuntu:full-latest}"
cache_dir="${ACT_CACHE_DIR:-${TMPDIR:-/tmp}/oncoseq-act}"

if ! command -v act >/dev/null 2>&1; then
    echo "Required command not found: act" >&2
    exit 1
fi

if [[ ! -f "$workflow" ]]; then
    echo "Compatibility workflow not found: $workflow" >&2
    exit 1
fi

commit="$(git -C "$repo_root" rev-parse --verify HEAD)"
branch="$(git -C "$repo_root" branch --show-current)"
mkdir -p "$cache_dir"

printf 'Running %s on %s (%s)\n' "$workflow" "${branch:-detached HEAD}" "$commit"
printf 'Event: %s | Job: %s | Runner image: %s\n' "$event_name" "$job_name" "$runner_image"

cd "$repo_root"
exec act "$event_name" \
    --workflows "$workflow" \
    --job "$job_name" \
    --bind \
    --platform "ubuntu-latest=$runner_image" \
    --action-cache-path "$cache_dir" \
    "$@"
