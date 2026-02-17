#!/bin/bash
# Run FEAT once, and if FLAME fails with excessive outliers, retry with robust_yn=0
# using a fresh output directory to avoid stale intermediate state.

set -u

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
  echo "Usage: $0 <feat_binary> <fsf_file> [njobs]" >&2
  exit 2
fi

feat_binary="$1"
fsf="$2"
njobs="${3:-}"

if [ ! -x "$feat_binary" ]; then
  echo "run_feat_with_outlier_retry: FEAT binary is not executable: $feat_binary" >&2
  exit 2
fi

if [ ! -r "$fsf" ]; then
  echo "run_feat_with_outlier_retry: FSF file not readable: $fsf" >&2
  exit 2
fi

run_feat() {
  local this_fsf="$1"
  local this_njobs="${2:-}"
  if [ -n "$this_njobs" ]; then
    "$feat_binary" "$this_fsf" -P "$this_njobs"
  else
    "$feat_binary" "$this_fsf"
  fi
}

extract_outputdir() {
  local fsf_file="$1"
  local out
  out=$(sed -n 's/^set fmri(outputdir)[[:space:]]*"\(.*\)"[[:space:]]*$/\1/p' "$fsf_file" | head -n 1)
  if [ -z "$out" ]; then
    out="${fsf_file%.fsf}"
  fi
  printf '%s\n' "$out"
}

extract_robust_yn() {
  local fsf_file="$1"
  awk '/^set fmri\(robust_yn\)[[:space:]]+/ { print $3; exit }' "$fsf_file"
}

has_excessive_outliers_error() {
  local feat_dir="$1"
  if [ ! -d "$feat_dir" ]; then
    return 1
  fi
  grep -Rqsi "excessive.*outliers detected" "$feat_dir"
}

canonical_outputdir=$(extract_outputdir "$fsf")
canonical_gfeat="${canonical_outputdir}.gfeat"

run_feat "$fsf" "$njobs"
exit_code=$?

if [ $exit_code -eq 0 ]; then
  exit 0
fi

robust_yn=$(extract_robust_yn "$fsf")
if [ "${robust_yn:-}" = "0" ]; then
  exit "$exit_code"
fi

if ! has_excessive_outliers_error "$canonical_gfeat"; then
  exit "$exit_code"
fi

retry_outputdir="${canonical_outputdir}.retry_no_outlier"
retry_gfeat="${retry_outputdir}.gfeat"
retry_fsf="${fsf%.fsf}__retry_no_outlier.fsf"
retry_started_at=$(date)

echo "WARNING: FEAT detected excessive FLAME outliers for ${canonical_gfeat}." >&2
echo "WARNING: Automatically retrying with robust outlier deweighting disabled (robust_yn=0)." >&2
echo "WARNING: Retry output directory: ${retry_gfeat}" >&2

if [ -e "$retry_fsf" ]; then
  rm -f "$retry_fsf"
fi

if [ -d "$retry_gfeat" ]; then
  stale_suffix=$(date +%Y%m%d-%H%M%S)
  mv "$retry_gfeat" "${retry_gfeat}.stale_${stale_suffix}"
fi

awk -v retry_out="$retry_outputdir" '
  BEGIN { set_robust = 0; set_output = 0 }
  {
    if ($0 ~ /^set fmri\(robust_yn\)[[:space:]]+/) {
      print "set fmri(robust_yn) 0"
      set_robust = 1
      next
    }
    if ($0 ~ /^set fmri\(outputdir\)[[:space:]]+/) {
      print "set fmri(outputdir) \"" retry_out "\""
      set_output = 1
      next
    }
    print
  }
  END {
    if (set_robust == 0) print "set fmri(robust_yn) 0"
    if (set_output == 0) print "set fmri(outputdir) \"" retry_out "\""
  }
' "$fsf" > "$retry_fsf"

run_feat "$retry_fsf" "$njobs"
retry_code=$?

if [ $retry_code -ne 0 ]; then
  echo "WARNING: Retry with robust_yn=0 also failed for ${fsf}." >&2
  exit "$retry_code"
fi

archived_initial=""
if [ -d "$canonical_gfeat" ]; then
  archived_initial="${canonical_gfeat}.failed_robust_$(date +%Y%m%d-%H%M%S)"
  mv "$canonical_gfeat" "$archived_initial"
fi

if [ ! -d "$retry_gfeat" ]; then
  echo "run_feat_with_outlier_retry: retry succeeded but output directory missing: ${retry_gfeat}" >&2
  exit 1
fi

mv "$retry_gfeat" "$canonical_gfeat"

if [ -d "$canonical_gfeat" ]; then
  cp "$retry_fsf" "${canonical_gfeat}/design_retry_no_outlier.fsf" 2>/dev/null || true
  {
    echo "Automatic FEAT retry triggered."
    echo "Reason: Excessive number of FLAME outliers detected in initial run."
    echo "Original FSF: ${fsf}"
    echo "Initial robust_yn: ${robust_yn:-NA}"
    echo "Retry FSF: ${retry_fsf}"
    echo "Retry started: ${retry_started_at}"
    if [ -n "$archived_initial" ]; then
      echo "Archived initial output: ${archived_initial}"
    fi
  } > "${canonical_gfeat}/.feat_auto_retry_warning"
fi

# Clean up retry artifacts now that the retry succeeded and provenance is
# preserved in .feat_auto_retry_warning inside the canonical .gfeat.
rm -f "$retry_fsf" 2>/dev/null || true
if [ -n "$archived_initial" ] && [ -d "$archived_initial" ]; then
  rm -rf "$archived_initial"
fi

echo "WARNING: FEAT retry with robust_yn=0 succeeded. Final output written to ${canonical_gfeat}." >&2
exit 0
