# Guarded recursive delete for the test harnesses.  Source this; do not execute it.
#
#     source "$repo/tests/lib/safe_rm.sh"
#     safe_rm_rf "$repo/tests" "$out_py" "$out_rs"
#
# WHY THIS EXISTS
#
# `rm -rf "$x"` is only as safe as `$x`.  `set -u` catches an *unset* variable but NOT an
# empty one, and every deletion path in this harness is built by concatenating other
# variables with a user-supplied dataset name:
#
#     out_py="$build/out/$name/py"
#
# If `$build` were ever empty that becomes `rm -rf "/out//py"`, and if `$name` were `..` it
# escapes the build directory entirely.  Neither is caught by `set -euo pipefail`.  The
# harness deletes directories on essentially every run, so the blast radius of getting this
# wrong once is the developer's working tree.
#
# `safe_rm_rf` refuses rather than guesses.  The first argument is a mandatory ROOT: the
# directory every target must live under.  That is what turns "delete this path" into
# "delete this path *inside the sandbox I expect*".
#
# NOT WRAPPED, DELIBERATELY: the `tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT` sites.
# `mktemp -d` either prints a fresh absolute path or fails, and under `set -e` a failure
# aborts before the trap is installed, so `$tmp` cannot be empty or attacker-shaped there.

# safe_rm_rf ROOT TARGET...
#
# Removes each TARGET recursively.  Every TARGET is validated BEFORE any is deleted, so a
# bad third argument cannot leave the first two already gone.  A TARGET whose parent does
# not exist is skipped (there is nothing to delete), not treated as an error.
safe_rm_rf() {
    local root="$1"; shift || {
        echo "safe_rm_rf: no ROOT given" >&2; return 2
    }
    if [[ -z "$root" ]]; then
        echo "safe_rm_rf: ROOT is empty (unset or failed command substitution?)" >&2
        return 2
    fi

    local root_abs
    root_abs="$(cd "$root" 2>/dev/null && pwd -P)" || {
        echo "safe_rm_rf: ROOT does not exist: $root" >&2; return 2
    }

    local p parent base parent_abs abs
    local -a approved=()

    # Pass 1 -- validate everything.
    for p in "$@"; do
        if [[ -z "$p" ]]; then
            echo "safe_rm_rf: refusing to delete an empty path" >&2
            return 2
        fi
        if [[ "$p" == *..* ]]; then
            echo "safe_rm_rf: refusing a path containing '..': $p" >&2
            return 2
        fi

        parent="$(dirname -- "$p")"
        base="$(basename -- "$p")"
        # A missing parent means the target cannot exist. Nothing to do, and resolving it
        # would fail -- so skip rather than refuse.
        parent_abs="$(cd "$parent" 2>/dev/null && pwd -P)" || continue
        abs="$parent_abs/$base"

        # `cd`+`pwd -P` above already resolved any symlinked ANCESTOR, so this comparison
        # cannot be fooled by a link partway up the path.
        if [[ "$abs" != "$root_abs"/* ]]; then
            echo "safe_rm_rf: refusing $abs -- outside permitted root $root_abs" >&2
            return 2
        fi
        if [[ "$abs" == "$root_abs" ]]; then
            echo "safe_rm_rf: refusing to delete the root itself: $root_abs" >&2
            return 2
        fi
        if [[ -L "$abs" ]]; then
            echo "safe_rm_rf: refusing to delete a symlink: $abs" >&2
            return 2
        fi
        if [[ -e "$abs" && ! -d "$abs" ]]; then
            echo "safe_rm_rf: refusing to delete a non-directory: $abs" >&2
            return 2
        fi
        approved+=("$abs")
    done

    # Pass 2 -- delete.
    local a
    for a in "${approved[@]+"${approved[@]}"}"; do
        rm -rf -- "$a"
    done
}
