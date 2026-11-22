
#!/usr/bin/env bash

set -euo pipefail

N="${1:-5000}"

shuf -n "$N" --random-source=<(yes 42) reference/annotation.bed12 > reference/genebody_subset.bed12

wc -l reference/genebody_subset.bed12

