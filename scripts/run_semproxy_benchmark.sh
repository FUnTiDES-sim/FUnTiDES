#!/usr/bin/env bash
set -euo pipefail

# Defaults (override with flags below)
EXS="100"
EYS="100"
EZS="100"
IMPLEMS="geos shiva"
MESHES="cartesian ucartesian"
ORDERS="1 2 3"
SEMPROXY_BIN="bin/semproxy"
CSV_OUT="results.csv"
DRY_RUN=0

usage() {
  cat <<EOF
Usage: $0 [options]

Options (space-separated lists allowed; quoted):
  --ex        "100"            # values for --ex
  --ey        "100"            # values for --ey
  --ez        "100"            # values for --ez
  --implem    "geos"        # {classic, optim, geos, shiva}
  --mesh      "cartesian"      # {cartesian, ucartesian}
  -o|--order  "1 2 3"          # polynomial orders
  --bin       PATH             # path to semproxy (default: bin/semproxy)
  --out       FILE             # CSV output file (default: results.csv)
  --dry-run                    # print commands but don't execute
  -h|--help                    # show this help

Example:
  $0 --ex "100" --ey "100" --ez "100" --implem "shiva" --mesh "cartesian" -o "2 3"
EOF
}

# Parse long options (GNU getopt)
if ! PARSED=$(getopt -o ho: --long help,ex:,ey:,ez:,implem:,mesh:,order:,bin:,out:,dry-run -- "$@"); then
  usage; exit 2
fi
eval set -- "$PARSED"

while true; do
  case "$1" in
    --ex)       EXS="$2"; shift 2 ;;
    --ey)       EYS="$2"; shift 2 ;;
    --ez)       EZS="$2"; shift 2 ;;
    --implem)   IMPLEMS="$2"; shift 2 ;;
    --mesh)     MESHES="$2"; shift 2 ;;
    -o|--order) ORDERS="$2"; shift 2 ;;
    --bin)      SEMPROXY_BIN="$2"; shift 2 ;;
    --out)      CSV_OUT="$2"; shift 2 ;;
    --dry-run)  DRY_RUN=1; shift ;;
    -h|--help)  usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

# Check binary exists
if [[ ! -x "$SEMPROXY_BIN" ]]; then
  echo "Error: semproxy binary not found or not executable at: $SEMPROXY_BIN" >&2
  exit 1
fi

# Prepare CSV header
echo "ex,ey,ez,implem,mesh,order,time_s" > "$CSV_OUT"

# Iterate over the Cartesian product of options
# shellcheck disable=SC2086
for order in $ORDERS; do
  for mesh in $MESHES; do
    for ex in $EXS; do
      for ey in $EYS; do
        for ez in $EZS; do
          for implem in $IMPLEMS; do
            cmd=( "$SEMPROXY_BIN"
                  --ex "$ex" --ey "$ey" --ez "$ez"
                  --implem "$implem"
                  --mesh "$mesh"
                  -o "$order"
                  --timemax "1.5"
                  --dt "0.001" )
            echo ">> ${cmd[*]}"

            if [[ "$DRY_RUN" -eq 1 ]]; then
              # Record as NA when dry-running
              echo "$ex,$ey,$ez,$implem,$mesh,$order,NA" >> "$CSV_OUT"
              continue
            fi

            # Run and capture the “Elapsed Kernel Time” (number only)
            # Expected line: "---- Elapsed Kernel Time : 2.87449 seconds."
            output="$("${cmd[@]}" | sed -nE 's/.*Elapsed Kernel Time *: *([0-9.]+) seconds\..*/\1/p' | tail -n1)"

            if [[ -z "$output" ]]; then
              echo "Warning: could not parse time for config ex=$ex ey=$ey ez=$ez implem=$implem mesh=$mesh o=$order" >&2
              time_val="NA"
            else
              time_val="$output"
            fi

            echo "$ex,$ey,$ez,$implem,$mesh,$order,$time_val" >> "$CSV_OUT"
          done
        done
      done
    done
  done
done


# Pretty-print wide table (one column per order) to stdout
# Pretty-print wide table (mesh before implem, rows grouped by mesh then implem)
echo
echo "Results (wide format):"

awk -F, -v ord_list="$ORDERS" '
NR==1 { next }  # skip header
{
  # key = ex,ey,ez,mesh,implem   (note mesh before implem)
  key = $1 FS $2 FS $3 FS $5 FS $4
  data[key "," $6] = $7
  keys[key] = 1
}
END {
  # Parse order list so we control column order
  n = split(ord_list, ords, /[[:space:]]+/)

  # Header
  printf "ex\tey\tez\tmesh\timplem"
  for (i = 1; i <= n; i++) printf "\torder_%s", ords[i]
  print ""

  # Sort rows: ex,ey,ez,mesh,implem
  c = asorti(keys, skeys)
  for (i = 1; i <= c; i++) {
    k = skeys[i]
    split(k, a, FS)
    # a[1]=ex, a[2]=ey, a[3]=ez, a[4]=mesh, a[5]=implem
    printf "%s\t%s\t%s\t%s\t%s", a[1], a[2], a[3], a[4], a[5]
    for (j = 1; j <= n; j++) {
      val = data[k "," ords[j]]
      printf "\t%s", val
    }
    print ""
  }
}
' "$CSV_OUT" | column -t

echo
echo "CSV (long format) saved to: $CSV_OUT"
