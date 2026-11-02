# Shared shell settings. Launch scripts are run from the repository root.
# Values already supplied in the environment override the versioned defaults.
while IFS=$'\t' read -r key value; do
  [[ "$key" == key || -z "$key" ]] && continue
  [[ "$key" =~ ^CH1_[A-Z_]+$ ]] || { echo "Invalid configuration key: $key" >&2; return 1; }
  export "$key=${!key:-$value}"
done < scripts/2-hmsc/run-config.tsv
for key in CH1_SAMPLES CH1_THIN CH1_CHAINS CH1_NFOLDS CH1_POSTERIOR_THIN; do
  [[ "${!key}" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid integer: $key" >&2; return 1; }
done
[[ "$CH1_TRANSIENT" =~ ^[0-9]+$ ]] || { echo 'Invalid transient' >&2; return 1; }
