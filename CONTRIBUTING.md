# Contributing

Run `bash tests/run_all.sh` and the documented synthetic example before
submitting a change. Include regression tests for fixes and deterministic
output checks for changes to scientific computation.

Keep input identities, coordinate conventions, signed SVLEN, complete-background
weights and `0`/`1`/`NA` distinctions explicit. Input and scientific checks must
remain strict, with errors reported to the caller.

For a bug report include the software and dependency versions, exact command,
resource limits, error text and a minimal synthetic reproducer. Do not upload
credentials, identifiable human data or private paths. Report sensitive issues
privately to the maintainer through an appropriate private channel.

Contributions are accepted under the Apache License, Version 2.0.
