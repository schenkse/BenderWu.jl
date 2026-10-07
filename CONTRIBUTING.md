# Contributing

## Feature branches and pull requests

Start each feature or fix from the latest `main` on its own branch.
Open a draft pull request targeting `main` and mark it ready when the change is complete.
The `dev` branch is retired from the development workflow.

```bash
git fetch origin
git switch -c feat/my-feature origin/main
```

Use conventional commit titles, such as `feat: add a new potential helper` or `fix: validate perturbation orders`.
Keep each pull request focused on one change.
Describe the problem, the solution, and how it was verified.
Disclose substantial agent-generated changes in the pull request description.

Run the tests locally before requesting review:

```bash
julia --project=. --threads=2 -e 'using Pkg; Pkg.test()'
```

CI runs on branch pushes and pull requests into `main`.
It runs the package tests, including Aqua checks, on the latest stable Julia version and LTS across Linux, macOS, and Windows.
Julia nightly runs on Linux and is allowed to fail.
The `gates` check passes when all supported Julia jobs pass.
GitHub requires a pull request, an up-to-date branch, and passing `gates` before merging into `main`, including for admins.
Reviewer approval is optional.
Benchmarks remain a manual task.

Before opening a pull request, rebase the branch onto the latest `origin/main`.
After merging, start the next feature from `main` on a new branch.

## Releases

Prepare releases through a pull request into `main`.
Set `version` in `Project.toml` to the release version without the `-DEV` suffix and add `.github/release-notes/vX.Y.Z.md` with the user-facing changes.
Merge the release pull request, then tag that commit on `main`:

```bash
git fetch origin
git switch main
git pull --ff-only origin main
git tag vX.Y.Z
git push origin vX.Y.Z
```

Use a stable version tag such as `v1.0.6`.
The release workflow reruns CI, verifies that the tag matches `Project.toml`, and requires a nonempty release-notes file for that tag.
It creates a draft GitHub release for you to review and publish.
GitHub supplies the source archives; this Julia package has no separate build artifacts.

After a release, bump `Project.toml` to the next development version, such as `1.0.7-DEV`, through another pull request.
