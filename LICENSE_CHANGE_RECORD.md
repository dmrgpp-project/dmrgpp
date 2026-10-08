# GPL-3.0-or-later License Change Record

## Approval

The project owner confirmed on 2026-09-14 that UT-Battelle, L.L.C. has fully
approved the relicensing and remains the copyright holder
for the covered first-party code.

No public approval document or reference identifier was supplied for inclusion
in this repository. This file records the confirmation provided for the
repository change; it does not replace any underlying non-public approval
record.

## Approved license and scope

The approved license expression is:

`GPL-3.0-or-later`

It applies to all historical and current first-party contributions embodied in
the tracked repository, including first-party files outside the principal
component directories. The principal components are:

- root DMRG++ files;
- `dmrg/`;
- `LanczosPlusPlus/`;
- `PsimagLite/`; and
- `cincuenta/`.

Independently licensed third-party material is excluded from this blanket
statement and retains its applicable copyright notices and license terms.

## Approved notice and source-header policy

For first-party files in scope:

- use concise, comment-style-appropriate `SPDX-FileCopyrightText` declarations;
- preserve accurate copyright holders and existing year expressions;
- use `SPDX-License-Identifier: GPL-3.0-or-later`;
- remove the superseded full UT-Battelle license boilerplate from active
  first-party file headers; and
- preserve the complete former UT-Battelle license—including its copyright
  notice, conditions, DOE contract acknowledgment, disclaimer, and
  non-infringement statement—in a repository-level `NOTICE` file as a
  historical notice, without paraphrasing it.

Third-party and mixed-provenance files must be classified before editing. Their
notices must not be removed or replaced with a first-party-only SPDX
declaration.

## Third-party handling decision

The approved policy is to:

1. preserve the verified Loki and René Nyffenegger Base64 terms verbatim in
   `THIRD_PARTY_NOTICES`;
2. remove stale notice entries for files absent from the tracked source tree;
3. audit all tracked files for additional copied, adapted, translated,
   generated, or separately licensed material; and
4. block the license migration on material whose license is unresolved or
   incompatible with GPL-3.0-or-later until that material is removed, replaced,
   excluded from the covered distribution, or separately relicensed.

## Canonical source and verification summary

- Authoritative GPL text: <https://www.gnu.org/licenses/gpl-3.0.txt>
- Observed SHA-256 for `LICENSE`:
  `3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986`
- The root `LICENSE` was verified byte-for-byte against a fresh authoritative
  download.
- The tracked-file fingerprint and SPDX audits found no stale first-party
  declaration; all intended declarations use `GPL-3.0-or-later`.
- All paths and representative source markers recorded in
  `THIRD_PARTY_NOTICES` were verified.
- Fresh Clang builds and runtime-banner checks passed for `dmrg`,
  `toolboxdmrg`, and `lanczos`.
- The configured full build, focused migration-related tests, manual generation,
  path/dependency edge-case checks, and temporary staged installation passed.
- The broader CTest run passed 275 of 277 tests. Two reproducible Nightly
  numerical failures—`u0_gImpEqGcluster_dmrg_check` and
  `neqGBEKFig3L3NearAtomic`—are in cincuenta/GBEK paths not changed by the
  licensing migration and are recorded in the Task 8 result.
- Independent automated reviews of the completed Task 8 changes passed after
  all reported documentation-build defects were corrected. These reviews do
  not substitute for the required human authority and component-maintainer
  reviews.

Detailed audit evidence is in `LICENSE_AUDIT.md`; task-by-task implementation
and verification results are in `LicenseChangePlan.md`.

## Final human-review gate

Task 9 remains open. On 2026-09-15, the requester directed that the focused
commit series be retained and that this consolidated review record be prepared,
but that Task 9 remain pending for human review.

Before landing the series, obtain and record:

1. final review by the approving copyright-holder/legal authority identified
   above; and
2. final review by maintainers familiar with every included first-party
   component: root DMRG++, `dmrg/`, `LanczosPlusPlus/`, `PsimagLite/`, and
   `cincuenta/`.

Do not mark Task 9 complete or land the migration until those reviews are
confirmed. No final human-review approval is claimed by this record.
