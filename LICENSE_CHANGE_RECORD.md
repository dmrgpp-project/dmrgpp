# GPL-3.0-or-later License Change Record

## Approval

The project owner confirmed on 2026-09-14 that UT-Battelle, L.L.C. has fully
approved the relicensing described in
[`LicenseChangePlan.md`](LicenseChangePlan.md) and remains the copyright holder
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

## Implementation gate

The approval, component scope, and source-header policy decisions are resolved.
Bulk relicensing remains gated on completion of the tracked-source provenance
and compatibility audit in Task 2 of `LicenseChangePlan.md`.
