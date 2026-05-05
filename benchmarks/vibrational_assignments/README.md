# ORCAVEDA vibrational assignment benchmark seed

This directory contains a small curated seed benchmark for checking ORCAVEDA
normal-mode assignment labels against traceable molecule-specific sources.

ORCAVEDA Stage 3D is a geometric and weighted independent-coordinate assignment
audit. These rows are not a claim that ORCAVEDA implements strict VEDA PED or a
full Wilson GF normal-coordinate analysis.

## Confidence levels

- `gold`: explicit molecule-specific vibrational assignment in a trusted
  compilation or peer-reviewed assignment table.
- `silver`: molecule-specific frequency evidence plus a literature-supported
  assignment that is plausible but not explicitly assigned in the same source
  row.
- `weak`: generic functional-group frequency heuristic only.

Rows should stay `silver` or `weak` until the exact assignment text is backed by
a source that reports that molecule and mode.

## Source notes

NIST Chemistry WebBook / Shimanouchi rows for water, ammonia, and acetaldehyde
include molecule-specific approximate mode labels and selected frequencies.

The aniline rows come from Ilic et al., Langmuir 2000. That article reports
calculated B3LYP/6-31G* frequencies, liquid IR, liquid Raman, gas IR reference
bands, and assignments for characteristic aniline modes. Mixed-mode labels are
preserved as mixed.

The benzene rows come from NIST Chemistry WebBook / Shimanouchi and include
IR/Raman activity and symmetry information. Benzene modes are useful for
checking aromatic C-H bends, aromatic ring stretches, ring deformation, and
ring-breathing labels without substituent effects.

The benzoic acid rows come from Bakker et al., J. Chem. Phys. 2003. They use
jet-cooled gas-phase monomer IR ion-dip spectroscopy and B3LYP/D95(d,p)
assignments. Mixed carboxylic-acid/ring assignments are preserved as mixed
rather than forced to a single functional group.

The acetophenone rows come from Attia and Schauermann, J. Phys. Chem. C 2020.
They use multilayer acetophenone IRAS on Pt(111) plus scaled gas-phase
MP2/aug-cc-pVQZ assignments. The surface/multilayer context is recorded in
each row because some bands, especially C=O, are shifted relative to gas-phase
NIST values.

The pyridine rows come from Wong and Colson, J. Mol. Spectrosc. 1984. They are
high-resolution gas-phase FT-IR fundamental assignments with symmetry labels.
Many rows intentionally use `pyridine ring system` instead of atomistic labels
because the source reports term values and symmetry assignments rather than
full displacement/PED descriptions.

The phenol rows come from Billes and Mohammed-Ziegler, Applied Spectroscopy
Reviews 2007. They currently cover O-H and aromatic C-H stretching frequencies
from compiled gas-phase IR, diluted-solution IR, density-functional, and DFT
sources. Phenol still needs lower-frequency ring, C-O, and O-H bending rows.

Acetamide rows are intentionally `silver`. NIST WebBook provides gas and
condensed-phase IR spectra but not assignment labels. CCCBDB provides
experimental fundamental frequencies for acetamide; the current labels use
amide-region literature context and must not be promoted to `gold` until a
source table with explicit acetamide assignments is added.

The acetamide far-IR rows from Kydd and Dunham are molecule-specific
vapor-phase assignments and may be `gold` where the authors make direct
assignments. Rows with hedging language such as `probably` or `may also be
seen` stay `silver` or `weak`.
