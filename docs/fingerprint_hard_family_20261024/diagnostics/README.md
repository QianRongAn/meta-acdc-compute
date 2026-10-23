# Retrospective hard-family error analysis

All 331 graphs and model outputs were exposed previously. Findings motivate developmental ablations; no independent evaluation is claimed.

| Family | Fold | Eligible queries | Composition | Model | Delta |
|---|---:|---:|---:|---:|---:|
| F293012a14526fbc3 | 2 | 2 | 1.0000 | 0.0000 | -1.0000 |
| Fb9d23d7d026efae5 | 0 | 2 | 1.0000 | 0.1667 | -0.8333 |
| F333328a21e5ae331 | 1 | 8 | 0.7375 | 0.5125 | -0.2250 |
| Fef92facf2b304718 | 3 | 7 | 0.9286 | 0.7500 | -0.1786 |
| F0465ecc03c008240 | 1 | 8 | 0.5923 | 0.4570 | -0.1353 |
| F19ed30cfb37ad7f2 | 2 | 17 | 0.4416 | 0.4432 | 0.0015 |
| Fde4a656c257aa9b3 | 0 | 6 | 0.3889 | 0.4444 | 0.0556 |
| F10427fd8efdc9dae | 4 | 2 | 0.0000 | 0.1667 | 0.1667 |
| F82d0ec93a01b1081 | 1 | 5 | 0.2667 | 0.5222 | 0.2556 |
| F5a0a5df565e94dcb | 0 | 2 | 0.2500 | 0.5833 | 0.3333 |
| F17560cbffb9e57b9 | 3 | 3 | 0.4167 | 0.7778 | 0.3611 |

Active positive anchors / all: 199827/476640.
Active anchors also seeing same-family negatives / active: 35710/199827.
Maximum identity-column mean after z-score: 1.07856025e-07.
Exact identical graph tuples across different receptor labels: 0 groups.
Zeroing one-hot means removes the direct composition mean, but does not prove all identity information is erased.
Edit distances refer to complete recorded chains; interface observability requires node-to-residue mapping and is not inferred here.
