import requests, pkg_resources
from typing import Callable


def remove_altloc(pdbblock: str) -> str:
    """Removes all altlocs and actually segi duplicates.
    Test line:

    ATOM      1  N   MET A   1       0.000  0.000  0.000  1.00 60.69           N
    """
    lines: List[str] = []
    seen: List[Tuple[str, str, int]] = []
    for line in pdbblock.split('\n'):
        if 'ANISOU' in line:
            continue  # skip
        if line[:4] != 'ATOM' and line[:6] != 'HETATM':
            lines.append(line)
            continue
        atom_info = line[12:16].strip(), line[21].strip(), int(line[22:26].strip())
        if atom_info not in seen:
            lines.append(f'{line[:16]} {line[17:]}')
            seen.append(atom_info)
        else:  # skip
            pass
    return '\n'.join(lines)


def fetch(pdb_acc: str) -> str:
    response: requests.Response = requests.get(f'https://files.rcsb.org/download/{pdb_acc}.pdb')
    response.raise_for_status()
    return response.text

def get_data(resource: str) -> str:
    return pkg_resources.resource_string('DTC_compchem_practical', f'data/{resource}').decode()
