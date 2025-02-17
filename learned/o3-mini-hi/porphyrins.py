"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins – natural pigments whose core structure is made of four pyrrole units connected by four methine bridges forming a macrocycle.
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is defined as a natural pigment whose core structure consists of 
    four pyrrole nuclei united through the alpha-positions by four methine groups
    to form a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a pyrrole-like ring.
    # Here we use a simple aromatic 5-membered ring containing an aromatic nitrogen.
    # Note: In many porphyrins the pyrrole nitrogen may be deprotonated or metal-coordinated
    # so we use "n" (aromatic nitrogen) without requiring an explicit hydrogen.
    pyrrole_pattern = Chem.MolFromSmarts("n1cccc1")
    if pyrrole_pattern is None:
        return False, "Error creating SMARTS pattern"

    # Find all substructure matches of the pyrrole unit.
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) < 4:
        return False, f"Found only {len(pyrrole_matches)} pyrrole-like rings; need at least 4 for a porphyrin core"
    
    # Heuristic step 1:
    # In a porphyrin the four pyrrole rings are fused by sharing four methine (carbon) atoms.
    # We iterate over all pairs of pyrrole ring matches and record the shared atom indices if they are carbons.
    bridging_carbons = set()
    n_matches = len(pyrrole_matches)
    for i in range(n_matches):
        for j in range(i+1, n_matches):
            # Convert tuples to sets of atom indices
            common_atoms = set(pyrrole_matches[i]).intersection(pyrrole_matches[j])
            for idx in common_atoms:
                # Check if the common atom is a carbon (atomic num 6).
                if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6:
                    bridging_carbons.add(idx)
    # In a porphyrin we expect to find four bridging methine carbons connecting the four pyrrole units.
    if len(bridging_carbons) < 4:
        return False, f"Only found {len(bridging_carbons)} bridging methine group(s); expected at least 4 for a porphyrin core"
    
    # Heuristic step 2 (optional):
    # One may check for overall conjugation or a large (macrocyclic) ring.
    # However, because the porphyrin macrocycle is a fused system, standard ring detection (SSSR)
    # may not return a single ring for the full periphery. We therefore rely on the connectivity of the pyrrole units.

    return True, "Molecule has at least four pyrrole-like rings fused via at least four bridging methine groups, consistent with a porphyrin core"