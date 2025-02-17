"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone
A prenylquinone is defined as a quinone substituted by a polyprenyl‐derived side chain.
This version improves upon the previous heuristics by broadening the quinone core SMARTS 
to capture common substituted benzoquinones and naphthoquinones.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    This routine requires that the molecule contains:
      1. A quinone core (detected by one or more SMARTS patterns that allow for substitution)
      2. A prenyl-derived side chain directly attached to an atom in the quinone core.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a prenylquinone, False otherwise.
        str: Reason explaining the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for quinone cores.
    # We include a more flexible set of patterns to capture benzoquinones and naphthoquinones.
    quinone_smarts = [
        "c1cc(=O)cc(=O)c1",               # Simple 1,4-benzoquinone pattern (aromatic ring with two C=O groups)
        "[#6]1:[#6]:[#6](=O):[#6]:[#6](=O):1", # Alternate benzoquinone representation
        "c1ccc2C(=O)c(c1)C(=O)cc2",        # Naphthoquinone pattern variant 1
        "c1cc2cc(=O)cc(=O)c2c1"            # Naphthoquinone pattern variant 2
    ]
    
    quinone_atom_indices = set()
    quinone_core_found = False
    # Try each quinone SMARTS pattern; if any pattern matches, we record the atom indices.
    for qs in quinone_smarts:
        pattern = Chem.MolFromSmarts(qs)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            quinone_core_found = True
            # Union all atoms from all matches of this pattern to capture the full core.
            for match in matches:
                for idx in match:
                    quinone_atom_indices.add(idx)
            # We break once we have at least one quinone core match.
            break
    if not quinone_core_found:
        return False, "No quinone core detected"
    
    # Define a SMARTS pattern for a typical prenyl fragment (isoprene-derived).
    # This pattern can capture the [CH2]C=C([CH3]) fragment.
    prenyl_pattern = Chem.MolFromSmarts("[CH2]C=C([CH3])")
    if prenyl_pattern is None:
        return False, "Error in prenyl SMARTS pattern"
    
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if not prenyl_matches:
        return False, "No prenyl fragment detected"
    
    # Look for at least one prenyl fragment that is directly attached to the quinone core.
    attached_prenyl_found = False
    for match in prenyl_matches:
        # For each atom in the prenyl fragment match, check neighboring atoms
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in quinone_atom_indices:
                    attached_prenyl_found = True
                    break
            if attached_prenyl_found:
                break
        if attached_prenyl_found:
            break
    
    if not attached_prenyl_found:
        return False, "Prenyl fragment not directly attached to the quinone core"
    
    # If both a quinone core and an attached prenyl fragment are present, classify as prenylquinone.
    return True, "Molecule contains a quinone core with an attached prenyl-derived side chain"

# For simple testing, you might run:
if __name__ == "__main__":
    # Example test: ubiquinone-2 from the provided examples.
    test_smiles = "COC1=C(OC)C(=O)C(C\\C=C(/C)CCC=C(C)C)=C(C)C1=O"
    result, reason = is_prenylquinone(test_smiles)
    print(f"Result: {result}, Reason: {reason}")