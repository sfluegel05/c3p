"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone
A prenylquinone is defined as a quinone substituted by a polyprenyl‐derived side chain.
This version improves upon the previous heuristics by:
  - Using more relaxed SMARTS for quinone cores (allowing some substitutions)
  - Requiring that the prenyl fragment be directly attached to the quinone core.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    This routine requires that the molecule contains:
      1. A quinone core (detected by one of two SMARTS patterns that allow for substituents)
      2. A prenyl-derived side chain that is directly attached to the quinone core. 
         We look for a typical isoprene fragment "[CH2]C=C([CH3])", and verify that at least one
         matching fragment is bonded to an atom in the quinone core.
         
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a prenylquinone, False otherwise.
        str: Reason for classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for quinone cores.
    # The first is a relaxed pattern for a substituted 1,4-benzoquinone.
    # The second is a simplified naphthoquinone pattern.
    quinone_smarts = [
        "[#6]1:[#6]:[#6](=O):[#6]:[#6](=O):1",  # substituted 1,4-benzoquinone pattern
        "c1ccc2C(=O)c(c1)C(=O)cc2"              # simplified naphthoquinone pattern
    ]
    
    quinone_core_found = False
    quinone_atom_indices = set()
    for qs in quinone_smarts:
        pattern = Chem.MolFromSmarts(qs)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            # Take the first match and record its atom indices
            quinone_core_found = True
            # (We could union all, but one clear attached match is enough.)
            quinone_atom_indices.update(matches[0])
            break
    if not quinone_core_found:
        return False, "No quinone core detected"
    
    # Define a SMARTS pattern for a typical prenyl fragment:
    prenyl_pattern = Chem.MolFromSmarts("[CH2]C=C([CH3])")
    if prenyl_pattern is None:
        return False, "Error in prenyl SMARTS pattern"
    
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if not prenyl_matches:
        return False, "No prenyl fragment detected"
    
    # Require that at least one prenyl fragment is directly attached to the quinone core.
    attached_prenyl_found = False
    for match in prenyl_matches:
        # For each atom in the prenyl fragment match, check its neighbors
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
    
    # If both a quinone core and an attached prenyl fragment are present, then classify as prenylquinone.
    return True, "Molecule contains a quinone core with an attached prenyl-derived side chain"

# For simple testing, you might run:
if __name__ == "__main__":
    # Example: ubiquinone-2 from the provided examples.
    test_smiles = "COC1=C(OC)C(=O)C(C\\C=C(/C)CCC=C(C)C)=C(C)C1=O"
    result, reason = is_prenylquinone(test_smiles)
    print(f"Result: {result}, Reason: {reason}")