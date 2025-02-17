"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone

A prenylquinone is defined as a quinone substituted by a polyprenyl‐derived side chain.
This implementation first tries to detect a quinone core by scanning for aromatic carbons
bearing a carbonyl (C=O) group. It then “grows” the core by including aromatic neighbors.
Finally, it requires that a prenyl fragment (here modeled by the isoprene fragment [CH2]C=C([CH3]))
is directly attached to this quinone core. 
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    
    The routine requires that the molecule contains:
      1. A quinone core: an aromatic ring (or fused ring system) having at least two carbons
         bearing a carbonyl group (C=O) as part of the aromatic system.
      2. A prenyl-derived side chain directly attached to the quinone core.
      
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a prenylquinone, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Identify candidate quinone core by finding aromatic carbons
    # that are directly double-bonded to an oxygen (i.e. a carbonyl group).
    aromatic_carbonyls = set()
    for atom in mol.GetAtoms():
        # Check for aromatic carbons.
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            # Look for a double bond to oxygen.
            for bond in atom.GetBonds():
                # Check that bond is a double bond.
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8:
                        aromatic_carbonyls.add(atom.GetIdx())
                        break

    quinone_core_indices = set()
    if len(aromatic_carbonyls) >= 2:
        # “Grow” the core by including aromatic neighbors of these carbonyl atoms.
        quinone_core_indices = set(aromatic_carbonyls)
        for idx in list(aromatic_carbonyls):
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic():
                    quinone_core_indices.add(neighbor.GetIdx())
    else:
        # As a fallback, try a set of predefined quinone SMARTS patterns.
        quinone_smarts = [
            "c1cc(=O)cc(=O)c1",                # Simple 1,4-benzoquinone
            "[#6]1:[#6]:[#6](=O):[#6]:[#6](=O):1",  # Alternate benzoquinone representation
            "c1ccc2C(=O)c(c1)C(=O)cc2",         # Naphthoquinone pattern variant 1
            "c1cc2cc(=O)cc(=O)c2c1"             # Naphthoquinone pattern variant 2
        ]
        for qs in quinone_smarts:
            pattern = Chem.MolFromSmarts(qs)
            if pattern is None:
                continue
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                # Use the first match to define the quinone core.
                quinone_core_indices = set(matches[0])
                break

    if not quinone_core_indices:
        return False, "No quinone core detected"
    
    # Step 2: Identify the prenyl fragment. We use a SMARTS to capture the common 
    # isoprene motif ([CH2]C=C([CH3])). This may need further extension to capture more prenyl variations.
    prenyl_pattern = Chem.MolFromSmarts("[CH2]C=C([CH3])")
    if prenyl_pattern is None:
        return False, "Error in prenyl SMARTS pattern"
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if not prenyl_matches:
        return False, "No prenyl fragment detected"
    
    # Step 3: Ensure that at least one prenyl fragment is directly attached (bonded) to the quinone core.
    attached_prenyl_found = False
    for match in prenyl_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in quinone_core_indices:
                    attached_prenyl_found = True
                    break
            if attached_prenyl_found:
                break
        if attached_prenyl_found:
            break

    if not attached_prenyl_found:
        return False, "Prenyl fragment not directly attached to the quinone core"

    return True, "Molecule contains a quinone core with an attached prenyl-derived side chain"

# For simple testing, you might run:
if __name__ == "__main__":
    # Test with one of the examples: ubiquinone-2.
    test_smiles = "COC1=C(OC)C(=O)C(C\\C=C(/C)CCC=C(C)C)=C(C)C1=O"
    result, reason = is_prenylquinone(test_smiles)
    print(f"Result: {result}, Reason: {reason}")