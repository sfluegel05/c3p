"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: A cyclitol carboxylic acid (quinic acid and derivatives)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid derivative based on its SMILES string.
    For our purposes, quinic acid is defined as a cyclitol (here, a non‐aromatic cyclohexane)
    that carries a carboxylic acid group. In many natural derivatives the cyclohexane ring may 
    be decorated (e.g., by caffeoyl, feruloyl or coumaroyl esters) but the core remains a cyclohexanecarboxylate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as quinic acid derivative, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for the carboxylic acid group.
    # One for the protonated acid and one for the deprotonated form.
    acid_smarts = Chem.MolFromSmarts("C(=O)[O;H1]")   # C(=O)OH
    acid_smarts2 = Chem.MolFromSmarts("C(=O)[O-]")      # C(=O)O- (deprotonated)

    acid_matches = mol.GetSubstructMatches(acid_smarts)
    acid_matches += mol.GetSubstructMatches(acid_smarts2)
    
    if not acid_matches:
        return False, "No carboxylic acid group found"

    # Get ring information for use later
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    # Now, for each carboxylic acid match, verify that the acid group is connected to a cyclohexane ring.
    # In our SMARTS match, the carboxyl carbon is at index 0.
    for match in acid_matches:
        acid_carbon_idx = match[0]
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
        # Check the neighbors of the acid carbon (exclude O atoms from acid group)
        for neighbor in acid_carbon.GetNeighbors():
            # We expect the acid group to be attached to an sp3 carbon (atomic num 6) in a ring.
            if neighbor.GetAtomicNum() != 6:
                continue

            # Check if this neighboring carbon is in any ring of size 6
            for ring in atom_rings:
                if neighbor.GetIdx() in ring and len(ring) == 6:
                    # Verify that every atom in the ring is a non-aromatic carbon (sp3)
                    ring_ok = True
                    for idx in ring:
                        atom = mol.GetAtomWithIdx(idx)
                        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
                            ring_ok = False
                            break
                    if ring_ok:
                        # Found a cyclohexane ring attached to a carboxylic acid group.
                        return True, "Found cyclohexane carboxylic acid core (quinic acid derivative)"
    return False, "No cyclohexane ring attached to a carboxylic acid was detected"
    
# Uncomment below lines for quick testing:
# test_smiles = "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O"  # (+)-quinic acid
# result, reason = is_quinic_acid(test_smiles)
# print(result, reason)