"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as any six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize SMARTS patterns for ketone and any double bond
    ketone_pattern = Chem.MolFromSmarts("C=O")
    double_bond_pattern = Chem.MolFromSmarts("C=C")

    # Get ring information
    ring_info = mol.GetRingInfo()

    # Check each ring to see if it matches the cyclohexenone criteria
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            # Extract the sub-molecule corresponding to the ring
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            ring_submol = Chem.PathToSubmol(mol, ring)
            
            # Check for non-aromatic ring
            if not rdMolDescriptors.CalcNumAromaticRings(ring_submol):
                # Check for ketone group
                if ring_submol.HasSubstructMatch(ketone_pattern):
                    # Check for exactly one additional C=C double bond
                    double_bond_count = 0
                    for atom in ring_atoms:
                        if atom.HasSubstructMatch(double_bond_pattern):
                            double_bond_count += 1

                    # Return True if have one ketone and one double bond
                    if double_bond_count == 1:
                        return True, "Contains a six-membered alicyclic ring with a ketone group and one other double bond"
    
    return False, "Does not meet cyclohexenone criteria"