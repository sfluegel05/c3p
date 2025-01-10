"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: CHEBI:48477 alkanethiol
An alkanethiol is a compound in which a sulfanyl group, -SH, is attached to an alkyl group.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group, -SH, is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for at least one sulfanyl group (-SH)
    sulfanyl_pattern = Chem.MolFromSmarts("[SH]")
    sulfanyl_matches = mol.GetSubstructMatches(sulfanyl_pattern)
    if len(sulfanyl_matches) == 0:
        return False, "No sulfanyl group (-SH) found"

    # Check if the sulfanyl group is attached to an alkyl group (carbon chain)
    # We look for a carbon atom connected to the sulfur atom
    for match in sulfanyl_matches:
        sulfur_idx = match[0]
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        neighbors = sulfur_atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                # Check if the carbon is part of a simple alkyl chain
                # We allow the carbon to be part of a chain or ring, but we need to ensure the molecule is not too complex
                # Exclude molecules with too many heteroatoms or functional groups
                heteroatom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in {1, 6})
                if heteroatom_count > 1:  # Allow only one heteroatom (the sulfur in -SH)
                    return False, "Molecule contains too many heteroatoms or functional groups"
                
                # Check for simple alkyl chain: no rings, no double/triple bonds, no other functional groups
                if rdMolDescriptors.CalcNumRings(mol) > 0:
                    return False, "Molecule contains rings, which are not typical for simple alkanethiols"
                
                # Check for double or triple bonds
                if any(bond.GetBondType() not in {Chem.BondType.SINGLE} for bond in mol.GetBonds()):
                    return False, "Molecule contains double or triple bonds, which are not typical for simple alkanethiols"
                
                return True, "Sulfanyl group (-SH) attached to a simple alkyl group"

    return False, "Sulfanyl group (-SH) not attached to a simple alkyl group"