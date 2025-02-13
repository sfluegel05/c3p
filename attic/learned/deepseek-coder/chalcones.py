"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:23003 chalcone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone that is 1,3-diphenylpropenone (benzylideneacetophenone), ArCH=CH(=O)Ar, and its derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core chalcone pattern: Ar-CH=CH-C(=O)-Ar
    chalcone_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CH]=[CH]-[C](=O)-[c]2[c][c][c][c][c]2")
    if not mol.HasSubstructMatch(chalcone_pattern):
        return False, "Core chalcone structure (Ar-CH=CH-C(=O)-Ar) not found"

    # Check for the presence of the propenone (CH=CH(=O)) structure
    propenone_pattern = Chem.MolFromSmarts("[CH]=[CH]-[C](=O)")
    if not mol.HasSubstructMatch(propenone_pattern):
        return False, "Propenone (CH=CH(=O)) structure not found"

    # Check for two aromatic rings attached to the propenone structure
    aromatic_rings = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1")
    aromatic_matches = mol.GetSubstructMatches(aromatic_rings)
    if len(aromatic_matches) < 2:
        return False, "Less than two aromatic rings found"

    # Check that the aromatic rings are connected to the propenone structure
    connected_aromatics = 0
    for match in aromatic_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in neighbor.GetBonds()):
                    connected_aromatics += 1
                    break
    if connected_aromatics < 2:
        return False, "Aromatic rings not properly connected to the propenone structure"

    return True, "Contains core chalcone structure (Ar-CH=CH-C(=O)-Ar) with possible substitutions"