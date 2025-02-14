"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: CHEBI:33567 aromatic amino acid
An amino acid whose structure includes an aromatic ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for amino acid pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NX3H2][CX4H]([CX3](=O)[OX1-,OX2H1])[CX4H]")
    amino_acid_match = mol.GetSubstructMatches(amino_acid_pattern)
    if not amino_acid_match:
        return False, "No amino acid pattern found"
    
    # Check if any atom in the amino acid backbone is aromatic
    amino_acid_atoms = [mol.GetAtomWithIdx(idx) for idx in amino_acid_match[0]]
    aromatic_atoms = [atom for atom in amino_acid_atoms if atom.GetIsAromatic()]
    if not aromatic_atoms:
        return False, "No aromatic atoms in amino acid backbone"
    
    # Additional checks for increased confidence
    has_carboxyl = any(atom.GetDegree() == 2 and atom.GetFormalCharge() == -1 for atom in mol.GetAtoms())
    has_amino = any(atom.GetDegree() == 3 and atom.GetFormalCharge() == 0 and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 3 for atom in mol.GetAtoms())
    
    if has_carboxyl and has_amino:
        return True, "Contains an aromatic ring in the amino acid backbone, with carboxyl and amino groups"
    else:
        return False, "Missing carboxyl or amino group"