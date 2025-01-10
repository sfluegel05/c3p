"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    Polar amino acids have side chains that can form hydrogen bonds, such as hydroxyl, amides, carboxyl, or basic nitrogen groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
        
        # Broadened amino acid backbone pattern
        backbone_patterns = [
            Chem.MolFromSmarts("N[C@@H]([C;R0&!R])[C](=O)O"), # Typical amino acids
            Chem.MolFromSmarts("N[C@H]([C;R0&!R])[C](=O)O")  # Include D-amino stereochemistry as well
        ]
        
        # Check for amino acid backbone
        if not any(mol.HasSubstructMatch(pattern) for pattern in backbone_patterns):
            return False, "No amino acid backbone found"

        # Expanded polar side chain groups
        polar_patterns = [
            Chem.MolFromSmarts("[NX3H2]"),  # Primary Amine
            Chem.MolFromSmarts("[OX2H]"),   # Hydroxyl
            Chem.MolFromSmarts("C(=O)[NX3]"), # Amide
            Chem.MolFromSmarts("C(=O)[OH]"), # Carboxylate
            Chem.MolFromSmarts("[nX2]"),    # Pyrrole-type Nitrogen
            Chem.MolFromSmarts("C=[NX3+]"), # Iminium Nitrogen
            Chem.MolFromSmarts("[SX2H]")    # Thiol
        ]

        # Check for any polar pattern match
        for pattern in polar_patterns:
            if mol.HasSubstructMatch(pattern):
                return True, "Contains a polar side chain capable of forming hydrogen bonds"

        return False, "No polar side chain found capable of forming hydrogen bonds"
    
    except Exception as e:
        return None, f"Error in processing SMILES: {str(e)}"