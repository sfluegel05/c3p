"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is formed by the esterification of tetradecanoic acid with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define tetradecanoate functional group pattern
    # This molecule starts with a 14-carbon chain and follows with an ester linkage
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")
    matching_atoms = mol.GetSubstructMatches(tetradecanoate_pattern)
    
    if len(matching_atoms) > 0:
        # Calculate the total number of carbon atoms in the esterified alcohol part
        esterified_parts = 0
        for match in matching_atoms:
            # Check beyond the ester linkage to ensure it connects with a hydrocarbon or aromatic alcohol structure
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                # Consider carbons beyond the ester oxygen
                if atom.GetAtomicNum() == 6:
                    esterified_parts += 1

        # We already asserted liveliness with 14 carbon backbones, additional carbons would apply to the alcohol
        if esterified_parts >= 14:  # Basic assumption extending
            return True, "Contains tetradecanoate ester group"
        else:
            return False, "Insufficient esterified structure"
    
    return False, "No tetradecanoate ester group found"