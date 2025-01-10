"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    Lysophosphatidic acids are characterized by one fatty acid chain, a glycerol backbone, 
    and a terminal phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with a terminal phosphate group
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[OX2H][CX4H][OX1C][P](=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "Glycerol-3-phosphate backbone not found"

    # Check specifically for one ester linkage
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester linkage, found {len(ester_matches)}"

    # Ensure only one chain is bonded to glycerol at the ester position
    # Lysophosphatidic acids will not have an additional fatty acid esterified
    fatty_acid_length = 0
    for bond in mol.GetBonds():
        if bond.GetEndAtom().GetSymbol() == 'C' and bond.GetBeginAtom().GetSymbol() == 'O':
            start_atom = bond.GetEndAtom()
            # Spread out to count carbon chain length
            carbon_seen = set()
            carbon_queue = [start_atom]
            while carbon_queue:
                next_carbon = carbon_queue.pop()
                if next_carbon.GetIdx() not in carbon_seen:
                    carbon_seen.add(next_carbon.GetIdx())
                    carbon_queue.extend(
                        [n for n in next_carbon.GetNeighbors() if n.GetSymbol() == 'C' and n.GetIdx() not in carbon_seen]
                    )
            fatty_acid_length = max(fatty_acid_length, len(carbon_seen))

    if fatty_acid_length < 8 or fatty_acid_length > 24:
        return False, "Fatty acid carbon count is out of typical LPA range"

    return True, "Structure consistent with lysophosphatidic acid: glycerol-3-phosphate backbone with one ester-linked fatty acid"