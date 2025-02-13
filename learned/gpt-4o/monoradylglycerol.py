"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol has a glycerol backbone with one acyl, alkyl, or alk-1-enyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define less restrictive glycerol backbone pattern allowing some variations
    glycerol_pattern = Chem.MolFromSmarts("[O][CH2]CO")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    
    if not glycerol_matches:
        return False, "No glycerol backbone found or pattern too restrictive"

    # Check for the presence of one acyl/alkyl/alk-1-enyl substituent
    # Look for ester linkage pointing to a single substituent
    acyl_pattern = Chem.MolFromSmarts("C(=O)[O]")
    alkyl_pattern = Chem.MolFromSmarts("C[C,C][C,C]")
    alk1enyl_pattern = Chem.MolFromSmarts("C=C")

    # While ensuring correct connection to glycerol
    substituent_count = 0
    reasons = []
    for atom in glycerol_matches[0]:
        if mol.GetAtomWithIdx(atom).GetSymbol() == 'O':  # Focus on substitution-prone oxygen
            acyl_matches = mol.GetAtomWithIdx(atom).GetNeighbors()[0].HasSubstructMatch(acyl_pattern)
            alkyl_matches = mol.GetAtomWithIdx(atom).GetNeighbors()[0].HasSubstructMatch(alkyl_pattern)
            alk1enyl_matches = mol.GetAtomWithIdx(atom).GetNeighbors()[0].HasSubstructMatch(alk1enyl_pattern)

            if acyl_matches:
                substituent_count += 1
                reasons.append("Acyl substituent linked to glycerol")
            if alkyl_matches:
                substituent_count += 1
                reasons.append("Alkyl substituent linked to glycerol")
            if alk1enyl_matches:
                substituent_count += 1
                reasons.append("Alk-1-enyl substituent linked to glycerol")

    if substituent_count == 1:
        return True, f"Contains glycerol backbone with a single substituent: {', '.join(reasons)}"
    
    return False, f"Glycerol backbone has {substituent_count} substitutions; require exactly one"